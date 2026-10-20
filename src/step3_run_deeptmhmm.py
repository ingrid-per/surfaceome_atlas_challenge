"""Step 3b: run the locally installed DeepTMHMM over the proteome, in batches.

Runs the academic package in tools/, unmodified, from a separate conda env
(`deeptmhmm`) because it needs torch 2.2 and numpy<2 while the analysis env is
on current versions. See D26 for why this runs locally rather than hosted.

Three things make this practical:

* **Batching with resume.** Each batch writes its own parsed .tsv. A batch that
  already has one is skipped, so an interrupted run continues instead of
  restarting a multi-hour job.

* **Embeddings are deleted after each batch.** DeepTMHMM caches an ESM
  embedding per residue, about 5 KB each. Over 11.4M residues that is ~58 GB if
  nothing is cleaned up. Per batch it is a few hundred MB.

* **Batches are never size 1.** predict.py draws a per-sequence matplotlib
  figure only when exactly one sequence is submitted, so any batch larger than
  one skips the plotting path. This is why the package needs no patching.

Usage:
    python src/step3_run_deeptmhmm.py --pilot
    python src/step3_run_deeptmhmm.py --batch-size 200
"""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
import time
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"
PKG = ROOT / "tools" / "DeepTMHMM-Academic-License-v1.0"
ENV = "deeptmhmm"

# predict.py writes its own working files; we keep only the parsed result.
WORK = PROC / "step3_deeptmhmm_work"
OUTDIR = PROC / "step3_deeptmhmm"


def read_fasta(path: Path) -> list[tuple[str, str]]:
    recs, acc, chunks = [], None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if acc is not None:
                recs.append((acc, "".join(chunks)))
            acc, chunks = line[1:].strip(), []
        else:
            chunks.append(line.strip())
    if acc is not None:
        recs.append((acc, "".join(chunks)))
    return recs


def parse_gff3(path: Path) -> list[dict]:
    """TMRs.gff3 -> one row per annotated region.

    Region labels are 'outside' (extracellular), 'inside' (cytoplasmic),
    'TMhelix', 'Beta sheet' and 'signal'.
    """
    rows = []
    for line in path.read_text().splitlines():
        if not line or line.startswith("#") or line.startswith("##"):
            continue
        f = line.split("\t")
        if len(f) < 4:
            continue
        try:
            rows.append({"accession": f[0], "region": f[1].strip(),
                         "start": int(f[2]), "end": int(f[3])})
        except ValueError:
            continue
    return rows


def parse_3line(path: Path) -> dict[str, str]:
    """predicted_topologies.3line -> accession -> protein type ('TM', 'GLOB', ...)."""
    types = {}
    lines = path.read_text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith(">"):
            head = line[1:]
            acc, _, ptype = head.partition("|")
            types[acc.strip()] = ptype.strip()
    return types


def summarise(gff_rows: list[dict], types: dict[str, str]) -> list[dict]:
    """Collapse regions into one row per protein, with the D20 length measure."""
    per: dict[str, list[dict]] = {}
    for r in gff_rows:
        per.setdefault(r["accession"], []).append(r)

    out = []
    for acc, regs in per.items():
        outside = [r["end"] - r["start"] + 1 for r in regs if r["region"] == "outside"]
        tm = [r for r in regs if r["region"] in ("TMhelix", "Beta sheet")]
        sig = [r for r in regs if r["region"] == "signal"]
        out.append({
            "accession": acc,
            "dtm_type": types.get(acc, ""),
            "dtm_n_tm": len(tm),
            "dtm_tm_spans": ";".join(f"{r['start']}-{r['end']}" for r in tm),
            "dtm_n_outside": len(outside),
            # D20: longest SINGLE contiguous extracellular segment
            "dtm_longest_outside": max(outside) if outside else 0,
            "dtm_total_outside": sum(outside),
            "dtm_has_signal": bool(sig),
            "dtm_signal_end": max((r["end"] for r in sig), default=0),
            "dtm_surface_call": len(tm) > 0,
        })
    return out


def run_batch(batch: list[tuple[str, str]], idx: int) -> list[dict]:
    WORK.mkdir(parents=True, exist_ok=True)
    fasta = WORK / f"batch_{idx:04d}.fasta"
    with fasta.open("w") as fh:
        for acc, seq in batch:
            fh.write(f">{acc}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    out = WORK / f"out_{idx:04d}"
    if out.exists():
        shutil.rmtree(out)   # predict.py refuses to write into an existing dir

    proc = subprocess.run(
        ["conda", "run", "--no-capture-output", "-n", ENV,
         "python", "predict.py", "--fasta", str(fasta), "--output-dir", str(out)],
        cwd=PKG, capture_output=True, text=True, timeout=60 * 60 * 6,
    )
    gff, three = out / "TMRs.gff3", out / "predicted_topologies.3line"
    if not gff.exists():
        raise RuntimeError(
            f"batch {idx}: no TMRs.gff3 written.\n"
            f"stdout tail: {proc.stdout[-800:]}\nstderr tail: {proc.stderr[-800:]}")

    rows = summarise(parse_gff3(gff),
                     parse_3line(three) if three.exists() else {})

    # Drop the embeddings cache before the next batch; it is the bulk of the
    # disk use and nothing downstream reads it.
    shutil.rmtree(out, ignore_errors=True)
    fasta.unlink(missing_ok=True)
    return rows


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true",
                    help="run the 50-sequence pilot instead of the proteome")
    ap.add_argument("--batch-size", type=int, default=200,
                    help="sequences per batch (must be >1 to skip plotting)")
    ap.add_argument("--limit", type=int, default=0,
                    help="stop after this many sequences (for timing)")
    args = ap.parse_args()

    if args.batch_size < 2:
        print("batch size must be at least 2, or predict.py enters its "
              "per-sequence plotting path")
        return 1
    if not PKG.exists():
        print(f"DeepTMHMM package not found at {PKG}")
        return 1

    src = PROC / ("step3_pilot.fasta" if args.pilot else "step3_proteome.fasta")
    tag = "pilot" if args.pilot else "full"
    records = read_fasta(src)
    if args.limit:
        records = records[:args.limit]

    outdir = OUTDIR / tag
    outdir.mkdir(parents=True, exist_ok=True)

    batches = [records[i:i + args.batch_size]
               for i in range(0, len(records), args.batch_size)]
    print(f"{len(records):,} sequences in {len(batches)} batches of "
          f"{args.batch_size} ({tag})")

    done = t0 = 0
    t0 = time.time()
    for i, batch in enumerate(batches):
        part = outdir / f"batch_{i:04d}.tsv"
        if part.exists():
            done += len(batch)
            continue
        started = time.time()
        rows = run_batch(batch, i)
        pd.DataFrame(rows).to_csv(part, sep="\t", index=False)
        done += len(batch)
        el = time.time() - started
        rate = done / (time.time() - t0)
        remaining = (len(records) - done) / rate if rate else 0
        print(f"  batch {i + 1}/{len(batches)}  {len(rows)} proteins  "
              f"{el:.0f}s  |  {done:,}/{len(records):,}  "
              f"{rate * 60:.1f} seq/min  eta {remaining / 3600:.1f} h", flush=True)

    parts = sorted(outdir.glob("batch_*.tsv"))
    df = pd.concat([pd.read_csv(p, sep="\t") for p in parts], ignore_index=True)
    df = df.drop_duplicates("accession").sort_values("accession", kind="mergesort")
    dest = PROC / f"step3_deeptmhmm_{tag}.csv"
    df.to_csv(dest, index=False)

    missing = {a for a, _ in records} - set(df.accession)
    lines = [
        f"STEP 3b — DeepTMHMM ({tag})",
        "Generated by src/step3_run_deeptmhmm.py — do not edit by hand.",
        "",
        f"sequences submitted                      : {len(records):,}",
        f"proteins with a prediction               : {len(df):,}",
        f"missing (no prediction, must be flagged) : {len(missing):,}",
        "",
        f"called transmembrane                     : {int(df.dtm_surface_call.sum()):,}",
        f"predicted type counts                    : "
        f"{df.dtm_type.value_counts().to_dict()}",
        f"with a predicted signal peptide          : {int(df.dtm_has_signal.sum()):,}",
    ]
    tmdf = df[df.dtm_surface_call]
    if len(tmdf):
        lines += [
            "",
            "longest extracellular segment, TM proteins:",
            f"  median                                 : "
            f"{tmdf.dtm_longest_outside.median():.0f} aa",
            f"  >= 20 aa (D20 default)                 : "
            f"{int((tmdf.dtm_longest_outside >= 20).sum()):,}",
            f"  <  20 aa (Flagged)                     : "
            f"{int((tmdf.dtm_longest_outside < 20).sum()):,}",
        ]
    report = "\n".join(lines)
    print("\n" + report)
    (RESULTS / f"step3b_deeptmhmm_{tag}.txt").write_text(report + "\n")
    print(f"\nWrote {dest}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
