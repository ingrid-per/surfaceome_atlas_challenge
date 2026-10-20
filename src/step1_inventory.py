"""Step 1 inventory: open every downloaded file and report what is actually in it.

The point is to see the shape and the columns rather than assume them, and to
print sample counts per group so an impossible number is visible immediately.

The 1.26 GB expression matrix is never loaded into memory — it is streamed, so
this runs in a few minutes and a few megabytes.

Writes results/step1_inventory.csv and results/step1_counts.txt.
"""

from __future__ import annotations

import csv
import gzip
import io
import sys
import zipfile
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step1_resources import CPTAC_COHORTS, CPTAC_NO_NORMAL

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
RESULTS = ROOT / "results"

rows: list[dict] = []
report: list[str] = []


def out(line: str = "") -> None:
    print(line)
    report.append(line)


def record(path: Path, nrows, ncols, note: str) -> None:
    rows.append({"file": path.name, "rows": nrows, "columns": ncols,
                 "bytes": path.stat().st_size, "note": note})


def peek_zip_tsv(path: Path, nshow: int = 2) -> tuple[int, int, list[str]]:
    """Row count, column count and header of the single TSV inside a zip."""
    with zipfile.ZipFile(path) as z:
        name = z.namelist()[0]
        with z.open(name) as fh:
            text = io.TextIOWrapper(fh, encoding="utf8", errors="replace")
            header = next(text).rstrip("\n").split("\t")
            n = 0
            shown = []
            for line in text:
                if n < nshow:
                    shown.append(line.rstrip("\n").split("\t"))
                n += 1
    out(f"  inner file : {name}")
    out(f"  shape      : {n:,} rows x {len(header)} columns")
    out(f"  columns    : {header[:12]}{' ...' if len(header) > 12 else ''}")
    for s in shown:
        out(f"  sample row : {s[:6]}")
    return n, len(header), header


def peek_text(path: Path, opener=open, nshow: int = 2) -> tuple[int, int]:
    with opener(path, "rt", encoding="utf8", errors="replace") as fh:
        header = next(fh).rstrip("\n").split("\t")
        n, shown = 0, []
        for line in fh:
            if n < nshow:
                shown.append(line.rstrip("\n").split("\t"))
            n += 1
    out(f"  shape      : {n:,} rows x {len(header)} columns")
    out(f"  columns    : {header[:10]}{' ...' if len(header) > 10 else ''}")
    for s in shown:
        out(f"  sample row : {s[:5]}")
    return n, len(header)


# ---------------------------------------------------------------- UniProt ----
def inventory_uniprot() -> None:
    out("\n" + "=" * 78)
    out("UNIPROT")
    out("=" * 78)

    p = RAW / "uniprot_human_reviewed.fasta.gz"
    n = 0
    with gzip.open(p, "rt", encoding="utf8") as fh:
        first = ""
        for line in fh:
            if line.startswith(">"):
                n += 1
                if n == 1:
                    first = line.strip()
    out(f"\n{p.name}")
    out(f"  sequences  : {n:,}")
    out(f"  first entry: {first[:90]}")
    record(p, n, 1, "reviewed human protein sequences")

    p = RAW / "uniprot_human_reviewed_topology.tsv.gz"
    out(f"\n{p.name}")
    nrows, ncols = peek_text(p, gzip.open)
    # How many entries actually carry the annotations Step 3 depends on?
    counts = Counter()
    with gzip.open(p, "rt", encoding="utf8") as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            for col in ("Transmembrane", "Topological domain", "Signal peptide",
                        "Lipidation"):
                i = idx.get(col)
                if i is not None and i < len(f) and f[i].strip():
                    counts[col] += 1
    out(f"  entries with annotation: {dict(counts)}")
    record(p, nrows, ncols, f"topology annotations: {dict(counts)}")


# -------------------------------------------------------------------- HPA ----
def inventory_hpa() -> None:
    out("\n" + "=" * 78)
    out("HUMAN PROTEIN ATLAS 25.1")
    out("=" * 78)
    for name, note in [
        ("proteinatlas.tsv.zip", "gene table (Step 2 spine)"),
        ("rna_tissue_hpa.tsv.zip", "normal-tissue nTPM (Step 5 HPA arm)"),
        ("normal_ihc_data.tsv.zip", "antibody staining in normal tissue (Step 6)"),
        ("subcellular_location.tsv.zip", "subcellular location calls (Step 3)"),
    ]:
        p = RAW / name
        out(f"\n{name}")
        nrows, ncols, header = peek_zip_tsv(p)
        record(p, nrows, ncols, note)
        # For the tissue table, show how many distinct tissues are covered.
        if "rna_tissue" in name or "normal_ihc" in name:
            tissues = set()
            with zipfile.ZipFile(p) as z, z.open(z.namelist()[0]) as fh:
                t = io.TextIOWrapper(fh, encoding="utf8", errors="replace")
                hdr = next(t).rstrip("\n").split("\t")
                ti = hdr.index("Tissue") if "Tissue" in hdr else None
                if ti is not None:
                    for line in t:
                        f = line.rstrip("\n").split("\t")
                        if ti < len(f):
                            tissues.add(f[ti])
            if tissues:
                out(f"  tissues    : {len(tissues)} distinct")


# ------------------------------------------------------------- Xena Toil ----
def inventory_xena() -> None:
    out("\n" + "=" * 78)
    out("UCSC XENA TOIL (TCGA + TARGET + GTEx)")
    out("=" * 78)

    pheno = RAW / "TcgaTargetGTEX_phenotype.txt.gz"
    out(f"\n{pheno.name}")
    meta: dict[str, dict] = {}
    with gzip.open(pheno, "rt", encoding="latin-1") as fh:
        hdr = next(fh).rstrip("\n").split("\t")
        i = {c: j for j, c in enumerate(hdr)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f or not f[0]:
                continue
            meta[f[0]] = {
                "study": f[i["_study"]] if i.get("_study", 99) < len(f) else "",
                "type": f[i["_sample_type"]] if i.get("_sample_type", 99) < len(f) else "",
                "disease": f[i["primary disease or tissue"]]
                if i.get("primary disease or tissue", 99) < len(f) else "",
            }
    out(f"  samples described: {len(meta):,}")
    record(pheno, len(meta), len(hdr), "per-sample study / type / disease")

    probemap = RAW / "gencode.v23.annotation.gene.probemap"
    out(f"\n{probemap.name}")
    nrows, ncols = peek_text(probemap)
    record(probemap, nrows, ncols, "Ensembl gene id -> symbol")

    matrix = RAW / "TcgaTargetGtex_rsem_gene_tpm.gz"
    out(f"\n{matrix.name}  (streamed, not loaded)")
    with gzip.open(matrix, "rt", encoding="utf8") as fh:
        header = next(fh).rstrip("\n").split("\t")
        samples = header[1:]
        ngenes = sum(1 for _ in fh)
    out(f"  shape      : {ngenes:,} genes x {len(samples):,} samples")
    record(matrix, ngenes, len(samples) + 1,
           "log2(TPM+0.001) — back-transform per D11")

    out("\n  samples by study and sample type:")
    combo = Counter((meta.get(s, {}).get("study", "?"),
                     meta.get(s, {}).get("type", "?")) for s in samples)
    for (study, stype), n in sorted(combo.items(), key=lambda kv: -kv[1]):
        out(f"    {n:>6,}  {study:<8} {stype}")

    # The groups D10 actually keeps.
    keep_tcga = {"Primary Tumor", "Solid Tissue Normal"}
    tcga_t = sum(n for (s, t), n in combo.items() if s == "TCGA" and t == "Primary Tumor")
    tcga_n = sum(n for (s, t), n in combo.items() if s == "TCGA" and t == "Solid Tissue Normal")
    gtex_n = sum(n for (s, t), n in combo.items() if s == "GTEX" and t == "Normal Tissue")
    out(f"\n  included by D10: TCGA tumor {tcga_t:,} | TCGA normal {tcga_n:,} | "
        f"GTEx normal {gtex_n:,}  = {tcga_t + tcga_n + gtex_n:,}")
    out(f"  excluded by D10: {len(samples) - tcga_t - tcga_n - gtex_n:,}")

    out("\n  TCGA primary tumors per cancer type:")
    per = Counter(meta[s]["disease"] for s in samples
                  if meta.get(s, {}).get("study") == "TCGA"
                  and meta[s]["type"] == "Primary Tumor")
    for d, n in sorted(per.items(), key=lambda kv: -kv[1]):
        out(f"    {n:>5,}  {d}")
    out(f"  -> {len(per)} cancer types")

    out("\n  TCGA solid tissue normals per cancer type (thin ones matter for Step 4):")
    pern = Counter(meta[s]["disease"] for s in samples
                   if meta.get(s, {}).get("study") == "TCGA"
                   and meta[s]["type"] == "Solid Tissue Normal")
    for d, n in sorted(pern.items(), key=lambda kv: -kv[1]):
        out(f"    {n:>5,}  {d}")
    out(f"  -> {len(pern)} cancer types have any matched normal")

    out("\n  GTEx normal samples per tissue:")
    perg = Counter(meta[s]["disease"] for s in samples
                   if meta.get(s, {}).get("study") == "GTEX"
                   and meta[s]["type"] == "Normal Tissue")
    for d, n in sorted(perg.items(), key=lambda kv: -kv[1]):
        out(f"    {n:>5,}  {d}")
    out(f"  -> {len(perg)} tissues")

    return per, pern, perg


# ------------------------------------------------------------------ CPTAC ----
def inventory_cptac() -> None:
    out("\n" + "=" * 78)
    out("CPTAC (LinkedOmicsKB freeze v1.2_reorganized)")
    out("=" * 78)
    prot = "proteomics_gene_abundance_log2_reference_intensity_normalized"
    out(f"\n  {'cohort':8} {'genes':>8} {'tumors':>8} {'normals':>8}")
    for c in CPTAC_COHORTS:
        counts = {}
        genes = None
        for kind in ("Tumor", "Normal"):
            p = RAW / f"{c}_{prot}_{kind}.txt"
            if not p.exists():
                counts[kind] = None
                continue
            with p.open(encoding="utf8", errors="replace") as fh:
                hdr = next(fh).rstrip("\n").split("\t")
                n = sum(1 for _ in fh)
            counts[kind] = len(hdr) - 1
            genes = n
            record(p, n, len(hdr), f"{c} {kind.lower()} protein abundance")
        note = "  <- no matched normal (D12)" if c in CPTAC_NO_NORMAL else ""
        out(f"  {c:8} {genes or 0:>8,} {counts.get('Tumor') or 0:>8,} "
            f"{counts.get('Normal') or 0:>8,}{note}")


# ------------------------------------------- GTEx proteome / CSPA / SURFY ----
def inventory_spreadsheets() -> None:
    out("\n" + "=" * 78)
    out("SPREADSHEET RESOURCES")
    out("=" * 78)
    for name, note in [
        ("Jiang2020_mmc3.xlsx", "GTEx proteome, Cell 2020 Table S3"),
        ("CSPA_S2_File.xlsx", "CSPA validated surfaceome"),
        ("SURFY_table_S3_surfaceome.xlsx", "SURFY predicted surfaceome"),
    ]:
        p = RAW / name
        out(f"\n{name}")
        with zipfile.ZipFile(p) as z:
            import re
            wb = z.read("xl/workbook.xml").decode("utf8", "replace")
            sheets = re.findall(r'<sheet[^>]*?name="([^"]+)"', wb)
        out(f"  sheets ({len(sheets)}):")
        for s in sheets:
            out(f"    - {s}")
        record(p, len(sheets), 0, f"{note}; sheets: {'; '.join(sheets)}")


def main() -> int:
    if not RAW.exists() or not any(RAW.iterdir()):
        print("data/raw is empty — run step1_download.py --download first")
        return 1
    RESULTS.mkdir(exist_ok=True)

    out("STEP 1 INVENTORY")
    out("Generated by src/step1_inventory.py — do not edit by hand.")

    inventory_uniprot()
    inventory_hpa()
    inventory_xena()
    inventory_cptac()
    inventory_spreadsheets()

    with (RESULTS / "step1_inventory.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["file", "rows", "columns", "bytes", "note"])
        w.writeheader()
        w.writerows(rows)
    (RESULTS / "step1_counts.txt").write_text("\n".join(report) + "\n")

    out("\n" + "=" * 78)
    out(f"{len(rows)} files inventoried. "
        f"Wrote results/step1_inventory.csv and results/step1_counts.txt")
    (RESULTS / "step1_counts.txt").write_text("\n".join(report) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
