"""Step 4a: pull the candidate genes out of the 1.26 GB expression matrix.

One streaming pass. The matrix is 60,498 genes x 19,131 samples; we need 3,912
rows of it, so it is read line by line and only matching rows are kept.

Two things this gets right that are easy to get wrong:

  D11  The stored values are log2(TPM + 0.001), NOT TPM and NOT log2(TPM+1).
       The floor value -9.9658 is exactly log2(0.001). Treating these as TPM
       would corrupt every fold change downstream. The back-transform lives in
       one function, `to_log2_tpm1`, with a test in step4_check.

  D17  Row ids carry a version suffix (ENSG00000146648.22). Matching uses the
       same strip_version() helper as Step 2, so the rule is applied by one
       function everywhere rather than re-implemented per script.

ALL 19,131 sample columns are kept, not just the ones Step 4 compares. Step 5
scores specificity across every vital tissue, so keeping the full panel here
means Step 5 never has to re-read the 1.26 GB file.

Writes data/processed/step4_candidate_tpm.parquet  (genes x samples, float32)
   and data/processed/step4_samples.parquet        (sample metadata)
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step2_build_gene_map import strip_version   # D17: one shared helper

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

MATRIX = RAW / "TcgaTargetGtex_rsem_gene_tpm.gz"
PHENO = RAW / "TcgaTargetGTEX_phenotype.txt.gz"

# log2(0.001); the matrix's floor for "not detected".
XENA_FLOOR = -9.9658
XENA_OFFSET = 0.001


def to_log2_tpm1(values: np.ndarray) -> np.ndarray:
    """log2(TPM + 0.001)  ->  log2(TPM + 1).   (D11)

    TPM = 2^x - 0.001, clipped at 0 because floating point can leave the floor
    value a hair below zero, and a negative TPM is meaningless.
    """
    tpm = np.power(2.0, values) - XENA_OFFSET
    np.clip(tpm, 0.0, None, out=tpm)
    return np.log2(tpm + 1.0)


def wanted_genes() -> tuple[set[str], pd.DataFrame]:
    """Version-stripped Ensembl genes for every surface candidate."""
    cand = pd.read_csv(PROC / "surface_candidates.csv")
    cand = cand[cand.is_candidate]
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)
    link = gmap[gmap.accession.isin(set(cand.accession))][
        ["accession", "ensembl_gene", "symbol"]].dropna(subset=["ensembl_gene"])
    link = link.drop_duplicates()
    return set(link.ensembl_gene), link


def load_samples() -> pd.DataFrame:
    rows = []
    with gzip.open(PHENO, "rt", encoding="latin-1") as fh:
        header = next(fh).rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if not f or not f[0]:
                continue

            def get(col: str) -> str:
                i = idx.get(col)
                return f[i].strip() if i is not None and i < len(f) else ""

            rows.append({
                "sample": f[0],
                "study": get("_study"),
                "sample_type": get("_sample_type"),
                "disease": get("primary disease or tissue"),
                "primary_site": get("_primary_site"),
            })
    return pd.DataFrame(rows)


def main() -> int:
    PROC.mkdir(parents=True, exist_ok=True)
    RESULTS.mkdir(exist_ok=True)

    want, link = wanted_genes()
    print(f"candidate Ensembl genes sought : {len(want):,}")

    with gzip.open(MATRIX, "rt") as fh:
        header = next(fh).rstrip("\n").split("\t")
        samples = header[1:]
        kept_ids, kept_rows = [], []
        n_in = 0
        for line in fh:
            tab = line.find("\t")
            gene = strip_version(line[:tab])
            n_in += 1
            if gene in want:
                vals = np.fromstring(line[tab + 1:], sep="\t", dtype=np.float64)
                if vals.size != len(samples):
                    raise RuntimeError(
                        f"{gene}: {vals.size} values for {len(samples)} samples")
                kept_ids.append(gene)
                kept_rows.append(to_log2_tpm1(vals).astype(np.float32))
            if n_in % 20000 == 0:
                print(f"  ... {n_in:,} rows read, {len(kept_ids):,} kept", flush=True)

    mat = pd.DataFrame(np.vstack(kept_rows), index=kept_ids, columns=samples)
    mat.index.name = "ensembl_gene"
    # A gene can appear more than once if the matrix carries two versioned ids
    # for the same stripped gene; keep the first and report it.
    dupes = mat.index.duplicated()
    if dupes.any():
        print(f"  {int(dupes.sum())} duplicate gene rows collapsed (first kept)")
        mat = mat[~dupes]

    mat.to_parquet(PROC / "step4_candidate_tpm.parquet", compression="zstd")

    samp = load_samples()
    samp = samp[samp["sample"].isin(samples)]
    samp.to_parquet(PROC / "step4_samples.parquet", index=False)

    missing = sorted(want - set(mat.index))
    miss_tbl = link[link.ensembl_gene.isin(missing)].drop_duplicates("ensembl_gene")
    miss_tbl.to_csv(RESULTS / "step4_genes_not_in_matrix.csv", index=False)

    vals = mat.to_numpy()
    lines = [
        "STEP 4a — candidate expression matrix",
        "Generated by src/step4_extract_matrix.py — do not edit by hand.",
        "",
        f"matrix rows read                         : {n_in:,}",
        f"sample columns                           : {len(samples):,}",
        f"candidate genes sought                   : {len(want):,}",
        f"candidate genes extracted                : {len(mat):,}",
        f"candidate genes absent from the matrix   : {len(missing):,}  "
        f"(listed in results/step4_genes_not_in_matrix.csv; D19)",
        "",
        "values after the D11 back-transform, log2(TPM+1):",
        f"  min / median / max                     : {vals.min():.3f} / "
        f"{np.median(vals):.3f} / {vals.max():.3f}",
        f"  exactly zero (not detected)            : {(vals == 0).mean():.1%}",
        f"  negative values (must be 0)            : {int((vals < 0).sum())}",
        "",
        f"parquet size                             : "
        f"{(PROC / 'step4_candidate_tpm.parquet').stat().st_size / 1e6:.0f} MB",
    ]
    report = "\n".join(lines)
    print("\n" + report)
    (RESULTS / "step4a_extract.txt").write_text(report + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
