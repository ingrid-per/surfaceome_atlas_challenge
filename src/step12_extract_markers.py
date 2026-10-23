"""Step 12c: pull marker-gene expression out of the raw Toil matrix.

Step 4 kept only the 3,912 surface candidates, and a lineage marker is almost
never a surface protein — TYR, MLANA, SOX10, ALB, KRT19 are all absent from
`step4_candidate_tpm.parquet`. Restricting the composite to genes that happen
to be stored would build the denominator out of surface proteins, the exact
class Step 12 is testing. So this makes its own pass over the 1.26 GB file.

Same rules as Step 4, applied by the same functions rather than copies:
version-stripped Ensembl ids (D17, `strip_version`) and the Toil back-transform
(D11, `to_log2_tpm1`). Output is log2(TPM+1), the scale everything downstream
already uses.

Writes data/processed/step12_marker_tpm.parquet (markers x samples).
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step2_build_gene_map import strip_version          # D17
from step4_extract_matrix import MATRIX, to_log2_tpm1   # D11

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

# Needed for checks 12.1/12.2 whether or not they are markers of a scored
# lineage. TYR and SOX10 are the held-out melanocyte positives; ALB and GFAP
# are negatives for melanocyte (hepatocyte and astrocyte genes).
CONTROL_GENES = ["TYR", "SOX10", "ALB", "GFAP", "MLANA", "PMEL"]


def wanted() -> tuple[set[str], pd.DataFrame]:
    mk = pd.read_csv(RESULTS / "step12_markers_used.csv")
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)
    sym2ens = (gmap.dropna(subset=["ensembl_gene"])
               .drop_duplicates("symbol").set_index("symbol").ensembl_gene)

    ids = set(mk.ensembl_gene.dropna())
    controls = {}
    for s in CONTROL_GENES:
        e = sym2ens.get(s)
        if isinstance(e, str):
            ids.add(e)
            controls[s] = e
        else:
            print(f"  WARNING: control gene {s} has no Ensembl id")
    print(f"marker genes with an Ensembl id : {mk.ensembl_gene.notna().sum()} "
          f"of {len(mk)} rows")
    dropped = mk[mk.ensembl_gene.isna()]
    if len(dropped):
        print(f"  unmapped, cannot be extracted: "
              f"{', '.join(sorted(set(dropped.symbol)))}")
    print(f"control genes                   : {len(controls)}")
    print(f"distinct genes to extract       : {len(ids)}")
    return ids, mk


def main() -> int:
    ids, mk = wanted()
    if not MATRIX.exists():
        print(f"{MATRIX} not found")
        return 1

    kept_rows, kept_ids, header = [], [], None
    n_lines = 0
    with gzip.open(MATRIX, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")[1:]
        for line in fh:
            n_lines += 1
            tab = line.find("\t")
            gene = strip_version(line[:tab])
            if gene in ids:
                vals = np.fromstring(line[tab + 1:], sep="\t", dtype=np.float32)
                kept_rows.append(to_log2_tpm1(vals))
                kept_ids.append(gene)
            if n_lines % 20000 == 0:
                print(f"  {n_lines:,} rows scanned, {len(kept_ids)} kept",
                      flush=True)

    mat = pd.DataFrame(np.vstack(kept_rows), index=kept_ids, columns=header,
                       dtype=np.float32)
    mat = mat[~mat.index.duplicated(keep="first")]
    mat.to_parquet(PROC / "step12_marker_tpm.parquet", compression="zstd")

    found = set(mat.index)
    absent = sorted(ids - found)
    print(f"\nrows in matrix                  : {n_lines:,}")
    print(f"marker/control genes extracted  : {mat.shape[0]} x {mat.shape[1]} samples")
    if absent:
        sym = mk.set_index("ensembl_gene").symbol.to_dict()
        print(f"requested but not in the matrix : {len(absent)} "
              f"({', '.join(sym.get(a, a) for a in absent[:8])})")
    print(f"values are log2(TPM+1); min {np.nanmin(mat.to_numpy()):.3f} "
          f"max {np.nanmax(mat.to_numpy()):.3f}")
    print(f"wrote {(PROC/'step12_marker_tpm.parquet').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
