"""Step 7c: is a candidate expressed by the tumour cells, or by what surrounds them?

The runbook lists "bulk tissue mixes tumour, stromal and immune cells" as a
known weakness. This turns it from a caveat into a measurement.

Across 9,185 TCGA primary tumours, correlate every candidate with markers of
the cell types that infiltrate tumours:

    PTPRC  (CD45)   all leukocytes
    CD3E            T cells
    CD68            macrophages
    PECAM1  (CD31)  endothelium

A gene that tracks CD45 across thousands of tumours is reporting immune
infiltrate, not tumour-cell expression. The separation is clean: HLA-DQA1
r=0.76 and FCGR3A r=0.71 against EPCAM r=0.08 and CA9 r=0.10.

This matters most for the pan-cancer view. Immune infiltration is universal
across cancers, so ranking by breadth preferentially surfaces immune markers —
which would be a bad target, since the drug would hit the patient's own
leukocytes.

Reported as a flag, never used to filter: a few real targets are genuinely
expressed on both tumour and immune cells.

Writes results/step7_lineage.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

# Immune markers come from the candidate matrix (they are surface proteins).
MARKERS = {"PTPRC": "leukocyte", "CD3E": "T cell", "CD68": "macrophage",
           "PECAM1": "endothelium"}
# Fibroblast and endothelial markers are NOT surface proteins, so they were
# never in the candidate set and had to be pulled from the raw matrix. Without
# a fibroblast axis the check missed LRRC15, which is a cancer-associated
# fibroblast protein that AbbVie's ABBV-085 targets in the stroma deliberately.
STROMAL = {"COL1A1": "fibroblast", "COL1A2": "fibroblast", "FAP": "fibroblast",
           "PDGFRB": "pericyte/fibroblast", "DCN": "fibroblast",
           "ACTA2": "smooth muscle", "VWF": "endothelium", "ROBO4": "endothelium"}
# Chosen from the observed separation, not from taste: among the 20 broadest
# candidates every known immune marker sits above 0.40 against CD45 and every
# known tumour-cell gene sits below 0.25.
IMMUNE_R = 0.40


def main() -> int:
    mat = pd.read_parquet(PROC / "step4_candidate_tpm.parquet")
    samp = pd.read_parquet(PROC / "step4_samples.parquet")
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)
    sym2g = (gmap.dropna(subset=["ensembl_gene"])
             .drop_duplicates("symbol").set_index("symbol").ensembl_gene)

    tum = samp[(samp.study == "TCGA") & (samp.sample_type == "Primary Tumor")]["sample"]
    ids = [s for s in tum if s in mat.columns]
    M = mat[ids]
    print(f"correlating {len(M):,} candidates across {len(ids):,} TCGA tumours")

    refs = {}
    for sym in MARKERS:
        g = sym2g.get(sym)
        if g in M.index:
            refs[sym] = M.loc[g].to_numpy(dtype=np.float64)
        else:
            print(f"  marker {sym} not in the candidate matrix — skipped")
    if not refs:
        print("no markers available")
        return 1

    # add the stromal axis
    sp = PROC / "step7_stromal_markers.parquet"
    if sp.exists():
        stro = pd.read_parquet(sp)
        stro = stro.loc[[s for s in ids if s in stro.index]]
        for sym in STROMAL:
            if sym in stro.columns and len(stro) == len(ids):
                refs[sym] = stro[sym].to_numpy(dtype=np.float64)
        print(f"  stromal markers added: "
              f"{[s for s in STROMAL if s in refs]}")

    X = M.to_numpy(dtype=np.float64)
    Xc = X - X.mean(axis=1, keepdims=True)
    Xn = np.sqrt((Xc ** 2).sum(axis=1))
    Xn[Xn == 0] = np.nan

    out = pd.DataFrame({"ensembl_gene": M.index})
    for sym, r in refs.items():
        rc = r - r.mean()
        rn = np.sqrt((rc ** 2).sum())
        out[f"r_{sym}"] = (Xc @ rc) / (Xn * rn)

    out["max_immune_r"] = out[[f"r_{s}" for s in ("PTPRC", "CD3E", "CD68")
                               if f"r_{s}" in out.columns]].max(axis=1)
    fib = [f"r_{s}" for s in ("COL1A1", "COL1A2", "FAP", "PDGFRB", "DCN")
           if f"r_{s}" in out.columns]
    out["max_stromal_r"] = out[fib].max(axis=1) if fib else np.nan
    endo = [f"r_{s}" for s in ("PECAM1", "VWF", "ROBO4") if f"r_{s}" in out.columns]
    out["max_endothelial_r"] = out[endo].max(axis=1) if endo else np.nan
    out["likely_infiltrate"] = out.max_immune_r >= IMMUNE_R
    out["likely_stromal"] = out.max_stromal_r >= IMMUNE_R
    out["microenvironment"] = np.select(
        [out.likely_infiltrate & out.likely_stromal,
         out.likely_infiltrate, out.likely_stromal,
         out.max_endothelial_r >= IMMUNE_R],
        ["immune + stromal", "immune infiltrate", "stroma / fibroblast",
         "endothelium"], default="tumour cell")
    out.to_csv(RESULTS / "step7_lineage.csv", index=False)

    full = pd.read_csv(RESULTS / "step7_ranking_full.csv")
    j = full.merge(out, on="ensembl_gene", how="left")
    genes = (j.groupby(["ensembl_gene", "symbol"])
             .agg(n_types=("cancer_type", "nunique"),
                  S_median=("S", "median"),
                  r=("max_immune_r", "first"),
                  flag=("likely_infiltrate", "first")).reset_index())

    print(f"\nmicroenvironment classes: "
          f"{out.microenvironment.value_counts().to_dict()}")
    print(f"\ncandidates flagged as likely immune infiltrate "
          f"(r >= {IMMUNE_R}): {int(genes.flag.sum()):,} of {len(genes):,}")
    broad = genes[genes.n_types >= 8]
    print(f"  among the {len(broad)} genes passing in 8+ cancer types: "
          f"{int(broad.flag.sum())} flagged ({broad.flag.mean():.0%})")
    narrow = genes[genes.n_types <= 2]
    print(f"  among the {len(narrow)} passing in 1-2 cancer types:  "
          f"{int(narrow.flag.sum())} flagged ({narrow.flag.mean():.0%})")
    print("\n  The flag concentrates in broad candidates, which is the point:")
    print("  immune infiltration is universal, so ranking by breadth finds it.")

    print("\ntop 12 by breadth, with the flag:")
    print(f"  {'gene':<11} {'types':>5} {'S med':>6} {'max r':>6}  verdict")
    for _, r in genes.sort_values(["n_types", "S_median"],
                                  ascending=False).head(12).iterrows():
        print(f"  {str(r.symbol):<11} {r.n_types:>5} {r.S_median:>6.2f} "
              f"{r.r:>6.2f}  {'LIKELY INFILTRATE' if r.flag else 'tumour-cell'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
