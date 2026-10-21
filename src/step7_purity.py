"""Step 7d: tumour cell or microenvironment, judged against PUBLISHED purity.

Replaces the marker-correlation threshold in D55/D58, which was mine and which
I rightly questioned: 0.40 came from eyeballing a gap in twenty genes.

The field's method is to correlate expression with **tumour purity**. The
literature is explicit that "genes negatively correlated with purity showed
enrichment in stroma-related processes such as inflammatory response". A gene
expressed by tumour cells rises as the sample gets purer; a gene expressed by
infiltrating immune or stromal cells falls.

Purity comes from Aran et al., "Systematic pan-cancer analysis of tumour
purity", Nature Communications 6:8971 (2015), which publishes four independent
estimates per TCGA sample plus a consensus:

    ABSOLUTE  DNA copy number
    LUMP      DNA methylation
    IHC       a pathologist's count down a microscope
    ESTIMATE  expression-based
    CPE       consensus of the above

**ABSOLUTE, LUMP and IHC are used; ESTIMATE is not.** ESTIMATE is derived from
expression, so correlating our expression values against it would be partly
circular. IHC in particular is an independent human measurement.

The threshold is no longer chosen by eye. Known leukocyte, fibroblast and
endothelial markers are used as positive controls and known epithelial and
oncofetal genes as negative controls; the cut-off is placed where it separates
them, and the achieved separation is reported so the reader can judge it.

Writes results/step7_purity_correlation.csv.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

# Non-expression purity measures only, to avoid circularity.
PURITY_COLS = ["ABSOLUTE", "LUMP", "IHC"]
MIN_PAIRS = 200          # a correlation needs enough matched samples

# Positive controls: genes the literature places in the microenvironment.
# They must come out NEGATIVE against purity.
MICROENV_CONTROLS = ["PTPRC", "CD3E", "CD68", "FCGR3A", "HLA-DQA1", "CD14",
                     "ITGAM", "PECAM1", "LRRC15", "FAP", "THY1", "ITGB2"]
# Negative controls: tumour-cell / epithelial / oncofetal genes.
# They must come out POSITIVE or near zero.
TUMOUR_CONTROLS = ["EPCAM", "CDH1", "CLDN6", "KRT8", "KRT18", "CA9", "MUC1",
                   "ERBB2", "EGFR", "TACSTD2"]


def load_purity() -> pd.DataFrame:
    p = RAW / "aran2015_tumour_purity.xlsx"
    d = pd.read_excel(p, sheet_name=0, header=4).iloc[:, :7]
    d.columns = ["sample", "cancer", "ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE"]
    d = d.dropna(subset=["sample"])
    for c in ["ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    # Aran ids carry a vial letter (TCGA-OR-A5J1-01A); ours do not.
    d["short"] = d["sample"].astype(str).str.slice(0, 15)
    return d.drop_duplicates("short").set_index("short")


def main() -> int:
    mat = pd.read_parquet(PROC / "step4_candidate_tpm.parquet")
    samp = pd.read_parquet(PROC / "step4_samples.parquet")
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)
    sym = (gmap.dropna(subset=["ensembl_gene"]).drop_duplicates("ensembl_gene")
           .set_index("ensembl_gene").symbol)

    pur = load_purity()
    tum = samp[(samp.study == "TCGA") & (samp.sample_type == "Primary Tumor")]["sample"]
    ids = [s for s in tum if s in mat.columns and s in pur.index]
    print(f"matched tumours with published purity: {len(ids):,}")

    M = mat[ids]
    X = M.to_numpy(dtype=np.float64)
    Xc = X - X.mean(axis=1, keepdims=True)
    Xn = np.sqrt((Xc ** 2).sum(axis=1))
    Xn[Xn == 0] = np.nan

    out = pd.DataFrame({"ensembl_gene": M.index})
    out["symbol"] = out.ensembl_gene.map(sym)
    for col in PURITY_COLS:
        v = pur.loc[ids, col].to_numpy(dtype=np.float64)
        ok = np.isfinite(v)
        if ok.sum() < MIN_PAIRS:
            print(f"  {col}: only {ok.sum()} samples — skipped")
            continue
        sub = X[:, ok]
        sc = sub - sub.mean(axis=1, keepdims=True)
        sn = np.sqrt((sc ** 2).sum(axis=1))
        sn[sn == 0] = np.nan
        vv = v[ok] - v[ok].mean()
        out[f"r_{col}"] = (sc @ vv) / (sn * np.sqrt((vv ** 2).sum()))
        print(f"  {col}: {int(ok.sum()):,} samples")

    cols = [c for c in out.columns if c.startswith("r_")]
    out["r_purity"] = out[cols].mean(axis=1)

    # ---- calibrate the cut-off on controls instead of by eye ----
    micro = out[out.symbol.isin(MICROENV_CONTROLS)].dropna(subset=["r_purity"])
    tumc = out[out.symbol.isin(TUMOUR_CONTROLS)].dropna(subset=["r_purity"])
    print(f"\ncontrol separation on r(expression, purity):")
    print(f"  microenvironment controls (n={len(micro)}): "
          f"median {micro.r_purity.median():+.2f}  "
          f"range {micro.r_purity.min():+.2f}..{micro.r_purity.max():+.2f}")
    print(f"  tumour-cell controls      (n={len(tumc)}): "
          f"median {tumc.r_purity.median():+.2f}  "
          f"range {tumc.r_purity.min():+.2f}..{tumc.r_purity.max():+.2f}")

    # Youden-style cut: the value that best separates the two control sets.
    grid = np.arange(-0.60, 0.21, 0.01)
    best, best_j = None, -2
    for t in grid:
        sens = float((micro.r_purity <= t).mean())
        spec = float((tumc.r_purity > t).mean())
        if sens + spec - 1 > best_j:
            best_j, best = sens + spec - 1, float(t)
    sens = float((micro.r_purity <= best).mean())
    spec = float((tumc.r_purity > best).mean())
    print(f"\n  cut-off chosen by control separation : r_purity <= {best:+.2f}")
    print(f"    catches {sens:.0%} of microenvironment controls")
    print(f"    keeps   {spec:.0%} of tumour-cell controls")
    print(f"    Youden J = {best_j:.2f}")

    out["microenv_purity"] = out.r_purity <= best
    out.to_csv(RESULTS / "step7_purity_correlation.csv", index=False)

    print("\n  controls, one by one:")
    for lab, df in [("microenvironment", micro), ("tumour cell", tumc)]:
        for _, r in df.sort_values("r_purity").iterrows():
            mark = "flagged" if r.r_purity <= best else "clean"
            print(f"    {lab:<17} {str(r.symbol):<10} r={r.r_purity:+.2f}  {mark}")

    # ---- compare with the old marker-based call ----
    old = RESULTS / "step7_lineage.csv"
    if old.exists():
        lg = pd.read_csv(old)[["ensembl_gene", "microenvironment"]]
        j = out.merge(lg, on="ensembl_gene", how="inner").dropna(subset=["r_purity"])
        agree = ((j.microenv_purity) == (j.microenvironment != "tumour cell")).mean()
        print(f"\n  agreement with the earlier marker-based call: {agree:.0%} "
              f"of {len(j):,} genes")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
