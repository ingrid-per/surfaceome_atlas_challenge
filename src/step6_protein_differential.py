"""Step 6a: tumour vs adjacent normal at PROTEIN level — the same analysis as
Step 4, run independently on the protein layer.

D48: this is not a confirmation filter. It runs on every measured protein, not
just the surface candidates, so the FDR is honest and the volcano shows the
whole landscape. Whether protein agrees with RNA is measured afterwards and
reported as information for whoever designs the next experiment. Nothing is
dropped for disagreeing.

Statistics deliberately mirror Step 4 so the two layers are comparable:
Wilcoxon rank-sum (D38), BH-FDR within each cohort (D39), minimum 10 samples
per group (D37).

One thing differs and matters (D49): CPTAC's normal is tumour-ADJACENT tissue.
The like-for-like RNA comparison is therefore the `tcga_nat` arm, not `gtex`.

Writes results/step6_protein_tumor_vs_normal.csv and results/step6_counts.txt.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

PROT = "proteomics_gene_abundance_log2_reference_intensity_normalized"
MIN_N = 10          # D37
MIN_OBSERVED = 5    # a protein needs this many non-missing values per group
FDR_MAX = 0.01      # D39
LOG2FC_MIN = 1.0    # see note below

# D12: the 8 cohorts with matched adjacent normals. D49: each maps to a TCGA
# cancer type already analysed at RNA level, so the layers can be joined.
COHORT_TO_TCGA = {
    "CCRCC": "Kidney Clear Cell Carcinoma",
    "COAD": "Colon Adenocarcinoma",
    "HNSCC": "Head & Neck Squamous Cell Carcinoma",
    "LSCC": "Lung Squamous Cell Carcinoma",
    "LUAD": "Lung Adenocarcinoma",
    "OV": "Ovarian Serous Cystadenocarcinoma",
    "PDAC": "Pancreatic Adenocarcinoma",
    "UCEC": "Uterine Corpus Endometrioid Carcinoma",
}


def strip_version(x: str) -> str:
    return str(x).split(".")[0]


def run_cohort(cohort: str) -> tuple[pd.DataFrame, dict]:
    t = pd.read_csv(RAW / f"{cohort}_{PROT}_Tumor.txt", sep="\t", index_col=0)
    n = pd.read_csv(RAW / f"{cohort}_{PROT}_Normal.txt", sep="\t", index_col=0)
    t.index = [strip_version(i) for i in t.index]
    n.index = [strip_version(i) for i in n.index]
    genes = t.index.intersection(n.index)
    a = t.loc[genes].to_numpy(dtype=float)
    b = n.loc[genes].to_numpy(dtype=float)

    info = {"cohort": cohort, "tcga_cancer_type": COHORT_TO_TCGA[cohort],
            "n_tumour": a.shape[1], "n_normal": b.shape[1],
            "genes_in_file": len(genes)}
    if a.shape[1] < MIN_N or b.shape[1] < MIN_N:
        info["status"] = f"SKIPPED: a group is below {MIN_N} (D37)"
        return pd.DataFrame(), info

    # Mass spectrometry leaves ~17% of values missing, and missingness is
    # abundance-dependent (MNAR), so a protein is only tested where enough
    # values were actually observed in BOTH groups.
    obs_a = (~np.isnan(a)).sum(1)
    obs_b = (~np.isnan(b)).sum(1)
    keep = (obs_a >= MIN_OBSERVED) & (obs_b >= MIN_OBSERVED)

    p = np.full(len(genes), np.nan)
    for i in np.flatnonzero(keep):
        x = a[i][~np.isnan(a[i])]
        y = b[i][~np.isnan(b[i])]
        p[i] = mannwhitneyu(x, y, alternative="two-sided",
                            method="asymptotic").pvalue
    ok = np.isfinite(p)
    fdr = np.full(len(genes), np.nan)
    fdr[ok] = multipletests(p[ok], method="fdr_bh")[1]   # D39, within cohort

    lfc = np.nanmedian(a, axis=1) - np.nanmedian(b, axis=1)
    out = pd.DataFrame({
        "ensembl_gene": genes,
        "cohort": cohort,
        "tcga_cancer_type": COHORT_TO_TCGA[cohort],
        "n_tumour": a.shape[1], "n_normal": b.shape[1],
        "n_observed_tumour": obs_a, "n_observed_normal": obs_b,
        "protein_log2fc": lfc, "protein_p": p, "protein_fdr": fdr,
    })
    out = out[ok].copy()
    out["protein_up"] = (out.protein_fdr < FDR_MAX) & (out.protein_log2fc >= LOG2FC_MIN)
    out["protein_down"] = (out.protein_fdr < FDR_MAX) & (out.protein_log2fc <= -LOG2FC_MIN)
    info.update({"status": "ok", "genes_tested": len(out),
                 "up": int(out.protein_up.sum()), "down": int(out.protein_down.sum())})
    return out, info


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    frames, infos = [], []
    for cohort in COHORT_TO_TCGA:
        df, info = run_cohort(cohort)
        infos.append(info)
        if len(df):
            frames.append(df)
        print(f"  {cohort:<6} T={info['n_tumour']:>4} N={info['n_normal']:>4}  "
              f"tested={info.get('genes_tested', 0):>6,}  "
              f"up={info.get('up', 0):>5,} down={info.get('down', 0):>5,}",
              flush=True)

    d = pd.concat(frames, ignore_index=True)

    # annotate with the surface candidate set, but do NOT filter by it (D48)
    cand = pd.read_csv(PROC / "surface_candidates.csv")
    cand = cand[cand.is_candidate]
    link = (cand.dropna(subset=["ensembl_gene"])
            .groupby("ensembl_gene")
            .agg(symbol=("symbol", "first"), surface_tier=("tier", "first")))
    d = d.merge(link, left_on="ensembl_gene", right_index=True, how="left")
    d["is_surface_candidate"] = d.symbol.notna()
    d.to_csv(RESULTS / "step6_protein_tumor_vs_normal.csv", index=False)

    log = pd.DataFrame(infos)
    log.to_csv(RESULTS / "step6_cohorts.csv", index=False)

    sc = d[d.is_surface_candidate]
    lines = [
        "STEP 6a — tumour vs adjacent normal, PROTEIN",
        "Generated by src/step6_protein_differential.py — do not edit by hand.",
        "Run on ALL measured proteins, not only surface candidates (D48), so the",
        "FDR is honest and the volcano shows the whole landscape.",
        f"Wilcoxon (D38); BH-FDR within cohort (D39); FDR < {FDR_MAX}, "
        f"|log2FC| >= {LOG2FC_MIN:g}.",
        "",
        f"cohorts run                              : {len(log[log.status == 'ok'])} of 8",
        f"result rows (gene x cohort)              : {len(d):,}",
        f"distinct proteins tested                 : {d.ensembl_gene.nunique():,}",
        f"  of which surface candidates            : {sc.ensembl_gene.nunique():,}",
        "",
        f"significantly up in tumour               : {int(d.protein_up.sum()):,}",
        f"significantly down in tumour             : {int(d.protein_down.sum()):,}",
        f"  up among surface candidates            : {int(sc.protein_up.sum()):,}",
        f"  down among surface candidates          : {int(sc.protein_down.sum()):,}",
        "",
        f"{'cohort':<8} {'cancer type':<38} {'T':>4} {'N':>4} {'tested':>7} "
        f"{'up':>6} {'down':>6}",
    ]
    for _, r in log.iterrows():
        lines.append(f"  {r.cohort:<6} {r.tcga_cancer_type[:38]:<38} "
                     f"{r.n_tumour:>4} {r.n_normal:>4} "
                     f"{int(r.get('genes_tested', 0) or 0):>7,} "
                     f"{int(r.get('up', 0) or 0):>6,} {int(r.get('down', 0) or 0):>6,}")
    report = "\n".join(lines)
    print("\n" + report)
    (RESULTS / "step6_counts.txt").write_text(report + "\n")
    print(f"\nWrote results/step6_protein_tumor_vs_normal.csv ({len(d):,} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
