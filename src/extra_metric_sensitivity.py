"""EXTRA: does the top 10 depend on ranking by S?

Not a runbook step. Step 9 asked whether the answer survives our THRESHOLDS;
this asks whether it survives our choice of METRIC. Motivated by
docs/metric_review.md, which found that no published framework ranks on a
single expression metric.

The same eligible set is ranked 21 ways. Everything except the ordering is held
fixed: eligibility (D52), the 2-per-family limit, and the baseline thresholds.
A gene enters each ranking through whichever cancer type maximises that metric,
exactly as S does.

  M1  S                       tumour median vs the worst healthy tissue (current)
  M2  S, prevalent only       as M1, restricted to genes above a prevalence gate
                              swept at 0.50 / 0.70 / 0.80
  M3  above-healthy fraction  share of TUMOUR SAMPLES above the worst healthy
                              tissue. S compares two point estimates; this uses
                              the whole tumour distribution
  M4  tau                     how concentrated the gene is in healthy tissue
                              (0 = everywhere, 1 = one tissue only)
  M5  composite               rank-rescaled S + prevalence + density + tau,
                              minus breadth, EQUAL weights
  M6  density                 median tumour TPM alone, a deliberate worst case

Each is computed against GTEx only, HPA only, and both (the higher, more
cautious value) wherever the metric depends on the healthy reference.

On M5's weights: ImmunoTar fits feature weights to known drug targets. Doing
that here would be circular, for the same reason D45 refused to calibrate a
cut-off on approved targets. Equal weights are arbitrary and are labelled as
such; the composite is a sensitivity probe, not a proposed ranking.

Writes results/extra_metric_runs.csv, extra_metric_stability.csv,
extra_metric_counts.txt and a figure.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import ACCENT, GREY, apply_style, save_figure
from step7_rank import apply_family_limit

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

TOP_N = 20
REFS = ["gtex", "hpa", "both"]
GATES = [0.50, 0.70, 0.80]


def healthy_max(prof: pd.DataFrame, ref: str) -> pd.Series:
    """Worst healthy tissue under one reference, in log2(TPM+1)."""
    if ref == "both":
        return prof[["gtex_max_log2", "hpa_max_log2"]].max(axis=1)
    return prof[f"{ref}_max_log2"]


def above_healthy_fraction(elig: pd.DataFrame, prof: pd.DataFrame) -> pd.DataFrame:
    """Share of a cancer type's tumour samples above the worst healthy tissue.

    The overlapping coefficient used in the ADC literature needs per-sample
    values on BOTH sides, and HPA publishes one consensus value per tissue, not
    samples. For 55% of eligible rows the worst tissue is an HPA tissue, so a
    true distribution overlap would have to drop those or silently substitute
    GTEx. This uses the full tumour distribution against a scalar healthy value
    instead, which answers the same question and works under every reference.
    """
    mat = pd.read_parquet(PROC / "step4_candidate_tpm.parquet")
    samp = pd.read_parquet(PROC / "step4_samples.parquet")
    tum = samp[(samp.study == "TCGA") & (samp.sample_type == "Primary Tumor")]
    col = {c: i for i, c in enumerate(mat.columns)}
    V = mat.to_numpy(dtype=np.float32)
    gidx = {g: i for i, g in enumerate(mat.index)}
    ids_by_cancer = {c: [col[s] for s in g["sample"] if s in col]
                     for c, g in tum.groupby("disease")}

    thresholds = {r: healthy_max(prof, r) for r in REFS}
    out = {f"above_{r}": np.full(len(elig), np.nan) for r in REFS}
    for k, (gene, cancer) in enumerate(zip(elig.ensembl_gene, elig.cancer_type)):
        ids = ids_by_cancer.get(cancer)
        gi = gidx.get(gene)
        if gi is None or not ids:
            continue
        x = V[gi, ids]
        for r in REFS:
            h = thresholds[r].get(gene, np.nan)
            if np.isfinite(h):
                out[f"above_{r}"][k] = float((x > h).mean())
    for c, v in out.items():
        elig[c] = v
    return elig


def rank_rescale(s: pd.Series) -> pd.Series:
    """Rank-rescale to [0,1]; ties averaged. Rank, not min-max, so one extreme
    value cannot dominate the composite."""
    return s.rank(pct=True, na_option="keep")


def main() -> int:
    apply_style()
    elig = pd.read_csv(RESULTS / "step7_ranking_full.csv")
    prof = pd.read_csv(RESULTS / "step5_healthy_profile.csv", index_col=0)
    top7 = pd.read_csv(RESULTS / "step7_top10.csv").symbol.tolist()

    n_hpa_worst = int((elig.healthy_max_source == "hpa").sum())
    print(f"{len(elig):,} eligible rows, {elig.ensembl_gene.nunique():,} genes")

    # ---- per-reference quantities ----
    for r in REFS:
        elig[f"S_{r}"] = elig.median_tumour_log2tpm1 - elig.ensembl_gene.map(
            healthy_max(prof, r))
    for r in ["gtex", "hpa"]:
        elig[f"tau_{r}"] = elig.ensembl_gene.map(prof[f"{r}_tau"])
    elig = above_healthy_fraction(elig, prof)

    # ---- the 21 rankings ----
    metrics = []
    for r in REFS:
        metrics.append((f"M1 S", r, elig[f"S_{r}"], None))
    for r in REFS:
        for g in GATES:
            metrics.append((f"M2 S, prevalence>={g:.2f}", r, elig[f"S_{r}"],
                            elig.frac_tumours_above_10tpm >= g))
    for r in REFS:
        metrics.append((f"M3 above-healthy fraction", r, elig[f"above_{r}"], None))
    for r in ["gtex", "hpa"]:
        metrics.append((f"M4 tau", r, elig[f"tau_{r}"], None))
    for r in REFS:
        comp = (rank_rescale(elig[f"S_{r}"])
                + rank_rescale(elig.frac_tumours_above_10tpm)
                + rank_rescale(np.log2(elig.tumour_tpm + 1))
                + rank_rescale(elig[f"tau_{r}"] if r != "both"
                               else elig[["tau_gtex", "tau_hpa"]].max(axis=1))
                - rank_rescale(elig.healthy_breadth_above_10tpm))
        metrics.append((f"M5 composite (equal weights)", r, comp, None))
    metrics.append((f"M6 density alone", "n/a", elig.tumour_tpm, None))

    runs, stability = [], {}
    for name, ref, score, mask in metrics:
        d = elig.assign(_m=score)
        if mask is not None:
            d = d[mask]
        d = d[d._m.notna()]
        if d.empty:
            continue
        best = (d.loc[d.groupby("ensembl_gene")._m.idxmax()]
                .sort_values("_m", ascending=False, kind="mergesort"))
        taken = list(apply_family_limit(best, TOP_N).symbol)
        runs.append({"metric": name, "reference": ref,
                     "n_genes": int(d.ensembl_gene.nunique()),
                     "top20": ";".join(map(str, taken))})
        for s in taken:
            stability[s] = stability.get(s, 0) + 1

    rdf = pd.DataFrame(runs)
    rdf.to_csv(RESULTS / "extra_metric_runs.csv", index=False)
    st = (pd.Series(stability, name="rankings_in_top20").rename_axis("symbol")
          .reset_index())
    st["stability"] = st.rankings_in_top20 / len(rdf)
    st["in_step7_top10"] = st.symbol.isin(top7)
    st = st.sort_values(["stability", "symbol"], ascending=[False, True])
    st.to_csv(RESULTS / "extra_metric_stability.csv", index=False)

    base = set(rdf[(rdf.metric == "M1 S") & (rdf.reference == "both")]
               .iloc[0].top20.split(";"))
    smap = dict(zip(st.symbol, st.stability))

    lines = [
        "EXTRA — does the top 10 depend on ranking by S?",
        "Generated by src/extra_metric_sensitivity.py — do not edit by hand.",
        "Not a runbook step. Step 9 varied the thresholds; this varies the metric.",
        "",
        f"rankings computed                        : {len(rdf)}",
        f"eligible rows                            : {len(elig):,}",
        f"  whose worst healthy tissue is HPA      : {n_hpa_worst:,} "
        f"({100*n_hpa_worst/len(elig):.0f}%)",
        "",
        "ranking rules (eligibility, family limit and thresholds held fixed):",
        "  M1 S                      tumour median vs worst healthy tissue",
        "  M2 S, prevalent only      as M1, gated at 0.50 / 0.70 / 0.80",
        "  M3 above-healthy fraction share of tumour SAMPLES above worst healthy",
        "  M4 tau                    concentration of the healthy profile",
        "  M5 composite              equal-weighted, arbitrary by design",
        "  M6 density alone          median tumour TPM, a deliberate worst case",
        "",
        "THE CURRENT TOP 10 under all " + str(len(rdf)) + " rankings:",
        "",
        f"  {'#':<3} {'gene':<10} {'survives':>9}  verdict",
    ]
    for i, s in enumerate(top7, 1):
        v = smap.get(s, 0.0)
        tag = ("metric-independent" if v >= .75 else
               "moderate" if v >= .40 else "DEPENDS ON S")
        lines.append(f"  {i:<3} {s:<10} {v:>8.0%}   {tag}")

    lines += ["", "genes reaching the top 20 under most rankings:", "",
              f"  {'gene':<10} {'survives':>9}  {'in top 10':>10}"]
    for _, r in st.head(22).iterrows():
        lines.append(f"  {str(r.symbol):<10} {r.stability:>8.0%}  "
                     f"{'yes' if r.in_step7_top10 else '':>10}")

    lines += ["", "how far each ranking departs from S (shared out of top 20):", ""]
    for _, r in rdf.iterrows():
        shared = len(base & set(r.top20.split(";")))
        lines.append(f"  {r.metric:<34} {r.reference:<5} {shared:>2}/20")

    # Saturating metrics need a health warning: if the top 20 are separated by
    # almost nothing, their ORDER carries little information even though the
    # ranking is well defined.
    lines += ["", "how well each metric separates its own top 20:", ""]
    for col, label in [("tau_gtex", "M4 tau (GTEx)"), ("tau_hpa", "M4 tau (HPA)"),
                       ("above_both", "M3 above-healthy fraction (both)"),
                       ("S_both", "M1 S (both)")]:
        # drop NA first: some genes have no HPA row at all, and a group that is
        # entirely NA makes idxmax raise rather than skip
        sub = elig[["ensembl_gene", col]].dropna()
        if sub.empty:
            continue
        b = sub.loc[sub.groupby("ensembl_gene")[col].idxmax()][col]
        if b.empty:
            continue
        hi, lo = b.max(), b.nlargest(TOP_N).min()
        lines.append(f"  {label:<34} top20 spans {lo:.4g} to {hi:.4g}"
                     f"  (range {hi-lo:.4g})")
    lines += [
        "",
        "  tau saturates near 1 for any tissue-restricted gene, so its top 20",
        "  are separated by thousandths and their ORDER is close to arbitrary.",
        "  Read M4 as 'is this gene tissue-restricted at all', not as a ranking.",
        "",
        "read the survival percentages with that in mind: a gene can miss the",
        "top 20 of a metric that barely discriminates, and that is weak evidence",
        "against it. M1-M3 discriminate; M4 and M6 are reference points.",
    ]

    report = "\n".join(lines)
    print("\n" + report)
    (RESULTS / "extra_metric_counts.txt").write_text(report + "\n")

    # ---- figure ----
    top = st.head(22).iloc[::-1]
    cols = [ACCENT if f else GREY for f in top.in_step7_top10]
    fig, ax = plt.subplots(figsize=(7.6, 7.0))
    ax.barh(range(len(top)), top.stability * 100, color=cols,
            edgecolor="white", linewidth=.6)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top.symbol, fontsize=9)
    ax.set_xlabel(f"% of {len(rdf)} ranking rules placing the gene in the top {TOP_N}")
    ax.set_title("Does the answer depend on ranking by S?\n"
                 "highlighted = currently in the top 10", loc="left")
    ax.set_xlim(0, 108)
    for i, v in enumerate(top.stability * 100):
        ax.text(v + 1.4, i, f"{v:.0f}%", va="center", fontsize=8)
    fig.tight_layout()
    print(f"\n  wrote {Path(save_figure(fig, 'extra_metric_stability')).relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
