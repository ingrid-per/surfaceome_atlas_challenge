"""Step 9: how much does the top 10 depend on the thresholds we chose?

Every cut-off in this pipeline is a convention. The current top 10 rests on one
setting of five of them. This re-runs the selection across a grid of reasonable
settings and reports, for each gene, the share of runs in which it reaches the
top 20 — its STABILITY. A gene in nearly every run is a robust finding; a gene
in a handful depends on one choice and must not lead a headline list.

Nothing is re-predicted and nothing is re-tested:

  - D20 stored the MEASURED longest extracellular segment rather than only the
    pass/fail verdict, so re-tiering at another cut-off is a re-read. The
    re-tiering is done by running src/step3_build_candidates.py itself, not by
    a copy of its rule — check 9.2 confirms the 20 aa re-run reproduces the
    committed candidate set.
  - BH-adjusted q-values do not depend on the threshold applied to them, and
    re-tiering does not change WHICH genes were tested, so the Step 4 q-values
    stay valid at every cell of the grid.
  - Step 5 stored each healthy summary per source, so every healthy setting is
    a column lookup.

Grid (the runbook's table, adapted where later decisions changed it):

    minimum extracellular length   10, 20, 50 residues          (D20)
    FDR threshold                  0.05, 0.01                   (D39)
    log2 fold change               1, 2, 3                      (D39)
    healthy reference              GTEx only, HPA only, both     (D2/D9)
    healthy summary                max; max excluding brain; max excluding
                                   testis+placenta; max excluding all three;
                                   90th percentile              (D44/D46)

The runbook's "vital-organ list, with and without testis / brain" row is
covered by the healthy-summary axis: D42 removed the vital-organ list, and D46
replaced it with exactly those exclusion variants.

Held fixed, because they are not conventions in the same sense: the 5 TPM
tumour floor (D43 — below it there is nothing for a drug to bind), the 2-per-
family limit (D52) and ranking on S (D52).

Writes results/step9_runs.csv, step9_stability.csv, step9_counts.txt and
results/figures/step9_stability.*
"""

from __future__ import annotations

import itertools
import subprocess
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import ACCENT, GREY, apply_style, save_figure
from step7_rank import MAX_PER_FAMILY, apply_family_limit, gene_groups

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

MIN_EXTRA = [10, 20, 50]
FDR = [0.05, 0.01]
LOG2FC = [1.0, 2.0, 3.0]
SOURCE = ["gtex", "hpa", "both"]
SUMMARY = ["max", "max_nobrain", "max_no_reproductive_fetal",
           "max_relaxed", "p90"]

TOP_N = 20              # the runbook asks for the share of runs in the top 20
MIN_TUMOUR_TPM = 5.0    # D43, held fixed
ELIGIBLE_TIERS = {"S1", "S2", "S3"}
TIER_RANK = {"S1": 0, "S2": 1, "S3": 2, "Flagged": 3}

# The default cell. Check 9.1: its top 10 must equal Step 7's exactly.
DEFAULT = dict(min_extra=20, fdr=0.01, log2fc=2.0, source="both", summary="max")


def retier(min_extra: int) -> Path:
    """Run Step 3's own tiering at this cut-off and return the output path."""
    out = PROC / f"surface_candidates_me{min_extra}.csv"
    cmd = [sys.executable, str(ROOT / "src" / "step3_build_candidates.py"),
           "--min-extracellular", str(min_extra), "--tag", f"me{min_extra}"]
    subprocess.run(cmd, check=True, capture_output=True)
    return out


def gene_tier(cand_path: Path, gmap: pd.DataFrame) -> pd.Series:
    """Ensembl gene -> its strongest surface tier, as Step 4 collapses it."""
    cand = pd.read_csv(cand_path)
    cand = cand[cand.is_candidate]
    link = gmap[gmap.accession.isin(set(cand.accession))][
        ["accession", "ensembl_gene"]].dropna().drop_duplicates()
    j = link.merge(cand[["accession", "tier"]], on="accession")
    j["rank"] = j.tier.map(TIER_RANK)
    j = j.sort_values(["ensembl_gene", "rank"], kind="mergesort")
    return j.groupby("ensembl_gene").tier.first()


def healthy_level(prof: pd.DataFrame, source: str, summary: str) -> pd.Series:
    """The healthy expression level under one (source, summary) setting."""
    stem = "p90_log2" if summary == "p90" else f"{summary}_log2"
    if source == "both":
        # D2/D9: two independent pipelines, take the more cautious (higher).
        return prof[[f"gtex_{stem}", f"hpa_{stem}"]].max(axis=1)
    return prof[f"{source}_{stem}"]


def main() -> int:
    apply_style()
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)
    s5 = pd.read_csv(RESULTS / "step5_specificity.csv")
    prof = pd.read_csv(RESULTS / "step5_healthy_profile.csv", index_col=0)
    grp, _ = gene_groups()

    # ---- re-tier at each cut-off, through Step 3's own code ----
    tiers = {}
    for me in MIN_EXTRA:
        tiers[me] = gene_tier(retier(me), gmap)
        n = int(tiers[me].isin(ELIGIBLE_TIERS).sum())
        print(f"  {me:>2} aa : {n:,} genes in S1/S2/S3")

    # check 9.2 — the 20 aa re-run must reproduce the committed candidate set
    committed = gene_tier(PROC / "surface_candidates.csv", gmap)
    rerun_ok = committed.equals(tiers[DEFAULT["min_extra"]].reindex(committed.index))
    stored = s5.drop_duplicates("ensembl_gene").set_index("ensembl_gene").tier
    collapse_ok = bool(
        (stored == committed.reindex(stored.index)).all())

    # ---- the pieces that do not change across the grid ----
    s5["_eligible_tpm"] = s5.tumour_tpm >= MIN_TUMOUR_TPM
    s5["gene_group"] = s5.symbol.map(grp)
    base = s5[["ensembl_gene", "symbol", "gene_group", "cancer_type",
               "median_tumour_log2tpm1", "fdr", "log2fc", "_eligible_tpm"]].copy()

    runs, stability = [], {}
    grid = list(itertools.product(MIN_EXTRA, FDR, LOG2FC, SOURCE, SUMMARY))
    print(f"\n{len(grid)} settings to run over "
          f"{base.ensembl_gene.nunique():,} genes x "
          f"{base.cancer_type.nunique()} cancer types")

    # cache the healthy level per (source, summary): 15 lookups, not 180
    hcache = {(src, summ): healthy_level(prof, src, summ)
              for src in SOURCE for summ in SUMMARY}

    for k, (me, fdr, lfc, src, summ) in enumerate(grid, 1):
        d = base.copy()
        d["tier"] = d.ensembl_gene.map(tiers[me])
        d["S"] = d.median_tumour_log2tpm1 - d.ensembl_gene.map(hcache[(src, summ)])

        elig = d[d.tier.isin(ELIGIBLE_TIERS) & d._eligible_tpm
                 & (d.fdr < fdr) & (d.log2fc >= lfc) & d.S.notna()]
        if elig.empty:
            continue
        # each gene's best cancer type, ranked by S, family limit applied (D52)
        best = (elig.loc[elig.groupby("ensembl_gene").S.idxmax()]
                .sort_values("S", ascending=False, kind="mergesort"))
        picked = apply_family_limit(best, TOP_N)
        taken = list(picked.symbol)

        is_default = (me == DEFAULT["min_extra"] and fdr == DEFAULT["fdr"]
                      and lfc == DEFAULT["log2fc"] and src == DEFAULT["source"]
                      and summ == DEFAULT["summary"])
        runs.append({"run": k, "min_extra": me, "fdr": fdr, "log2fc": lfc,
                     "source": src, "summary": summ, "is_default": is_default,
                     "n_genes_eligible": int(elig.ensembl_gene.nunique()),
                     "top20": ";".join(map(str, taken))})
        for s in taken:
            stability[s] = stability.get(s, 0) + 1
        if k % 30 == 0:
            print(f"  {k}/{len(grid)}", flush=True)

    rdf = pd.DataFrame(runs)
    rdf.to_csv(RESULTS / "step9_runs.csv", index=False)

    top7 = pd.read_csv(RESULTS / "step7_top10.csv").symbol.tolist()
    st = (pd.Series(stability, name="runs_in_top20").rename_axis("symbol")
          .reset_index())
    st["stability"] = st.runs_in_top20 / len(rdf)
    st["in_step7_top10"] = st.symbol.isin(top7)
    st = st.sort_values(["stability", "symbol"], ascending=[False, True])
    st.to_csv(RESULTS / "step9_stability.csv", index=False)

    # ---- check 9.1: the default cell reproduces Step 7 ----
    drow = rdf[rdf.is_default]
    default_top20 = drow.iloc[0].top20.split(";") if len(drow) else []
    default_top10 = default_top20[:10]
    reproduces = default_top10 == top7

    smap = dict(zip(st.symbol, st.stability))
    lines = [
        "STEP 9 — robustness of the top 10 to the thresholds we chose",
        "Generated by src/step9_robustness.py — do not edit by hand.",
        "",
        "Nothing is re-predicted. Re-tiering runs step3_build_candidates.py at",
        "each cut-off (the payoff of D20 storing the measured length); the BH",
        "q-values are threshold-independent and the tested gene set does not",
        "change, so Step 4 is not re-run.",
        "",
        f"settings run                             : {len(rdf)} of {len(grid)}",
        f"  minimum extracellular (D20)            : {MIN_EXTRA}",
        f"  FDR (D39)                              : {FDR}",
        f"  log2 fold change (D39)                 : {LOG2FC}",
        f"  healthy reference (D2/D9)              : {SOURCE}",
        f"  healthy summary (D44/D46)              : {SUMMARY}",
        f"  held fixed                             : 5 TPM tumour floor (D43),",
        f"                                           {MAX_PER_FAMILY} per family (D52),",
        "                                           ranking on S (D52)",
        "",
        "checks:",
        f"  9.1 default cell reproduces the Step 7 top 10 : {reproduces}",
        f"  9.2 20 aa re-tier reproduces the committed set: {rerun_ok}",
        f"  9.3 gene-level tier collapse matches Step 5    : {collapse_ok}",
        f"  9.4 stability values within [0, 1]             : "
        f"{bool(st.stability.between(0, 1).all())}",
        "",
        "genes eligible per setting:",
        f"  min {rdf.n_genes_eligible.min():,}   "
        f"median {int(rdf.n_genes_eligible.median()):,}   "
        f"max {rdf.n_genes_eligible.max():,}",
        "",
        f"genes reaching the top {TOP_N} in at least one run : {len(st):,}",
        f"  in every run (stability 1.00)                 : "
        f"{int((st.stability == 1).sum())}",
        f"  in >= 75% of runs                             : "
        f"{int((st.stability >= .75).sum())}",
        f"  in >= 50% of runs                             : "
        f"{int((st.stability >= .50).sum())}",
        f"  in <  10% of runs (threshold-dependent)       : "
        f"{int((st.stability < .10).sum())}",
        "",
        "THE CURRENT TOP 10, by how well each survives the grid:",
        "",
        f"  {'#':<3} {'gene':<10} {'stability':>9}  verdict",
    ]
    for i, s in enumerate(top7, 1):
        v = smap.get(s, 0.0)
        tag = ("robust" if v >= .75 else "moderate" if v >= .40
               else "threshold-dependent")
        lines.append(f"  {i:<3} {s:<10} {v:>8.0%}   {tag}")

    lines += ["", f"THE ROBUST CORE — top {TOP_N} under most settings:", "",
              f"  {'gene':<10} {'stability':>9}  {'in Step 7 top 10':>17}"]
    for _, r in st.head(25).iterrows():
        lines.append(f"  {str(r.symbol):<10} {r.stability:>8.0%}  "
                     f"{'yes' if r.in_step7_top10 else '':>17}")

    # which single choice moves the answer most
    dtop = set(default_top20)
    lines += ["", "which choice moves the answer most",
              f"(mean overlap of a run's top {TOP_N} with the default run's):", ""]
    spread, sens_rows = [], []
    for axis in ["min_extra", "fdr", "log2fc", "source", "summary"]:
        per = {}
        for val, g in rdf.groupby(axis):
            per[val] = float(np.mean([len(dtop & set(x.split(";"))) / len(dtop)
                                      for x in g.top20]))
        spread.append((max(per.values()) - min(per.values()), axis, per))
        rng = max(per.values()) - min(per.values())
        for val, ov in per.items():
            sens_rows.append({"axis": axis, "value": val,
                              "mean_overlap_with_default": round(ov, 4),
                              "axis_spread": round(rng, 4)})
    # Written out so the deliverable reads these rather than restating them.
    pd.DataFrame(sens_rows).sort_values(
        ["axis_spread", "axis", "mean_overlap_with_default"],
        ascending=[False, True, False]).to_csv(
        RESULTS / "step9_sensitivity.csv", index=False)
    for rng, axis, per in sorted(spread, reverse=True):
        lines.append(f"  {axis}  (spread {rng:.0%})")
        for val, ov in sorted(per.items(), key=lambda kv: -kv[1]):
            lines.append(f"      {str(val):<26} overlap {ov:>5.0%}")
        lines.append("")

    n_hpa_missing = int(prof.hpa_max_log2.isna().sum())
    lines += [
        "read this with two caveats:",
        f"  - the HPA-only runs score {n_hpa_missing:,} genes with no HPA row as",
        "    missing rather than clean, so they drop out of those runs instead",
        "    of ranking highly; that is the cautious direction.",
        "  - stability measures sensitivity to OUR thresholds only. It says",
        "    nothing about whether the underlying measurement is right, and a",
        "    gene can be perfectly stable and still a poor target (OR2I1 is",
        "    stable and is a putative pseudogene product).",
    ]

    report = "\n".join(lines)
    print("\n" + report)
    (RESULTS / "step9_counts.txt").write_text(report + "\n")

    # ---------------------------------------------------------- figure --
    top = st.head(25).iloc[::-1]
    cols = [ACCENT if f else GREY for f in top.in_step7_top10]
    fig, ax = plt.subplots(figsize=(7.6, 7.6))
    ax.barh(range(len(top)), top.stability * 100, color=cols,
            edgecolor="white", linewidth=.6)
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(top.symbol, fontsize=9)
    ax.set_xlabel(f"% of {len(rdf)} threshold settings in which the gene "
                  f"reaches the top {TOP_N}")
    ax.set_title("Which candidates survive our own threshold choices\n"
                 "highlighted = currently in the top 10", loc="left")
    ax.set_xlim(0, 108)
    for i, v in enumerate(top.stability * 100):
        ax.text(v + 1.4, i, f"{v:.0f}%", va="center", fontsize=8)
    fig.tight_layout()
    p = save_figure(fig, "step9_stability")
    print(f"\n  wrote {Path(p).relative_to(ROOT) if p else 'figure'}")

    ok = reproduces and rerun_ok and collapse_ok
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
