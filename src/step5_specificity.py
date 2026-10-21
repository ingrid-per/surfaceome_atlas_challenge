"""Step 5: how absent is each candidate from healthy tissue?

This is the safety half of the question. Step 4 asked "is it up in tumours";
this asks "is it anywhere it must not be".

Decisions implemented:

  D42  No vital-organ list. Following ImmunoTar: score on the MAXIMUM across
       tissues, count the BREADTH of tissues above 10 TPM, and compute a
       brain-excluded variant because antibodies and CAR-T cross the
       blood-brain barrier poorly. The tissue giving the maximum is reported.
  D44  Maximum, not a percentile. The 90th percentile is reported alongside
       but never scores: with ~50 tissues it discards the top five, and for
       MSLN those five include lung, the organ where mesothelin drugs are
       actually toxic.
  D43  The tumour must carry at least 5 TPM, or there is nothing to bind.
  D2/D9 Two independent healthy sources — Toil-GTEx (same pipeline as the
       tumours) and HPA nTPM (different pipeline entirely). The more cautious,
       i.e. HIGHER, of the two is used to score.

  S = log2(tumour median TPM + 1) - log2(max healthy TPM + 1)

Writes results/step5_specificity.csv and the per-gene healthy profile.
"""

from __future__ import annotations

import io
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

HIGH_TPM = 10.0     # D42, ImmunoTar's 90th-percentile cut-off for "high"
MIN_TUMOUR_TPM = 5.0  # D43

# Tissue sets that get their own S variant (D46). None of these is excluded
# from the headline S; each is reported in parallel so Step 7 can choose.
#   brain      - antibodies and CAR-T cross the blood-brain barrier poorly (D42)
#   testis     - immune-privileged; the cancer-testis antigen class exists
#                precisely because of this
#   placenta   - simply not present in the patient; oncofetal antigens such as
#                GPC3 score badly against it for no clinical reason
GTEX_BRAIN_PREFIX = "Brain -"
BRAIN = {"cerebral cortex", "choroid plexus", "Pituitary"}
REPRO_FETAL = {"placenta", "Testis", "testis"}


def tissue_sets(columns) -> dict:
    brain = {c for c in columns
             if c.startswith(GTEX_BRAIN_PREFIX) or c in BRAIN}
    rf = {c for c in columns if c in REPRO_FETAL}
    return {"nobrain": brain,
            "no_reproductive_fetal": rf,
            "relaxed": brain | rf}


def log2p1(tpm: np.ndarray) -> np.ndarray:
    return np.log2(tpm + 1.0)


def tau(values: np.ndarray) -> float:
    """Tissue-specificity index: 0 = everywhere, 1 = one tissue only."""
    v = values[~np.isnan(values)]
    if v.size < 2 or v.max() <= 0:
        return np.nan
    return float(np.sum(1.0 - v / v.max()) / (v.size - 1))


def summarise(profile: pd.DataFrame, tag: str) -> pd.DataFrame:
    """profile: genes x tissues, in log2(TPM+1). One row per gene."""
    vals = profile.to_numpy()
    tissues = np.array(profile.columns)
    tpm = np.power(2.0, vals) - 1.0
    imax = np.nanargmax(vals, axis=1)

    data = {
        "gene": profile.index,
        f"{tag}_max_log2": vals[np.arange(len(vals)), imax],
        f"{tag}_max_tissue": tissues[imax],
        f"{tag}_p90_log2": np.nanpercentile(vals, 90, axis=1),
        f"{tag}_median_log2": np.nanmedian(vals, axis=1),
        f"{tag}_n_tissues_above_10tpm": (tpm >= HIGH_TPM).sum(axis=1),
        f"{tag}_n_tissues": vals.shape[1],
        f"{tag}_tau": [tau(r) for r in tpm],
    }
    # D46: one maximum per exclusion variant, each reported separately
    for name, drop in tissue_sets(profile.columns).items():
        keep = [i for i, c in enumerate(tissues) if c not in drop]
        sub = vals[:, keep]
        j = np.nanargmax(sub, axis=1)
        data[f"{tag}_max_{name}_log2"] = sub[np.arange(len(sub)), j]
        data[f"{tag}_max_{name}_tissue"] = np.array(tissues[keep])[j]
    return pd.DataFrame(data).set_index("gene")


def gtex_profile() -> pd.DataFrame:
    mat = pd.read_parquet(PROC / "step4_candidate_tpm.parquet")
    samp = pd.read_parquet(PROC / "step4_samples.parquet")
    g = samp[(samp.study == "GTEX") & (samp.sample_type == "Normal Tissue")]
    cols = {}
    for tissue, grp in g.groupby("disease"):
        ids = [s for s in grp["sample"] if s in mat.columns]
        if ids:
            cols[tissue] = mat[ids].median(axis=1)
    return pd.DataFrame(cols)


def hpa_profile(genes: set[str]) -> pd.DataFrame:
    with zipfile.ZipFile(RAW / "rna_tissue_hpa.tsv.zip") as z, \
            z.open(z.namelist()[0]) as fh:
        h = pd.read_csv(io.TextIOWrapper(fh, encoding="utf8"), sep="\t")
    h = h[h.Gene.isin(genes)]
    wide = h.pivot_table(index="Gene", columns="Tissue", values="nTPM",
                         aggfunc="median")
    return np.log2(wide + 1.0)


def main() -> int:
    RESULTS.mkdir(exist_ok=True)

    gt = gtex_profile()
    print(f"GTEx healthy profile : {gt.shape[0]:,} genes x {gt.shape[1]} tissues")
    hp = hpa_profile(set(gt.index))
    print(f"HPA healthy profile  : {hp.shape[0]:,} genes x {hp.shape[1]} tissues")

    g_sum = summarise(gt, "gtex")
    h_sum = summarise(hp, "hpa")
    healthy = g_sum.join(h_sum, how="left")

    # D2/D9: two independent pipelines; take the more cautious (higher) value.
    pair = healthy[["gtex_max_log2", "hpa_max_log2"]]
    healthy["healthy_max_log2"] = pair.max(axis=1)
    healthy["healthy_max_source"] = np.where(
        healthy.hpa_max_log2.isna(), "gtex_only",
        np.where(healthy.hpa_max_log2 > healthy.gtex_max_log2, "hpa", "gtex"))
    healthy["healthy_max_tissue"] = np.where(
        healthy.healthy_max_source == "hpa",
        healthy.hpa_max_tissue, healthy.gtex_max_tissue)
    healthy["healthy_max_relaxed_tissue"] = np.where(
        healthy.hpa_max_relaxed_log2 > healthy.gtex_max_relaxed_log2,
        healthy.hpa_max_relaxed_tissue, healthy.gtex_max_relaxed_tissue)
    for name in ("nobrain", "no_reproductive_fetal", "relaxed"):
        healthy[f"healthy_max_{name}_log2"] = healthy[
            [f"gtex_max_{name}_log2", f"hpa_max_{name}_log2"]].max(axis=1)
    healthy["healthy_breadth_above_10tpm"] = healthy[
        ["gtex_n_tissues_above_10tpm", "hpa_n_tissues_above_10tpm"]].max(axis=1)
    healthy.to_csv(RESULTS / "step5_healthy_profile.csv")

    # ---- join to the tumour side (GTEx arm of Step 4) ----
    d = pd.read_csv(RESULTS / "step4_tumor_vs_normal.csv")
    d = d[d.arm == "gtex"].copy()
    d = d.merge(healthy, left_on="ensembl_gene", right_index=True, how="left")

    d["tumour_tpm"] = np.power(2.0, d.median_tumour_log2tpm1) - 1.0
    d["expressed_in_tumour"] = d.tumour_tpm >= MIN_TUMOUR_TPM          # D43
    # D44/D46: the headline S plus one variant per tissue exclusion. No
    # threshold is applied to any of them (D45) - Step 7 decides how to use them.
    d["S"] = d.median_tumour_log2tpm1 - d.healthy_max_log2
    for name in ("nobrain", "no_reproductive_fetal", "relaxed"):
        d[f"S_{name}"] = (d.median_tumour_log2tpm1
                          - d[f"healthy_max_{name}_log2"])
    d["S_p90_crosscheck"] = d.median_tumour_log2tpm1 - d[
        ["gtex_p90_log2", "hpa_p90_log2"]].max(axis=1)

    # D45: Step 5's only hard filter is the 5 TPM floor. It is not a judgment -
    # below it there is nothing for a drug to bind.
    d["carried_forward"] = d.passes_default & d.expressed_in_tumour

    cols = ["ensembl_gene", "symbol", "tier", "cancer_type", "comparator_grade",
            "median_tumour_log2tpm1", "tumour_tpm", "log2fc", "fdr",
            "passes_default", "expressed_in_tumour",
            "healthy_max_log2", "healthy_max_tissue", "healthy_max_source",
            "healthy_max_nobrain_log2", "healthy_max_no_reproductive_fetal_log2",
            "healthy_max_relaxed_log2", "healthy_max_relaxed_tissue",
            "healthy_breadth_above_10tpm",
            "gtex_max_log2", "gtex_max_tissue", "gtex_tau",
            "hpa_max_log2", "hpa_max_tissue",
            "S", "S_nobrain", "S_no_reproductive_fetal", "S_relaxed",
            "S_p90_crosscheck", "carried_forward",
            "low_power", "gpi", "cspa", "surfy"]
    out = d[cols].sort_values("S", ascending=False)
    out.to_csv(RESULTS / "step5_specificity.csv", index=False)

    # ---- P6b: pick the ubiquitous-surface control BY RULE ----
    cand = pd.read_csv(PROC / "surface_candidates.csv")
    pm = (cand[cand.is_candidate]
          .groupby("ensembl_gene")
          .agg(pm=("plasma_membrane_evidence", "any")))
    ctrl = healthy.join(
        d.drop_duplicates("ensembl_gene").set_index("ensembl_gene")[["symbol", "tier"]],
        how="inner").join(pm, how="left")
    # Rule: surface tier S1 AND expressed above 10 TPM in EVERY healthy tissue.
    # (tau was tried first but never falls below 0.36 across 51 tissues, so it
    # does not separate ubiquitous genes here; breadth does, cleanly.)
    n_tis = int(ctrl.gtex_n_tissues.max())
    rule = ctrl[(ctrl.tier == "S1")
                & ctrl.pm.fillna(False)          # positive plasma-membrane evidence
                & (ctrl.gtex_n_tissues_above_10tpm >= n_tis)]
    rule = rule.sort_values("gtex_tau")

    # ---- report ----
    passing = out[out.carried_forward]
    lines = [
        "STEP 5 — specificity against healthy tissue",
        "Generated by src/step5_specificity.py — do not edit by hand.",
        f"S = log2(tumour median TPM+1) - log2(max healthy TPM+1).  "
        f"Max, not percentile (D44). No vital-organ list (D42).",
        f"Tumour must reach {MIN_TUMOUR_TPM:g} TPM (D43). "
        f"Breadth counts tissues above {HIGH_TPM:g} TPM.",
        "",
        f"healthy tissues, GTEx                    : {gt.shape[1]}",
        f"healthy tissues, HPA                     : {hp.shape[1]}",
        f"genes with both sources                  : "
        f"{int(healthy.hpa_max_log2.notna().sum()):,} of {len(healthy):,}",
        f"  healthy maximum came from GTEx         : "
        f"{int((healthy.healthy_max_source == 'gtex').sum()):,}",
        f"  ...from HPA                            : "
        f"{int((healthy.healthy_max_source == 'hpa').sum()):,}",
        f"  ...GTEx only (no HPA row)              : "
        f"{int((healthy.healthy_max_source == 'gtex_only').sum()):,}",
        "",
        "funnel, gene x cancer-type rows (GTEx arm):",
        f"  tested                                 : {len(out):,}",
        f"  pass Step 4 thresholds                 : {int(out.passes_default.sum()):,}",
        f"  ...and reach {MIN_TUMOUR_TPM:g} TPM in tumour (D43)    : "
        f"{int(out.carried_forward.sum()):,}   <- carried to Step 6",
        "",
        f"distinct genes carried forward           : "
        f"{out.loc[out.carried_forward, 'ensembl_gene'].nunique():,}",
        "",
        "No threshold is applied to S (D45). For reference only, how many of the",
        "carried-forward rows would clear S > 0 under each tissue exclusion:",
        f"  S           all tissues                : "
        f"{int((out.carried_forward & (out.S > 0)).sum()):,}",
        f"  S_nobrain   brain excluded             : "
        f"{int((out.carried_forward & (out.S_nobrain > 0)).sum()):,}",
        f"  S_no_reproductive_fetal                : "
        f"{int((out.carried_forward & (out.S_no_reproductive_fetal > 0)).sum()):,}",
        f"  S_relaxed   brain+testis+placenta out  : "
        f"{int((out.carried_forward & (out.S_relaxed > 0)).sum()):,}",
        "",
        "effect of using the maximum rather than the 90th percentile:",
        f"  carried-forward rows with S <= 0 but S_p90 > 0 : "
        f"{int((out.carried_forward & (out.S <= 0) & (out.S_p90_crosscheck > 0)).sum()):,}",
        "    (what a 90th-percentile summary would have called clean)",
        "",
        "P6b — ubiquitous-surface control, selected by rule:",
        "  rule: surface tier S1, positive plasma-membrane evidence (D31), and",
        "        above 10 TPM in ALL 51 GTEx tissues; ranked by tau",
        "        (lowest tau = most uniformly expressed)",
    ]
    for _, r in rule.head(8).iterrows():
        lines.append(f"    {str(r.symbol):<10} tau={r.gtex_tau:.3f}  "
                     f"tissues>10TPM={int(r.healthy_breadth_above_10tpm):>2}  "
                     f"max={r.healthy_max_log2:.2f} ({r.healthy_max_tissue})")
    lines += [
        "",
        "top 15 by specificity score S:",
        f"  {'gene':<10} {'tier':<8} {'cancer type':<32} {'S':>6} {'tumTPM':>8} "
        f"{'worst healthy tissue':<26}",
    ]
    for _, r in passing.head(15).iterrows():
        lines.append(f"  {str(r.symbol):<10} {r.tier:<8} {r.cancer_type[:32]:<32} "
                     f"{r.S:>6.2f} {r.tumour_tpm:>8.1f} "
                     f"{str(r.healthy_max_tissue)[:26]:<26}")

    report = "\n".join(lines)
    print("\n" + report)
    (RESULTS / "step5_counts.txt").write_text(report + "\n")
    print(f"\nWrote results/step5_specificity.csv ({len(out):,} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
