"""Step 7: apply the pre-registered selection rule (D52).

The rule was written and dated in docs/decisions.md BEFORE this script was first
run, which is the point of the runbook's rule 3 and is verified by check 7.4.
Nothing here is tuned to the output.

  Eligible      surface tier S1/S2/S3 (not Flagged)
                passes the RNA defaults against GTEx healthy in >=1 cancer type
                tumour median >= 5 TPM
  Rank by       the specificity score S, highest first
  Per cancer    top N within each cancer type
  Global top 10 each gene's best cancer type, ranked by S, family limit applied
  Family limit  at most 2 per HGNC gene group; a gene in several groups is
                capped by its MOST SPECIFIC group (fewest members); a gene with
                no group is never capped
  Reported but never ranked on: protein concordance (D48), comparator grade,
                low_power, trajectory, CSPA/SURFY, the S variants

Writes results/step7_ranking_full.csv, step7_per_cancer_type.csv,
step7_top10.csv, step7_reserve.csv and step7_counts.txt.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
RESULTS = ROOT / "results"

TOP_N_PER_CANCER = 10
GLOBAL_TOP = 10
RESERVE = 5
MAX_PER_FAMILY = 2
ELIGIBLE_TIERS = {"S1", "S2", "S3"}


def gene_groups() -> tuple[dict[str, str], dict[str, int]]:
    """symbol -> its most specific HGNC gene group, plus group sizes."""
    h = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", low_memory=False)
    h = h[["symbol", "gene_group"]].dropna()
    pairs = []
    for sym, grp in zip(h.symbol, h.gene_group):
        for g in str(grp).split("|"):
            g = g.strip()
            if g:
                pairs.append((sym, g))
    df = pd.DataFrame(pairs, columns=["symbol", "group"])
    sizes = df.group.value_counts().to_dict()
    # A gene in several groups is capped by its most specific one, i.e. the
    # group with fewest members. CEACAM5 sits in four groups; without this the
    # cap would be driven by a vast bucket like "CD molecules".
    df["size"] = df.group.map(sizes)
    best = df.sort_values("size").drop_duplicates("symbol").set_index("symbol")
    return best.group.to_dict(), sizes


def apply_family_limit(df: pd.DataFrame, n: int) -> pd.DataFrame:
    """Take the top n rows, allowing at most MAX_PER_FAMILY per gene group."""
    taken, counts, dropped = [], {}, []
    for _, r in df.iterrows():
        g = r.gene_group
        if isinstance(g, str) and g:
            if counts.get(g, 0) >= MAX_PER_FAMILY:
                dropped.append((r.symbol, g))
                continue
            counts[g] = counts.get(g, 0) + 1
        taken.append(r)
        if len(taken) >= n:
            break
    out = pd.DataFrame(taken)
    out.attrs["dropped_by_family"] = dropped
    return out


def add_homogeneity(elig: pd.DataFrame) -> pd.DataFrame:
    """Fraction of tumours expressing the gene above 1 and 10 TPM, and above
    the gene's own worst healthy tissue.

    The 10 TPM line is a convention. The worst healthy tissue is the bar a drug
    actually has to clear, and it differs per gene, so the last column answers
    "in what share of patients would this target have any therapeutic window at
    all". Reported only — it never enters the ranking (D52 is pre-registered).
    """
    mat = pd.read_parquet(ROOT / "data" / "processed" / "step4_candidate_tpm.parquet")
    samp = pd.read_parquet(ROOT / "data" / "processed" / "step4_samples.parquet")
    tum = samp[(samp.study == "TCGA") & (samp.sample_type == "Primary Tumor")]
    cols = {c: i for i, c in enumerate(mat.columns)}
    V = mat.to_numpy(dtype=np.float32)
    gidx = {g: i for i, g in enumerate(mat.index)}

    f1, f10, q25, fh = [], [], [], []
    ids_by_cancer = {c: [cols[s] for s in g["sample"] if s in cols]
                     for c, g in tum.groupby("disease")}
    for gene, cancer, h in zip(elig.ensembl_gene, elig.cancer_type,
                               elig.healthy_max_log2):
        ids = ids_by_cancer.get(cancer, [])
        gi = gidx.get(gene)
        if gi is None or not ids:
            f1.append(np.nan); f10.append(np.nan); q25.append(np.nan)
            fh.append(np.nan); continue
        vals = V[gi, ids]                      # log2(TPM+1), per tumour
        tpm = np.power(2.0, vals) - 1.0
        f1.append(float((tpm > 1).mean()))
        f10.append(float((tpm > 10).mean()))
        q25.append(float(np.percentile(tpm, 25)))
        # compared on the same log2(TPM+1) scale the healthy maximum is stored on
        fh.append(float((vals > h).mean()) if np.isfinite(h) else np.nan)
    elig["frac_tumours_above_1tpm"] = f1
    elig["frac_tumours_above_10tpm"] = f10
    elig["tumour_tpm_q25"] = q25
    elig["frac_tumours_above_worst_healthy"] = fh
    return elig


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    s5 = pd.read_csv(RESULTS / "step5_specificity.csv")

    # ---- eligibility (D52) ----
    elig = s5[s5.tier.isin(ELIGIBLE_TIERS)
              & s5.passes_default
              & s5.expressed_in_tumour].copy()
    n_rows_in = len(s5)

    grp, sizes = gene_groups()
    elig["gene_group"] = elig.symbol.map(grp)
    elig["gene_group_size"] = elig.gene_group.map(sizes)

    # ---- context columns: reported, never ranked on ----
    try:
        traj = pd.read_csv(RESULTS / "step4_trajectory.csv")[
            ["ensembl_gene", "cancer_type", "trajectory", "field_step"]]
        elig = elig.merge(traj, on=["ensembl_gene", "cancer_type"], how="left")
    except FileNotFoundError:
        elig["trajectory"] = np.nan
    try:
        con = pd.read_csv(RESULTS / "step6_rna_protein_concordance.csv")
        con = con[["ensembl_gene", "tcga_cancer_type", "agreement",
                   "protein_log2fc", "protein_fdr"]].rename(
            columns={"tcga_cancer_type": "cancer_type"})
        elig = elig.merge(con, on=["ensembl_gene", "cancer_type"], how="left")
        elig["agreement"] = elig.agreement.fillna("protein_not_measured")
    except FileNotFoundError:
        elig["agreement"] = "protein_not_measured"

    # ---- expression homogeneity, reported not ranked on ----
    # A median hides whether a target is on most tumour cells or only a subset.
    # DCSTAMP in thyroid has a median of 94 TPM but a 25th percentile of 0.9:
    # expressed strongly in some tumours, absent in others. A drug reaches only
    # the first group. This is the "high AND homogeneous" criterion the
    # surfaceome literature asks for, and a median alone cannot show it.
    elig = add_homogeneity(elig)

    elig = elig.sort_values("S", ascending=False, kind="mergesort")
    elig.to_csv(RESULTS / "step7_ranking_full.csv", index=False)

    # ---- per cancer type ----
    per = (elig.sort_values("S", ascending=False, kind="mergesort")
           .groupby("cancer_type", sort=True)
           .head(TOP_N_PER_CANCER))
    per.to_csv(RESULTS / "step7_per_cancer_type.csv", index=False)

    # ---- global: each gene's best cancer type ----
    best = elig.loc[elig.groupby("ensembl_gene").S.idxmax()].copy()
    best["n_cancer_types_passing"] = best.ensembl_gene.map(
        elig.groupby("ensembl_gene").cancer_type.nunique())
    best = best.sort_values("S", ascending=False, kind="mergesort")

    picked = apply_family_limit(best, GLOBAL_TOP + RESERVE)
    top10 = picked.head(GLOBAL_TOP)
    reserve = picked.iloc[GLOBAL_TOP:GLOBAL_TOP + RESERVE]
    top10.to_csv(RESULTS / "step7_top10.csv", index=False)
    reserve.to_csv(RESULTS / "step7_reserve.csv", index=False)

    # ---- report ----
    lines = [
        "STEP 7 — ranking under the pre-registered rule (D52)",
        "Generated by src/step7_rank.py — do not edit by hand.",
        "The rule was recorded in docs/decisions.md before this script first ran.",
        "",
        "funnel:",
        f"  gene x cancer-type rows in             : {n_rows_in:,}",
        f"  eligible (tier S1/S2/S3, passes RNA,",
        f"            >= 5 TPM in tumour)          : {len(elig):,}",
        f"  distinct genes eligible                : {elig.ensembl_gene.nunique():,}",
        f"  cancer types represented               : {elig.cancer_type.nunique()}",
        "",
        f"per cancer type: top {TOP_N_PER_CANCER} each            : {len(per):,} rows",
        f"global top {GLOBAL_TOP} + reserve {RESERVE}               : "
        f"{len(top10)} + {len(reserve)}",
        "",
        "dropped by the 2-per-family limit:",
    ]
    dropped = picked.attrs.get("dropped_by_family", [])
    if dropped:
        for sym, g in dropped[:12]:
            lines.append(f"  {sym:<12} third or later member of '{g}'")
    else:
        lines.append("  none")

    lines += ["", "GLOBAL TOP 10", "",
              f"  {'#':<3} {'gene':<10} {'tier':<5} {'cancer type':<34} "
              f"{'S':>6} {'log2FC':>7} {'tumTPM':>8} {'grade':<13} "
              f"{'n_types':>7}  worst healthy tissue"]
    for i, (_, r) in enumerate(top10.iterrows(), 1):
        lines.append(
            f"  {i:<3} {str(r.symbol):<10} {r.tier:<5} {r.cancer_type[:34]:<34} "
            f"{r.S:>6.2f} {r.log2fc:>7.2f} {r.tumour_tpm:>8.1f} "
            f"{str(r.comparator_grade):<13} {int(r.n_cancer_types_passing):>7}  "
            f"{str(r.healthy_max_tissue)[:24]}")

    lines += ["", "RESERVE (ranks 11-15)", ""]
    for i, (_, r) in enumerate(reserve.iterrows(), GLOBAL_TOP + 1):
        lines.append(f"  {i:<3} {str(r.symbol):<10} {r.tier:<5} "
                     f"{r.cancer_type[:34]:<34} {r.S:>6.2f}")

    lines += ["", "context carried but NOT ranked on (D48, D53):", ""]
    lines.append(f"  protein agreement among the top 10: "
                 f"{top10.agreement.value_counts().to_dict()}")
    lines.append(f"  comparator grades among the top 10: "
                 f"{top10.comparator_grade.value_counts().to_dict()}")
    if "trajectory" in top10.columns:
        lines.append(f"  trajectory among the top 10      : "
                     f"{top10.trajectory.dropna().value_counts().to_dict()}")
    lines.append(f"  low_power flagged                : "
                 f"{int(top10.low_power.sum())} of {len(top10)}")
    lines += ["", "expression homogeneity of the top 10 "
              "(what fraction of tumours actually carry it):", ""]
    for i, (_, r) in enumerate(top10.iterrows(), 1):
        lines.append(f"  {i:<3} {str(r.symbol):<10} "
                     f">1 TPM in {r.frac_tumours_above_1tpm:>5.0%} of tumours, "
                     f">10 TPM in {r.frac_tumours_above_10tpm:>5.0%}, "
                     f"lower quartile {r.tumour_tpm_q25:>8.1f} TPM")

    report = "\n".join(lines)
    print(report)
    (RESULTS / "step7_counts.txt").write_text(report + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
