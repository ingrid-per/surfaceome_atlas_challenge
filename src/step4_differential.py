"""Step 4b: Wilcoxon tumour vs normal, per cancer type and per arm.

Decisions implemented here:

  D33  Comparator grade travels with every result (direct / imperfect /
       lineage_only). It is data, not a filter.
  D34  GTEx is the primary arm; TCGA adjacent-normal is a SEPARATE arm. Never
       pooled — pooling would confound a batch effect with the tumour signal,
       and Aran et al. show the two normals are biologically different states.
  D35  Ovarian runs a second GTEx arm (Fallopian Tube).
  D37  Minimum 10 samples per group; `low_power` under 30.
  D38  Wilcoxon rank-sum (Mann-Whitney U).
  D39  BH-FDR within each (cancer type, arm) group, never pooled. Defaults
       FDR < 0.01 and log2FC >= 2.

Every comparison is one vectorised scipy call over all genes at once, so this
runs in seconds rather than looping 250,000 times.

Writes results/step4_tumor_vs_normal.csv (long) and results/step4_counts.txt.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
DOCS = ROOT / "docs"
RESULTS = ROOT / "results"

MIN_N = 10        # D37
LOW_POWER_N = 30  # D37
FDR_MAX = 0.01    # D39
LOG2FC_MIN = 2.0  # D39

TIER_RANK = {"S1": 0, "S2": 1, "S3": 2, "Flagged": 3}


def gene_annotation() -> pd.DataFrame:
    """One row per Ensembl gene: its accessions, symbol and best surface tier.

    A gene can carry several accessions (D15 keeps one-to-many mappings), so the
    tier reported is the strongest of them and every accession is listed.
    """
    cand = pd.read_csv(PROC / "surface_candidates.csv")
    cand = cand[cand.is_candidate]
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)
    link = gmap[gmap.accession.isin(set(cand.accession))][
        ["accession", "ensembl_gene"]].dropna().drop_duplicates()
    j = link.merge(
        cand[["accession", "symbol", "tier", "comparator_unused"]
             if "comparator_unused" in cand.columns
             else ["accession", "symbol", "tier", "gpi_any", "cspa", "surfy"]],
        on="accession")
    j["rank"] = j.tier.map(TIER_RANK)
    j = j.sort_values(["ensembl_gene", "rank"], kind="mergesort")
    out = (j.groupby("ensembl_gene")
             .agg(symbol=("symbol", "first"),
                  tier=("tier", "first"),
                  accessions=("accession", lambda s: ";".join(sorted(set(s)))),
                  gpi=("gpi_any", "any"), cspa=("cspa", "any"),
                  surfy=("surfy", "any"))
             .reset_index())
    return out


def comparisons(samples: pd.DataFrame) -> list[dict]:
    """Build the comparison list from the reviewed tissue map (D33-D35)."""
    tmap = pd.read_csv(DOCS / "step4_tissue_map_DRAFT.tsv", sep="\t")
    tumour = samples[(samples.study == "TCGA")
                     & (samples.sample_type == "Primary Tumor")]
    gtex = samples[(samples.study == "GTEX")
                   & (samples.sample_type == "Normal Tissue")]
    nat = samples[(samples.study == "TCGA")
                  & (samples.sample_type == "Solid Tissue Normal")]

    out = []
    for _, r in tmap.iterrows():
        cancer = r.tcga_cancer_type
        t_ids = tumour.loc[tumour.disease == cancer, "sample"].tolist()

        def add(arm: str, ids: list[str], label: str) -> None:
            out.append({"cancer_type": cancer, "arm": arm,
                        "normal_label": label,
                        "comparator_grade": r.comparator_grade,
                        "tumour_ids": t_ids, "normal_ids": ids})

        for col, arm in (("gtex_tissues", "gtex"), ("gtex_alt_tissues", "gtex_alt")):
            spec = r[col]
            if isinstance(spec, str) and spec.strip():
                tissues = [x.strip() for x in spec.split("|") if x.strip()]
                ids = gtex.loc[gtex.disease.isin(tissues), "sample"].tolist()
                add(arm, ids, " | ".join(tissues))

        ids = nat.loc[nat.disease == cancer, "sample"].tolist()
        if ids:
            add("tcga_nat", ids, f"{cancer} adjacent normal")
    return out


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    mat = pd.read_parquet(PROC / "step4_candidate_tpm.parquet")
    samples = pd.read_parquet(PROC / "step4_samples.parquet")
    ann = gene_annotation()
    print(f"genes {len(mat):,}   samples {mat.shape[1]:,}   "
          f"annotated {len(ann):,}")

    col_pos = {c: i for i, c in enumerate(mat.columns)}
    values = mat.to_numpy(dtype=np.float32)
    genes = mat.index.to_numpy()

    results, log = [], []
    for c in comparisons(samples):
        ti = [col_pos[s] for s in c["tumour_ids"] if s in col_pos]
        ni = [col_pos[s] for s in c["normal_ids"] if s in col_pos]
        tag = f"{c['cancer_type']} [{c['arm']}]"

        if len(ti) < MIN_N or len(ni) < MIN_N:
            why = (f"tumour n={len(ti)}" if len(ti) < MIN_N
                   else f"normal n={len(ni)}") + f" < {MIN_N} (D37)"
            log.append({"cancer_type": c["cancer_type"], "arm": c["arm"],
                        "normal_label": c["normal_label"],
                        "grade": c["comparator_grade"],
                        "n_tumour": len(ti), "n_normal": len(ni),
                        "genes_tested": 0, "genes_passing": 0,
                        "status": "SKIPPED", "reason": why})
            print(f"  skip  {tag:<58} {why}")
            continue

        a = values[:, ti]
        b = values[:, ni]
        # asymptotic applies the tie correction that 20.5% exact zeros demands
        res = mannwhitneyu(a, b, axis=1, alternative="two-sided",
                           method="asymptotic")
        p = np.asarray(res.pvalue, dtype=np.float64)
        p = np.where(np.isfinite(p), p, 1.0)
        fdr = multipletests(p, method="fdr_bh")[1]   # D39: within this group

        med_t = np.median(a, axis=1)
        med_n = np.median(b, axis=1)
        lfc = med_t - med_n     # already log2(TPM+1), so a difference is a log2FC

        low_power = len(ti) < LOW_POWER_N or len(ni) < LOW_POWER_N
        passes = (fdr < FDR_MAX) & (lfc >= LOG2FC_MIN)

        results.append(pd.DataFrame({
            "ensembl_gene": genes,
            "cancer_type": c["cancer_type"],
            "arm": c["arm"],
            "normal_label": c["normal_label"],
            "comparator_grade": c["comparator_grade"],
            "n_tumour": len(ti),
            "n_normal": len(ni),
            "median_tumour_log2tpm1": med_t,
            "median_normal_log2tpm1": med_n,
            "log2fc": lfc,
            "p_value": p,
            "fdr": fdr,
            "low_power": low_power,
            "passes_default": passes,
        }))
        log.append({"cancer_type": c["cancer_type"], "arm": c["arm"],
                    "normal_label": c["normal_label"],
                    "grade": c["comparator_grade"],
                    "n_tumour": len(ti), "n_normal": len(ni),
                    "genes_tested": len(genes), "genes_passing": int(passes.sum()),
                    "status": "low_power" if low_power else "ok", "reason": ""})
        print(f"  ran   {tag:<58} T={len(ti):>4} N={len(ni):>4}  "
              f"pass={int(passes.sum()):>4}{'  (low power)' if low_power else ''}")

    df = pd.concat(results, ignore_index=True)
    df = df.merge(ann, on="ensembl_gene", how="left")
    cols = ["ensembl_gene", "symbol", "tier", "accessions", "cancer_type", "arm",
            "normal_label", "comparator_grade", "n_tumour", "n_normal",
            "median_tumour_log2tpm1", "median_normal_log2tpm1", "log2fc",
            "p_value", "fdr", "low_power", "passes_default",
            "gpi", "cspa", "surfy"]
    df = df[cols].sort_values(["cancer_type", "arm", "fdr"], kind="mergesort")
    df.to_csv(RESULTS / "step4_tumor_vs_normal.csv", index=False)

    logdf = pd.DataFrame(log)
    logdf.to_csv(RESULTS / "step4_comparisons.csv", index=False)

    ran = logdf[logdf.status != "SKIPPED"]
    lines = [
        "STEP 4b — tumour vs normal, RNA",
        "Generated by src/step4_differential.py — do not edit by hand.",
        f"Wilcoxon rank-sum (D38); BH-FDR within each cancer type and arm (D39);",
        f"defaults FDR < {FDR_MAX} and log2FC >= {LOG2FC_MIN} (D39); "
        f"minimum n = {MIN_N} per group (D37).",
        "",
        f"genes tested                             : {len(mat):,}",
        f"comparisons attempted                    : {len(logdf)}",
        f"  ran                                    : {len(ran)}",
        f"  skipped for too few samples            : "
        f"{int((logdf.status == 'SKIPPED').sum())}",
        f"  flagged low_power (a group under {LOW_POWER_N})    : "
        f"{int((logdf.status == 'low_power').sum())}",
        "",
        f"result rows (gene x cancer type x arm)   : {len(df):,}",
        f"rows passing the defaults                : {int(df.passes_default.sum()):,}",
        f"distinct genes passing anywhere          : "
        f"{df.loc[df.passes_default, 'ensembl_gene'].nunique():,}",
        "",
        "by arm:",
    ]
    for arm, grp in df.groupby("arm"):
        lines.append(f"  {arm:<10} rows={len(grp):>7,}  passing={int(grp.passes_default.sum()):>5,}"
                     f"  genes={grp.loc[grp.passes_default,'ensembl_gene'].nunique():>5,}")
    lines += ["", "by comparator grade:"]
    for g, grp in df.groupby("comparator_grade"):
        lines.append(f"  {g:<14} rows={len(grp):>7,}  "
                     f"passing={int(grp.passes_default.sum()):>5,}")
    lines += ["", "per comparison:", ""]
    lines.append(f"  {'cancer type':<40} {'arm':<9} {'T':>5} {'N':>5} "
                 f"{'tested':>7} {'pass':>5}  grade / status")
    for _, r in logdf.iterrows():
        lines.append(f"  {r.cancer_type[:40]:<40} {r.arm:<9} {r.n_tumour:>5} "
                     f"{r.n_normal:>5} {r.genes_tested:>7} {r.genes_passing:>5}  "
                     f"{r.grade} / {r.status} {r.reason}")

    report = "\n".join(lines)
    print("\n" + "\n".join(lines[:28]))
    (RESULTS / "step4_counts.txt").write_text(report + "\n")
    print(f"\nWrote results/step4_tumor_vs_normal.csv ({len(df):,} rows)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
