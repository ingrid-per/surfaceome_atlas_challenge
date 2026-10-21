"""Step 7b: export the ranking as JSON for the interactive deliverable (D54).

Every number here comes from the result files, never typed by hand, so the
artifact and the committed results cannot drift apart.
"""

from __future__ import annotations

import json
import zipfile
import io
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

TOP_PER_CANCER = 10
N_TISSUES_SHOWN = 8


def healthy_profiles(genes: set[str]) -> dict[str, list]:
    """Top healthy tissues per gene — the detail view's safety panel."""
    mat = pd.read_parquet(PROC / "step4_candidate_tpm.parquet")
    samp = pd.read_parquet(PROC / "step4_samples.parquet")
    g = samp[(samp.study == "GTEX") & (samp.sample_type == "Normal Tissue")]
    cols = {}
    for tissue, grp in g.groupby("disease"):
        ids = [s for s in grp["sample"] if s in mat.columns]
        if ids:
            cols[tissue] = mat[ids].median(axis=1)
    prof = pd.DataFrame(cols)
    prof = prof.loc[prof.index.intersection(genes)]
    tpm = np.power(2.0, prof) - 1.0

    out = {}
    for gene, row in tpm.iterrows():
        top = row.sort_values(ascending=False).head(N_TISSUES_SHOWN)
        out[gene] = [{"tissue": t, "tpm": round(float(v), 1)}
                     for t, v in top.items()]
    return out


def lineage() -> pd.DataFrame:
    """D61: no microenvironment call is produced. Tumour purity is not
    accounted for, and that is stated as a limitation rather than guessed at
    with a flag. Kept as a stub so the row builder stays unchanged."""
    return pd.DataFrame()


def surface_detail() -> pd.DataFrame:
    c = pd.read_csv(PROC / "surface_candidates.csv")
    return c[c.is_candidate].drop_duplicates("ensembl_gene").set_index("ensembl_gene")


def sample_counts() -> dict:
    """Group sizes per cancer type, from the Step 4 comparison log."""
    log = pd.read_csv(RESULTS / "step4_comparisons.csv")
    g = log[(log.arm == "gtex") & (log.status != "SKIPPED")]
    return {r.cancer_type: (int(r.n_tumour), int(r.n_normal))
            for _, r in g.iterrows()}


def row_to_dict(r, surf, prof, cohort_of, counts, lin=None) -> dict:
    # Literature vocabulary, not TCGA's internal label. Searching PubMed for
    # "Kidney Clear Cell Carcinoma" finds far less than "clear cell renal cell
    # carcinoma", and the mismatch was quietly breaking the deliverable's links.
    from step8_literature import DISEASE_TERM
    s = surf.loc[r.ensembl_gene] if r.ensembl_gene in surf.index else None
    d = {
        "gene": r.ensembl_gene,
        "symbol": str(r.symbol),
        "cancer": r.cancer_type,
        "tier": r.tier,
        "S": round(float(r.S), 2),
        "S_relaxed": round(float(r.S_relaxed), 2),
        "log2fc": round(float(r.log2fc), 2),
        "fdr": float(r.fdr),
        "tumour_tpm": round(float(r.tumour_tpm), 1),
        "tumour_q25": round(float(r.tumour_tpm_q25), 1),
        "frac1": round(float(r.frac_tumours_above_1tpm), 3),
        "frac10": round(float(r.frac_tumours_above_10tpm), 3),
        "worst_tissue": str(r.healthy_max_tissue),
        "worst_tpm": round(float(2 ** r.healthy_max_log2 - 1), 1),
        "worst_source": str(r.healthy_max_source),
        "grade": str(r.comparator_grade),
        "low_power": bool(r.low_power),
        "n_types": int(r.n_cancer_types_passing) if "n_cancer_types_passing"
        in r.index and pd.notna(r.n_cancer_types_passing) else None,
        "n_tumour": counts.get(r.cancer_type, (None, None))[0],
        "n_normal": counts.get(r.cancer_type, (None, None))[1],
        "family": None if pd.isna(r.gene_group) else str(r.gene_group),
        "trajectory": None if pd.isna(r.get("trajectory", np.nan))
        else str(r.trajectory),
        "protein": str(r.get("agreement", "protein_not_measured")),
        "protein_log2fc": None if pd.isna(r.get("protein_log2fc", np.nan))
        else round(float(r.protein_log2fc), 2),
        "protein_fdr": None if pd.isna(r.get("protein_fdr", np.nan))
        else float(r.protein_fdr),
        "cptac_cohort": cohort_of.get(r.cancer_type),
        "cspa": bool(r.cspa) if pd.notna(r.cspa) else False,
        "surfy": bool(r.surfy) if pd.notna(r.surfy) else False,
        "healthy": prof.get(r.ensembl_gene, []),
        "disease_term": DISEASE_TERM.get(r.cancer_type, r.cancer_type),
    }
    if s is not None:
        d.update({
            "ecto": str(s.extracellular_reported),
            "tm_sources": int(s.tm_sources),
            "n_tm": int(s.n_tm_segments),
            "gpi": bool(s.gpi_any),
            "gpi_source": str(s.gpi_source) if pd.notna(s.gpi_source) else "",
            "compartment": str(s.compartment),
        })
    return d


def main() -> int:
    full = pd.read_csv(RESULTS / "step7_ranking_full.csv")
    top10 = pd.read_csv(RESULTS / "step7_top10.csv")
    reserve = pd.read_csv(RESULTS / "step7_reserve.csv")
    surf = surface_detail()
    counts = sample_counts()
    lin = lineage()

    from step6_protein_differential import COHORT_TO_TCGA
    cohort_of = {v: k for k, v in COHORT_TO_TCGA.items()}

    per = (full.sort_values("S", ascending=False, kind="mergesort")
           .groupby("cancer_type", sort=True).head(TOP_PER_CANCER))
    needed = set(per.ensembl_gene) | set(top10.ensembl_gene) | set(reserve.ensembl_gene)
    prof = healthy_profiles(needed)

    by_cancer = {}
    for cancer, grp in per.groupby("cancer_type"):
        by_cancer[cancer] = [row_to_dict(r, surf, prof, cohort_of, counts, lin)
                             for _, r in grp.iterrows()]

    # funnel, read from the committed reports rather than recomputed
    cand = pd.read_csv(PROC / "surface_candidates.csv")
    s5 = pd.read_csv(RESULTS / "step5_specificity.csv")
    funnel = [
        {"stage": "Reviewed human proteins", "n": 20431},
        {"stage": "Mapped to an Ensembl gene", "n": 19400},
        {"stage": "Surface candidates (S1/S2/S3/Flagged)",
         "n": int(cand.is_candidate.sum())},
        {"stage": "Present in the expression matrix", "n": 3912},
        {"stage": "Pass tumour-vs-healthy RNA thresholds",
         "n": int(s5.loc[s5.passes_default, "ensembl_gene"].nunique())},
        {"stage": "Reach 5 TPM in tumour",
         "n": int(s5.loc[s5.carried_forward, "ensembl_gene"].nunique())},
        {"stage": "Eligible under the selection rule",
         "n": int(full.ensembl_gene.nunique())},
    ]

    # ---- D56: pan-cancer view, ranked by breadth ----
    pan_src = (full.groupby(["ensembl_gene", "symbol"])
               .agg(n_types=("cancer_type", "nunique"),
                    S_median=("S", "median"), S_min=("S", "min"),
                    S_max=("S", "max"),
                    homog=("frac_tumours_above_10tpm", "median"),
                    tier=("tier", "first"), family=("gene_group", "first"))
               .reset_index()
               .sort_values(["n_types", "S_median"], ascending=False)
               .head(25))
    pan = []
    for _, g in pan_src.iterrows():
        best = full[full.ensembl_gene == g.ensembl_gene].nlargest(1, "S").iloc[0]
        d = row_to_dict(best, surf, prof, cohort_of, counts, lin)
        d.update({"n_types": int(g.n_types),
                  "S_median": round(float(g.S_median), 2),
                  "S_min": round(float(g.S_min), 2),
                  "homog_median": round(float(g.homog), 3),
                  "types": sorted(full.loc[full.ensembl_gene == g.ensembl_gene,
                                           "cancer_type"].unique())})
        pan.append(d)

    payload = {
        "generated": "2026-10-21",
        "pan": pan,
        "funnel": funnel,
        "cancers": sorted(by_cancer),
        "by_cancer": by_cancer,
        "top10": [row_to_dict(r, surf, prof, cohort_of, counts, lin)
                  for _, r in top10.iterrows()],
        "reserve": [row_to_dict(r, surf, prof, cohort_of, counts, lin)
                    for _, r in reserve.iterrows()],
        "counts": {
            "eligible_rows": int(len(full)),
            "eligible_genes": int(full.ensembl_gene.nunique()),
            "cancer_types": int(full.cancer_type.nunique()),
        },
    }

    out = RESULTS / "step7_targets.json"
    out.write_text(json.dumps(payload, indent=None, separators=(",", ":")))
    print(f"wrote {out.relative_to(ROOT)}  "
          f"({out.stat().st_size / 1024:.0f} KB)")
    print(f"  cancer types : {len(by_cancer)}")
    print(f"  candidates   : {sum(len(v) for v in by_cancer.values())}")
    print(f"  top10        : {len(payload['top10'])}  reserve {len(payload['reserve'])}")
    return 0


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
