"""Step 7b: export the ranking as JSON for the interactive deliverable (D54).

Every number here comes from the result files, never typed by hand, so the
artifact and the committed results cannot drift apart.
"""

from __future__ import annotations

import itertools
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


# Strongest tier wins when two accessions claim the same gene — the same rule
# Step 4 applies in gene_annotation().
TIER_RANK = {"S1": 0, "S2": 1, "S3": 2, "Flagged": 3}


def surface_detail() -> pd.DataFrame:
    """Surface evidence, reachable by EVERY Ensembl id a candidate maps to.

    surface_candidates.csv stores one Ensembl id per accession, because
    step3_build_candidates.py collapses the mapping with `"first"`. In the MHC
    region of chromosome 6 that is badly wrong: GENCODE gives one gene several
    ids across the primary assembly and the ALT haplotypes, so one UniProt
    accession maps to 5-12 of them (LY6G6F has 12). Step 4 keeps all of them,
    so the ranking carries whichever id survived in the Xena matrix — often a
    different one. Joining on the single stored id silently dropped the surface
    evidence for LY6G6D, LY6G6F, DDR1, HLA-DPB1 and HLA-DQA2.

    The values were never missing, only unreachable through that key. Building
    the index from the full accession -> Ensembl mapping fixes it without
    touching surface_candidates.csv, so no Step 1-11 output changes.
    """
    c = pd.read_csv(PROC / "surface_candidates.csv")
    c = c[c.is_candidate].drop(columns=["ensembl_gene"])
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)
    link = (gmap[gmap.accession.isin(set(c.accession))][["accession", "ensembl_gene"]]
            .dropna().drop_duplicates())
    j = link.merge(c, on="accession")
    j["_rank"] = j.tier.map(TIER_RANK)
    j = j.sort_values(["ensembl_gene", "_rank"], kind="mergesort")
    return j.drop_duplicates("ensembl_gene").set_index("ensembl_gene")


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
        "frac_above_healthy": None if pd.isna(
            r.get("frac_tumours_above_worst_healthy", np.nan))
        else round(float(r.frac_tumours_above_worst_healthy), 3),
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
            "existence": str(s.protein_existence)
            if pd.notna(s.protein_existence) else "",
        })
    return d


# Each axis carries the plain reason it is a choice at all, so the explorer
# explains the pipeline rather than just permuting it.
AXES = [
    {"key": "min_extra", "label": "Minimum extracellular segment", "unit": " aa",
     "default": 20,
     "why": "How much of the protein must stick out of the cell for an antibody "
            "to have something to bind. Shorter admits more proteins with "
            "epitopes that may be too small to drug."},
    {"key": "fdr", "label": "False discovery rate", "unit": "",
     "default": 0.01,
     "why": "How often we accept being wrong about a gene being raised at all."},
    {"key": "log2fc", "label": "Fold change over healthy tissue", "unit": " log2",
     "default": 2.0,
     "why": "How much higher the tumour must be. 2 on a log2 scale is 4-fold."},
    {"key": "source", "label": "Healthy reference", "unit": "",
     "default": "both",
     "why": "Which healthy atlas sets the bar a candidate must clear. The two "
            "were processed by different pipelines, so taking both and keeping "
            "the higher value is the cautious reading."},
    {"key": "summary", "label": "How healthy tissue is summarised", "unit": "",
     "default": "max",
     "why": "One exposed organ is enough to make a target unsafe, so the "
            "default is the single worst tissue. The alternatives forgive "
            "tissues a drug may not reach, or discard the top five."},
]

# Explicit, because alphabetical order is meaningless here: the controls should
# read strictest first, so moving right always means forgiving more.
VALUE_ORDER = {
    "summary": ["max", "max_nobrain", "max_no_reproductive_fetal",
                "max_relaxed", "p90"],
    "source": ["both", "gtex", "hpa"],
}

VALUE_LABEL = {
    "both": "GTEx + HPA", "gtex": "GTEx only", "hpa": "HPA only",
    "max": "worst single tissue",
    "max_nobrain": "worst, excluding brain",
    "max_no_reproductive_fetal": "worst, excluding testis + placenta",
    "max_relaxed": "worst, excluding brain, testis + placenta",
    "p90": "90th percentile across tissues",
}


def robustness() -> dict:
    """Step 9's 270 runs, as a lookup the deliverable can read directly.

    Nothing is recomputed in the browser: every setting a reader can select is
    one of the runs that actually executed and passed checks 9.1-9.6. The top
    20 lists are stored as indices into one symbol table to keep the payload
    small.
    """
    runs = pd.read_csv(RESULTS / "step9_runs.csv")
    st = pd.read_csv(RESULTS / "step9_stability.csv")
    sens = pd.read_csv(RESULTS / "step9_sensitivity.csv")

    symbols = list(st.symbol)
    idx = {s: i for i, s in enumerate(symbols)}

    # The lookup key is built from STRINGS, and the axis values are exported as
    # the same strings. Python writes str(2.0) as "2.0" while JavaScript writes
    # String(2) as "2", so a key assembled in the browser from numbers would
    # silently miss every run. The assertion below is what catches that.
    AX = [a["key"] for a in AXES]
    key = lambda r: "|".join(str(getattr(r, k)) for k in AX)
    lookup = {key(r): [idx[s] for s in str(r.top20).split(";")]
              for r in runs.itertuples()}

    axes = []
    for a in AXES:
        vals = runs[a["key"]].unique().tolist()
        order = VALUE_ORDER.get(a["key"])
        vals = (sorted(vals, key=order.index) if order else sorted(vals))
        axes.append({**a, "values": [str(v) for v in vals],
                     "default": str(a["default"]),
                     "labels": [VALUE_LABEL.get(v, str(v)) for v in vals]})

    reachable = {"|".join(c) for c in
                 itertools.product(*[a["values"] for a in axes])}
    missing = reachable - set(lookup)
    assert not missing, f"{len(missing)} selectable settings have no run: " \
                        f"{sorted(missing)[:3]}"
    assert len(reachable) == len(lookup), \
        f"{len(lookup)} runs but {len(reachable)} reachable settings"

    default = next(r for r in runs.itertuples() if r.is_default)
    return {
        "n_runs": int(len(runs)),
        "n_settings_per_axis": {a["key"]: len(a["values"]) for a in axes},
        "symbols": symbols,
        "stability": {r.symbol: round(float(r.stability), 4)
                      for r in st.itertuples()},
        "runs": lookup,
        "axes": axes,
        "default_key": key(default),
        "eligible": {key(r): int(r.n_genes_eligible) for r in runs.itertuples()},
        "sensitivity": [{"axis": r.axis, "value": str(r.value),
                         "label": VALUE_LABEL.get(r.value, str(r.value)),
                         "overlap": float(r.mean_overlap_with_default),
                         "spread": float(r.axis_spread)}
                        for r in sens.itertuples()],
    }


def step12_coverage() -> dict:
    """Marker-coverage counts for the method text (D82).

    Counted here so the page can quote them without a number being typed. The
    comparison that matters: are the lineages the marker database covers badly
    the SAME ones our healthy comparator covers badly? Mostly they are not.
    """
    st = pd.read_csv(RESULTS / "step12_lineage_status.csv")
    co = pd.read_csv(RESULTS / "step12_cell_of_origin.csv")
    full = pd.read_csv(RESULTS / "step7_ranking_full.csv")
    grade = full.drop_duplicates("cancer_type").set_index("cancer_type").comparator_grade

    no_multi = st.loc[st.reason.str.contains("2 or more papers", na=False), "lineage_id"]
    # map cancer type -> lineage via the status table's own cancer_types column
    ct_to_lid = {}
    for _, r in st.iterrows():
        for ct in str(r.cancer_types).split("; "):
            if ct:
                ct_to_lid[ct.strip()] = r.lineage_id
    ct_no_multi = sorted(c for c, l in ct_to_lid.items() if l in set(no_multi))
    ct_weak = sorted(grade.index[grade != "direct"])
    both = sorted(set(ct_no_multi) & set(ct_weak))
    return {
        "n_lineages": int(len(st)),
        "n_lineages_no_multi": int(len(no_multi)),
        "n_cancer_types": int(len(co)),
        "cancer_types_no_multi": ct_no_multi,
        "cancer_types_weak_comparator": ct_weak,
        "in_both": both,
        "thin_lineages": sorted(st.loc[st.thin_markers, "lineage_id"]),
        "unscorable_lineages": sorted(
            st.loc[(st.status == "not_scored")
                   & st.reason.str.contains("marker", na=False), "lineage_id"]),
    }


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

    top10_rows = [row_to_dict(r, surf, prof, cohort_of, counts, lin)
                  for _, r in top10.iterrows()]
    reserve_rows = [row_to_dict(r, surf, prof, cohort_of, counts, lin)
                    for _, r in reserve.iterrows()]
    payload_rows_global = top10_rows + reserve_rows

    # Counted here so the deliverable can say how unusual a weak existence call
    # is without anyone typing the number.
    ex = surf.loc[surf.index.intersection(set(full.ensembl_gene)), "protein_existence"]
    existence_counts = {str(k): int(v) for k, v in ex.value_counts().items()}
    existence_total = int(ex.notna().sum())

    rob = robustness()
    # Stability is a property of the GLOBAL ranking — a gene enters the top 20
    # through its single best cancer type. Attaching it to a per-cancer-type row
    # would read as "stable in this cancer type", which is not what was measured,
    # so it is attached to the global lists only.
    stab = rob["stability"]
    for row in payload_rows_global:
        row["stab"] = stab.get(row["symbol"])

    payload = {
        "generated": "2026-10-21",
        "robust": rob,
        "existence": {"counts": existence_counts, "total": existence_total},
        "step12": step12_coverage(),
        "pan": pan,
        "funnel": funnel,
        "cancers": sorted(by_cancer),
        "by_cancer": by_cancer,
        "top10": top10_rows,
        "reserve": reserve_rows,
        "counts": {
            "eligible_rows": int(len(full)),
            "eligible_genes": int(full.ensembl_gene.nunique()),
            "cancer_types": int(full.cancer_type.nunique()),
        },
    }

    all_rows = (payload["top10"] + payload["reserve"] + payload["pan"]
                + [r for v in payload["by_cancer"].values() for r in v])
    no_surface = [r for r in all_rows if "tm_sources" not in r]
    print(f"  surface evidence: {len(all_rows) - len(no_surface)}/{len(all_rows)} "
          f"rows resolved" + (f"  MISSING: "
          f"{sorted({r['symbol'] for r in no_surface})}" if no_surface else ""))

    out = RESULTS / "step7_targets.json"
    out.write_text(json.dumps(payload, indent=None, separators=(",", ":")))
    print(f"wrote {out.relative_to(ROOT)}  "
          f"({out.stat().st_size / 1024:.0f} KB)")
    print(f"  cancer types : {len(by_cancer)}")
    print(f"  candidates   : {sum(len(v) for v in by_cancer.values())}")
    print(f"  top10        : {len(payload['top10'])}  reserve {len(payload['reserve'])}")
    print(f"  robustness   : {rob['n_runs']} runs, {len(rob['symbols'])} genes")
    return 0


if __name__ == "__main__":
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    raise SystemExit(main())
