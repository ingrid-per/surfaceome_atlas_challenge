"""Step 12d: does a high S mean more per cell, or a rarer cell of origin?

For every eligible gene x cancer row, in log2(TPM+1):

    composite(sample) = median over that lineage's markers
    delta_tumour      = median(gene, tumour)    - median(composite, tumour)
    delta_healthy     = median(gene, home)      - median(composite, home)
    shift             = delta_tumour - delta_healthy

If a gene simply rides its lineage, both deltas are the same and shift is 0 —
the gap S measures is then arithmetic: the lineage is a smaller share of the
healthy tissue than of the tumour. A clearly positive shift means the tumour
carries more of the gene per cell of that lineage.

The gene under test is always removed from its own composite, otherwise a
marker gene would be compared against itself and shift to 0 by construction.

Nothing here feeds back into S or the ranking (D73). Cut-offs are reported at
0.5 / 1.0 / 1.5 log2 with 1.0 as the headline (D80); 1.0 log2 is twofold.

Writes results/step12_lineage_ratio.csv, step12_enrichment.txt,
step12_summary.txt.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step12_markers import LINEAGES

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

CUTOFFS = [0.5, 1.0, 1.5]       # D80
HEADLINE = 1.0
MIN_MARKERS_AFTER_SELF = 3      # a composite needs 3 markers once the gene is out
TOP_N = 50                      # enrichment: top 50 by S
FLOOR_TPM1 = 1.0                # log2(TPM+1) below this is at the expression floor
BIG_COMPOSITE_SHIFT = 2.0       # |composite difference| beyond this is flagged


def load():
    co = pd.read_csv(RESULTS / "step12_cell_of_origin.csv")
    st = pd.read_csv(RESULTS / "step12_lineage_status.csv")
    mk = pd.read_csv(RESULTS / "step12_markers_used.csv")
    elig = pd.read_csv(RESULTS / "step7_ranking_full.csv")
    marker_mat = pd.read_parquet(PROC / "step12_marker_tpm.parquet")
    cand_mat = pd.read_parquet(PROC / "step4_candidate_tpm.parquet")
    samp = pd.read_parquet(PROC / "step4_samples.parquet")
    return co, st, mk, elig, marker_mat, cand_mat, samp


def main() -> int:
    co, st, mk, elig, marker_mat, cand_mat, samp = load()

    co["lineage_id"] = co.cancer_type.map(lambda c: LINEAGES[c][0] or "none")
    status = st.set_index("lineage_id")
    co = co.join(status[["status", "thin_markers", "reason"]], on="lineage_id")

    tum = samp[(samp.study == "TCGA") & (samp.sample_type == "Primary Tumor")]
    gtx = samp[(samp.study == "GTEX") & (samp.sample_type == "Normal Tissue")]
    tum_cols = {c: [s for s in g["sample"] if s in marker_mat.columns]
                for c, g in tum.groupby("disease")}
    home_cols = {t: [s for s in g["sample"] if s in marker_mat.columns]
                 for t, g in gtx.groupby("disease")}

    markers_by_lineage = {lid: [e for e in g.ensembl_gene.dropna()
                                if e in marker_mat.index]
                          for lid, g in mk[mk.used].groupby("lineage_id")}

    def gene_vector(ens: str, cols: list[str]) -> np.ndarray | None:
        """A gene's log2(TPM+1) over samples, from whichever matrix has it."""
        if ens in cand_mat.index:
            return cand_mat.loc[ens, cols].to_numpy(dtype=np.float32)
        if ens in marker_mat.index:
            return marker_mat.loc[ens, cols].to_numpy(dtype=np.float32)
        return None

    def composite(lid: str, cols: list[str], drop: str | None) -> tuple[float, int]:
        """Median across samples of the per-sample median across markers."""
        genes = [g for g in markers_by_lineage.get(lid, []) if g != drop]
        if len(genes) < MIN_MARKERS_AFTER_SELF or not cols:
            return float("nan"), len(genes)
        sub = marker_mat.loc[genes, cols].to_numpy(dtype=np.float32)
        per_sample = np.median(sub, axis=0)
        return float(np.median(per_sample)), len(genes)

    # The composite only changes when the gene under test is itself a marker,
    # which is rare, so compute it once per lineage and patch those few rows.
    base_cache: dict[tuple, tuple[float, int]] = {}

    def cached(lid, cols, drop, key):
        if drop is not None and drop in markers_by_lineage.get(lid, []):
            return composite(lid, cols, drop)
        if key not in base_cache:
            base_cache[key] = composite(lid, cols, None)
        return base_cache[key]

    rows = []
    co_by_ct = co.set_index("cancer_type")
    for r in elig.itertuples():
        info = co_by_ct.loc[r.cancer_type]
        lid = info.lineage_id
        out = {"ensembl_gene": r.ensembl_gene, "symbol": r.symbol,
               "cancer_type": r.cancer_type, "tier": r.tier, "S": r.S,
               "log2fc": r.log2fc, "tumour_tpm": r.tumour_tpm,
               "cell_of_origin": info.cell_of_origin, "lineage_id": lid,
               "home_tissue": info.home_tissue,
               "confidence": info.confidence,
               "lineage_share_expected": info.lineage_share_of_home,
               "thin_markers": bool(info.thin_markers)}

        if info.status != "scored":
            rows.append(out | {"scored": False,
                               "unscored_reason": f"lineage: {info.reason}"})
            continue
        if not isinstance(info.home_tissue, str) or not info.home_tissue:
            rows.append(out | {"scored": False,
                               "unscored_reason": "no home tissue in GTEx (D73/D74)"})
            continue

        tcols, hcols = tum_cols.get(r.cancer_type, []), home_cols.get(info.home_tissue, [])
        gt = gene_vector(r.ensembl_gene, tcols) if tcols else None
        gh = gene_vector(r.ensembl_gene, hcols) if hcols else None
        if gt is None or gh is None:
            rows.append(out | {"scored": False,
                               "unscored_reason": "gene absent from the matrix"})
            continue

        ct, n_used = cached(lid, tcols, r.ensembl_gene, (lid, "t", r.cancer_type))
        ch, _ = cached(lid, hcols, r.ensembl_gene, (lid, "h", info.home_tissue))
        if not np.isfinite(ct) or not np.isfinite(ch):
            rows.append(out | {"scored": False,
                               "unscored_reason": f"only {n_used} markers after "
                                                  "removing the gene itself"})
            continue

        dt = float(np.median(gt)) - ct
        dh = float(np.median(gh)) - ch
        rows.append(out | {
            "scored": True, "unscored_reason": "",
            "n_markers_used": n_used,
            "composite_tumour": round(ct, 4), "composite_home": round(ch, 4),
            "composite_diff": round(ct - ch, 4),
            "delta_tumour": round(dt, 4), "delta_healthy": round(dh, 4),
            "shift": round(dt - dh, 4),
            "home_composite_at_floor": ch < FLOOR_TPM1,
            "composite_diff_large": abs(ct - ch) > BIG_COMPOSITE_SHIFT,
        })

    df = pd.DataFrame(rows)
    # NOTE: the column is called "shift" per the spec, and DataFrame.shift is
    # a method — every access below must be by brackets, never by attribute.
    for c in CUTOFFS:
        df[f"tracks_lineage_{c}"] = (df["shift"].abs() < c
                                     if "shift" in df.columns else False)
    df.to_csv(RESULTS / "step12_lineage_ratio.csv", index=False)

    # ---------------------------------------------------------- enrichment --
    scored = df[df.scored].copy()
    clear = scored[scored.confidence == "clear"]
    well = clear[~clear.thin_markers]
    thin = clear[clear.thin_markers]

    def compare(d: pd.DataFrame, label: str) -> list[str]:
        if len(d) < TOP_N + 10:
            return [f"{label}: only {len(d)} scored rows — not enough to test", ""]
        d = d.sort_values("S", ascending=False, kind="mergesort")
        top, rest = d.head(TOP_N), d.iloc[TOP_N:]
        out = [f"{label}  (n={len(d):,}; top {TOP_N} vs the other {len(rest):,})", ""]
        out.append(f"  median shift, top {TOP_N}      : {top['shift'].median():+.2f} log2")
        out.append(f"  median shift, the rest      : {rest['shift'].median():+.2f} log2")
        u, p = mannwhitneyu(top["shift"], rest["shift"], alternative="two-sided")
        out.append(f"  Mann-Whitney                : p = {p:.3g}")
        out.append("")
        for c in CUTOFFS:
            a = float((top["shift"].abs() < c).mean())
            b = float((rest["shift"].abs() < c).mean())
            ratio = (a / b) if b > 0 else float("inf")
            out.append(f"  |shift| < {c:<4}  top {TOP_N}: {a:5.1%}   rest: {b:5.1%}"
                       f"   ratio {ratio:.2f}x")
        out.append("")
        return out

    # --------------------------------------------------------- diagnostic --
    # Does `shift` measure the GENE, or just the lineage? If the composite
    # moves a lot between tumour and healthy, shift moves with it by
    # construction. This correlation decides whether the enrichment below
    # means anything.
    diag = well.dropna(subset=["composite_diff", "shift"])
    rho_cs, p_cs = spearmanr(diag.composite_diff, diag["shift"])
    rho_sS, p_sS = spearmanr(diag.S, diag["shift"])

    lines = ["STEP 12 — lineage dilution, enrichment at the top of the ranking",
             "Generated by src/step12_lineage_ratio.py — do not edit by hand.",
             "",
             "!! READ THE DIAGNOSTIC FIRST. The enrichment numbers below are",
             "!! reported because they were computed, not because they are",
             "!! trustworthy.",
             "",
             "DIAGNOSTIC — is shift a property of the gene or of the lineage?",
             "",
             f"  shift vs composite difference : rho {rho_cs:+.3f}  (p {p_cs:.2g})",
             f"  shift vs S                    : rho {rho_sS:+.3f}  (p {p_sS:.2g})",
             "",
             "  shift is almost entirely determined by how far the lineage",
             "  composite moves between tumour and healthy tissue, and barely",
             "  related to S at all. Per-lineage medians:",
             ""]
    gl = (diag.groupby("lineage_id")
          .agg(n=("shift", "size"), comp=("composite_diff", "median"),
               sh=("shift", "median"), exp=("lineage_share_expected", "first"))
          .sort_values("comp"))
    lines.append(f"  {'lineage':<20}{'n':>5}{'composite diff':>16}{'median shift':>14}"
                 f"  expected share")
    for lid, r in gl.iterrows():
        lines.append(f"  {lid:<20}{int(r.n):>5}{r.comp:>+16.2f}{r.sh:>+14.2f}"
                     f"  {r.exp}")
    lines += [
        "",
        "  Every lineage's median shift is close to the negative of its",
        "  composite difference. That is arithmetic, not biology: the gene term",
        "  is swamped by the denominator. Two things move the denominator and",
        "  only one of them is cell fraction —",
        "",
        "    1. the lineage really is a different share of the two tissues;",
        "    2. tumours DEDIFFERENTIATE, so they lose the very markers that",
        "       define the lineage. Colonocyte markers (ALPI, FABP1, FABP2)",
        "       are LOWER in colon tumour than in healthy colon, which reads",
        "       as 'fewer colonocytes' when it is really 'less differentiated",
        "       colonocytes'.",
        "",
        "  The positive controls fail accordingly: TYR, held out of its own",
        "  composite, should track the melanocyte lineage and instead scores",
        "  +2.52 in uveal melanoma. See results/step12_controls.csv.",
        "",
        "  CONCLUSION: this implementation does not separate dilution from",
        "  per-cell expression. The numbers below stand as a record of what was",
        "  run; they are not a finding.",
        "",
             "shift = (gene - lineage composite) in tumour, minus the same in the",
             "lineage's healthy home tissue. Near zero means the gene tracks its",
             "lineage, so a high S reflects the cell of origin being rarer in the",
             "healthy tissue rather than more protein per cell.",
             "1.0 log2 means twofold.",
             ""]
    lines += compare(well, "WELL-SUPPORTED LINEAGES (headline)")
    lines += compare(thin, "THIN-MARKER LINEAGES (reported separately, D80)")

    head = well.sort_values("S", ascending=False).head(TOP_N)
    lines += ["", f"the top {TOP_N} by S in well-supported lineages:", "",
              f"  {'gene':<10} {'cancer type':<34} {'S':>5} {'shift':>7}  verdict"]
    for _, r in head.iterrows():
        v = ("tracks its lineage" if abs(r["shift"]) < HEADLINE
             else "more per cell" if r["shift"] > 0 else "less per cell")
        lines.append(f"  {str(r.symbol):<10} {r.cancer_type[:34]:<34} "
                     f"{r.S:>5.2f} {r['shift']:>+7.2f}  {v}")

    (RESULTS / "step12_enrichment.txt").write_text("\n".join(lines) + "\n")

    # ------------------------------------------------------------ controls --
    # 12.1 positive: a canonical melanocyte gene should TRACK the melanocyte
    #      composite. TYR is itself a marker, so it is held out of its own
    #      composite; SOX10 is canonical but was never selected, so it is an
    #      independent positive.
    # 12.2 negative: genes of another lineage must NOT track melanocytes.
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)
    sym2ens = (gmap.dropna(subset=["ensembl_gene"])
               .drop_duplicates("symbol").set_index("symbol").ensembl_gene)
    CONTROLS = [("TYR", "positive, held out of its own composite"),
                ("SOX10", "positive, canonical but never selected"),
                ("MLANA", "positive, held out of its own composite"),
                ("KRT5", "negative, keratinocyte"),
                ("ALB", "negative, hepatocyte")]
    crows = []
    for sym, kind in CONTROLS:
        ens = sym2ens.get(sym)
        for ct in ["Uveal Melanoma", "Skin Cutaneous Melanoma"]:
            home = co_by_ct.loc[ct].home_tissue
            tcols, hcols = tum_cols.get(ct, []), home_cols.get(home, [])
            gt = gene_vector(ens, tcols) if isinstance(ens, str) and tcols else None
            gh = gene_vector(ens, hcols) if isinstance(ens, str) and hcols else None
            if gt is None or gh is None:
                crows.append({"gene": sym, "kind": kind, "cancer_type": ct,
                              "shift": float("nan"), "note": "not in the matrix"})
                continue
            ctum, n_used = composite("melanocyte", tcols, ens)
            chome, _ = composite("melanocyte", hcols, ens)
            sh = (float(np.median(gt)) - ctum) - (float(np.median(gh)) - chome)
            crows.append({"gene": sym, "kind": kind, "cancer_type": ct,
                          "n_markers_used": n_used,
                          "median_tumour": round(float(np.median(gt)), 3),
                          "median_home": round(float(np.median(gh)), 3),
                          "shift": round(sh, 3),
                          "tracks_lineage": abs(sh) < HEADLINE})
    cdf = pd.DataFrame(crows)
    cdf.to_csv(RESULTS / "step12_controls.csv", index=False)

    # ------------------------------------------------------------- summary --
    grade = (elig.drop_duplicates("cancer_type")
             .set_index("cancer_type").comparator_grade)
    co2 = co.set_index("cancer_type")
    no_multi = set(st.loc[st.reason.str.contains("2 or more papers", na=False),
                          "lineage_id"])
    ct_no_multi = sorted(co2.index[co2.lineage_id.isin(no_multi)])
    ct_weak_comp = sorted(grade.index[grade != "direct"])

    both = sorted(set(ct_no_multi) & set(ct_weak_comp))
    sm = [
        "STEP 12 — summary",
        "Generated by src/step12_lineage_ratio.py — do not edit by hand.",
        "",
        "WHAT WAS MEASURED",
        "",
        "  For each eligible gene and cancer type, how far the gene sits above",
        "  its lineage's marker genes in the tumour, minus the same distance in",
        "  the healthy tissue where that lineage lives. A shift near zero means",
        "  the gene simply rides its lineage: the tumour looks high because the",
        "  cell of origin is a larger share of the tumour than of the healthy",
        "  tissue, not because each cell carries more. 1.0 log2 means twofold.",
        "",
        "COVERAGE",
        "",
        f"  eligible gene x cancer rows      : {len(df):,}",
        f"  scored                           : {int(df.scored.sum()):,}",
        f"  not scored                       : {int((~df.scored).sum()):,}",
        "",
    ]
    for reason, n in df[~df.scored].unscored_reason.value_counts().items():
        sm.append(f"    {reason:<52} {n:>5,}")
    sm += [
        "",
        f"  lineages defined                 : {len(st)}",
        f"    scored                         : {int((st.status=='scored').sum())}"
        f"  (thin markers: {int(st.thin_markers.sum())})",
        f"    no multi-paper markers (D77)   : {len(no_multi)}",
        "",
        "MARKER COVERAGE AGAINST COMPARATOR QUALITY",
        "",
        "  A fair question is whether the lineages the marker database covers",
        "  badly are the same ones our healthy comparator covers badly. If they",
        "  were, Step 12 would be blind exactly where Step 4 is weakest.",
        "",
        f"  cancer types with no multi-paper markers ({len(ct_no_multi)}):",
    ]
    for c in ct_no_multi:
        sm.append(f"    {c:<42} comparator: {grade.get(c, '—')}")
    sm += ["",
           f"  cancer types without a direct comparator ({len(ct_weak_comp)}):"]
    for c in ct_weak_comp:
        sm.append(f"    {c:<42} grade: {grade.get(c)}")
    sm += ["",
           f"  in both groups ({len(both)}): {', '.join(both) if both else 'none'}",
           ""]
    sm.append("  They are largely NOT the same lineages." if len(both) <= 1 else
              "  They overlap substantially — read both limitations together.")
    sm += [
        "",
        "CONTROLS",
        "",
        f"  {'gene':<7}{'cancer type':<28}{'shift':>7}  {'tracks':<7} {'expectation'}",
    ]
    for _, r in cdf.iterrows():
        sm.append(f"  {r.gene:<7}{r.cancer_type[:27]:<28}{r['shift']:>+7.2f}  "
                  f"{str(r.get('tracks_lineage','—')):<7} {r.kind}")
    sm += [
        "",
        "KNOWN BLIND SPOTS",
        "",
        "  Pheochromocytoma and thyroid carcinoma cannot be scored: chromaffin",
        "  cells have one CellMarker gene (NCAM1, pan-neural) and thyroid",
        "  follicular cells have two, both transcription factors that are low",
        "  and noisy in bulk RNA. DBH and DCSTAMP are therefore untested here.",
        "",
        "  DCSTAMP may not belong to the thyroid follicular lineage at all.",
        "  DC-STAMP is a myeloid gene, so the confound for that candidate is",
        "  immune infiltrate, not lineage dilution — the limitation D61 states",
        "  and deliberately does not compute.",
        "",
        "  Mesothelioma, ovarian and DLBCL have no healthy tissue that",
        "  represents their lineage (D73, D74), so UPK3B is also untested.",
    ]
    (RESULTS / "step12_summary.txt").write_text("\n".join(sm) + "\n")

    print("\n".join(lines))
    print("\n" + "\n".join(sm[sm.index("CONTROLS"):]))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
