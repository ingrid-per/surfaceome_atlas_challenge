"""Step 12b: pick lineage marker genes from CellMarker 2.0 (Gate 2).

Implements D76-D79.

  rank candidates by the number of DISTINCT PMIDs citing them, take the top 8
  fewer than 3 markers            -> lineage not scored
  zero markers cited by 2+ papers -> lineage not scored (D77)
  fewer than 3 such markers       -> scored, carries thin_markers

Paper count is the only quality signal the database offers that does not just
reflect how large one study was. Gastric epithelium is the cautionary case: 75
candidate genes, none cited twice, topped by ACTC1 and SBSPON, which are not
gastric markers. A composite built on that is a meaningless denominator, and a
meaningless denominator still produces a number.

CellMarker's own vocabulary does not match the names a pathologist would use,
and some near-matches are traps — "Follicular cell" sits beside T and B
follicular immune cells. So each lineage below is defined by an explicit list
of accepted `cell_name` values plus a tissue restriction, never by a
substring match.

The candidate gene under test is NOT excluded here; that is done per row at
scoring time, because the same marker set serves every gene.

Writes results/step12_markers_used.csv and results/step12_lineage_status.csv.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

CELLMARKER = RAW / "CellMarker2_Cell_marker_Human.xlsx"
TOP_N = 8          # D76
MIN_MARKERS = 3    # D76
MIN_MULTI = 1      # D77: at least one marker cited by 2+ papers

# cancer type -> (lineage id, accepted CellMarker cell_name values, tissue filter)
#
# Keyed by cancer type rather than by cell-of-origin text because the same
# words mean different lineages in different organs: "luminal epithelial cell"
# is breast in one row and prostate in another, and they need different tissue
# restrictions. Cancer types sharing a lineage share its id, so they share one
# marker set.
LINEAGES: dict[str, tuple[str, list[str], str | None]] = {
    "Adrenocortical Cancer": ("adrenal_cortical", ["Definitive zone cell"], "Adrenal"),
    "Brain Lower Grade Glioma": ("astrocyte", ["Astrocyte"], None),
    "Glioblastoma Multiforme": ("astrocyte", ["Astrocyte"], None),
    "Breast Invasive Carcinoma": ("luminal_breast",
        ["Luminal cell", "Luminal epithelial cell", "Mature luminal cell"], "Breast"),
    "Prostate Adenocarcinoma": ("luminal_prostate",
        ["Luminal cell", "Luminal epithelial cell"], "Prostate"),
    "Cervical & Endocervical Cancer": ("squamous",
        ["Squamous epithelial cell", "Squamous cell"], None),
    "Esophageal Carcinoma": ("squamous",
        ["Squamous epithelial cell", "Squamous cell"], None),
    "Head & Neck Squamous Cell Carcinoma": ("squamous",
        ["Squamous epithelial cell", "Squamous cell"], None),
    "Cholangiocarcinoma": ("cholangiocyte", ["Cholangiocyte"], None),
    "Colon Adenocarcinoma": ("colonocyte", ["Enterocyte"], "Intestine|Colon"),
    "Rectum Adenocarcinoma": ("colonocyte", ["Enterocyte"], "Intestine|Colon"),
    "Kidney Clear Cell Carcinoma": ("proximal_tubule",
        ["Proximal tubular cell", "Pan Proximal tubular cell",
         "Proximal convoluted tubular cell", "Proximal straight tubular cell",
         "Proximal tubule cell"], "Kidney"),
    "Kidney Papillary Cell Carcinoma": ("proximal_tubule",
        ["Proximal tubular cell", "Pan Proximal tubular cell",
         "Proximal convoluted tubular cell", "Proximal straight tubular cell",
         "Proximal tubule cell"], "Kidney"),
    "Kidney Chromophobe": ("intercalated",
        ["Intercalated cell", "Pan-Collecting duct-intercalated cell",
         "Collecting duct-intercalated cell (type A)",
         "Collecting duct-intercalated cell (type B)",
         "Type A intercalated cell"], "Kidney"),
    "Liver Hepatocellular Carcinoma": ("hepatocyte", ["Hepatocyte"], "Liver"),
    "Lung Adenocarcinoma": ("alveolar_type2",
        ["Alveolar type II (ATII) cell", "Type II pneumocyte",
         "Alveolar cell Type 2", "Alveolar pneumocyte Type II",
         "Pulmonary alveolar cell Type 2"], "Lung"),
    "Lung Squamous Cell Carcinoma": ("bronchial_basal", ["Basal cell"], "Lung"),
    "Mesothelioma": ("mesothelial", ["Mesothelial cell"], None),
    "Ovarian Serous Cystadenocarcinoma": ("tubal_secretory",
        ["Secretory cell"], "Fallopian|Oviduct|Ovary"),
    "Pancreatic Adenocarcinoma": ("pancreatic_ductal", ["Ductal cell"], "Pancreas"),
    "Pheochromocytoma & Paraganglioma": ("chromaffin", ["Chromaffin cell"], "Adrenal"),
    "Skin Cutaneous Melanoma": ("melanocyte", ["Melanocyte"], None),
    "Uveal Melanoma": ("melanocyte", ["Melanocyte"], None),
    "Stomach Adenocarcinoma": ("gastric_epithelial",
        ["Pit cell", "Pit progenitor cell", "Parietal cell", "Chief cell",
         "Zymogenic chief cell", "Epithelial cell"], "Stomach"),
    "Testicular Germ Cell Tumor": ("germ_cell", ["Germ cell"], "Gonad|Testis"),
    "Thyroid Carcinoma": ("thyroid_follicular", ["Follicular cell"], "Thyroid"),
    "Uterine Corpus Endometrioid Carcinoma": ("endometrial_epithelial",
        ["Epithelial cell", "Ciliated cell"], "Uterus|Endometrium"),
    "Uterine Carcinosarcoma": ("endometrial_epithelial",
        ["Epithelial cell", "Ciliated cell"], "Uterus|Endometrium"),
    "Diffuse Large B-Cell Lymphoma": ("b_cell", ["B cell"], None),
    "Sarcoma": ("", [], None),          # D73: no single lineage applies
}

# Lineages accepted despite thin citation support, with the reason (D78).
OVERRIDES = {
    "germ_cell": "canonical germ-cell markers (SYCP3, DDX4, DMRT1); single "
                 "citations reflect CellMarker's curation depth, not the biology",
}


def load_cellmarker() -> pd.DataFrame:
    if not CELLMARKER.exists():
        raise SystemExit(f"{CELLMARKER} not found — run "
                         "src/step1_download.py --only M1 --download")
    df = pd.read_excel(CELLMARKER)
    # Normal cells only: we are measuring how much of a HEALTHY lineage is
    # present, so markers defined inside tumours would be circular.
    df = df[(df.cell_type == "Normal cell") & df.Symbol.notna()].copy()
    df["Symbol"] = df.Symbol.astype(str).str.strip()
    return df


def main() -> int:
    cm = load_cellmarker()
    co = pd.read_csv(RESULTS / "step12_cell_of_origin.csv")
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)
    sym2ens = (gmap.dropna(subset=["ensembl_gene"])
               .drop_duplicates("symbol").set_index("symbol").ensembl_gene)

    missing = set(co.cancer_type) - set(LINEAGES)
    assert not missing, f"cancer types with no lineage definition: {sorted(missing)}"

    marker_rows, status_rows, seen = [], [], {}
    for _, row in co.iterrows():
        lid, names, tissue = LINEAGES[row.cancer_type]
        if lid in seen:                      # shared lineage, already resolved
            seen[lid]["cancer_types"].append(row.cancer_type)
            continue

        if not lid or not names:
            lid = lid or "none"          # sarcoma: keep one consistent key
            status_rows.append({"lineage_id": lid,
                                "cell_of_origin": row.cell_of_origin,
                                "n_markers": 0, "n_multi_paper": 0,
                                "status": "not_scored", "thin_markers": False,
                                "reason": "no single lineage applies"})
            seen[lid] = status_rows[-1] | {"cancer_types": [row.cancer_type]}
            continue

        sub = cm[cm.cell_name.isin(names)]
        if tissue:
            sub = sub[sub.tissue_class.str.contains(tissue, case=False, na=False)]

        support = (sub.groupby("Symbol").PMID.nunique()
                   .sort_values(ascending=False, kind="mergesort"))
        # Deterministic order: papers descending, then symbol, so two runs agree.
        support = support.reset_index().sort_values(
            ["PMID", "Symbol"], ascending=[False, True]).set_index("Symbol").PMID
        picked = support.head(TOP_N)
        n_multi = int((picked >= 2).sum())

        if len(picked) < MIN_MARKERS:
            status, thin, reason = ("not_scored", False,
                                    f"only {len(picked)} marker(s) in CellMarker")
        elif n_multi < MIN_MULTI:
            status, thin, reason = ("not_scored", False,
                                    "no marker cited by 2 or more papers (D77)")
        else:
            thin = n_multi < 3
            status = "scored"
            reason = (f"override: {OVERRIDES[lid]}" if thin and lid in OVERRIDES
                      else (f"only {n_multi} of {len(picked)} markers cited twice"
                            if thin else ""))

        for sym, n_pmid in picked.items():
            marker_rows.append({
                "lineage_id": lid, "cell_of_origin": row.cell_of_origin,
                "symbol": sym, "ensembl_gene": sym2ens.get(sym),
                "n_papers": int(n_pmid),
                "cell_names": " | ".join(names),
                "tissue_filter": tissue or "",
                "used": status == "scored"})

        rec = {"lineage_id": lid, "cell_of_origin": row.cell_of_origin,
               "n_markers": int(len(picked)), "n_multi_paper": n_multi,
               "status": status, "thin_markers": thin, "reason": reason}
        status_rows.append(rec)
        seen[lid] = rec | {"cancer_types": [row.cancer_type]}

    st = pd.DataFrame(status_rows)
    st["cancer_types"] = st.lineage_id.map(
        lambda l: "; ".join(sorted(seen[l]["cancer_types"])))
    st = st.sort_values(["status", "thin_markers", "lineage_id"])
    st.to_csv(RESULTS / "step12_lineage_status.csv", index=False)

    mk = pd.DataFrame(marker_rows)
    unmapped = mk[mk.ensembl_gene.isna()]
    mk.to_csv(RESULTS / "step12_markers_used.csv", index=False)

    scored = st[st.status == "scored"]
    print(f"lineages defined                : {len(st)}")
    print(f"  scored                        : {len(scored)}"
          f"  (of which thin: {int(scored.thin_markers.sum())})")
    print(f"  not scored                    : {int((st.status=='not_scored').sum())}")
    nomulti = st[st.reason.str.contains("2 or more papers", na=False)]
    print(f"    for lack of multi-paper markers: {len(nomulti)}")
    print(f"markers selected                : {len(mk)} "
          f"({mk.symbol.nunique()} distinct genes)")
    if len(unmapped):
        print(f"  WITHOUT an Ensembl id         : {len(unmapped)} "
              f"({', '.join(unmapped.symbol)})")

    print("\nper lineage:")
    print(f"  {'lineage':<24}{'n':>3}{'multi':>6}  {'status':<11}{'markers'}")
    for _, r in st.iterrows():
        m = mk[(mk.lineage_id == r.lineage_id)]
        tag = " THIN" if r.thin_markers else ""
        print(f"  {r.lineage_id:<24}{r.n_markers:>3}{r.n_multi_paper:>6}  "
              f"{r.status + tag:<11}{', '.join(m.symbol.head(8))}")
        if r.reason:
            print(f"  {'':<24}{'':>9}  -> {r.reason}")

    print(f"\nwrote {(RESULTS/'step12_markers_used.csv').relative_to(ROOT)} "
          f"and {(RESULTS/'step12_lineage_status.csv').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
