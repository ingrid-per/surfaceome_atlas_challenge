"""Step 12a: cell of origin and lineage home tissue, one row per cancer type.

Gate 1. This writes the table and nothing else — no expression is read, no
score is computed. The mapping is a set of judgements about tumour biology and
has to be approved before any number depends on it (runbook rule 3).

Two things are being recorded per cancer type, and they are different:

  cell_of_origin   the cell the tumour is thought to arise from
  home_tissue      the GTEx tissue in OUR data where that cell actually lives

They come apart. Uveal melanoma arises from melanocytes; there is no eye in
GTEx, and melanocytes live in skin. Cholangiocarcinoma arises from
cholangiocytes, which live inside liver but make up a small part of it. The
second case is the whole point of Step 12: when the lineage is a minority of
its home tissue, the healthy median is diluted and S is inflated for every gene
that lineage carries.

`lineage_share_of_home` records the EXPECTATION (majority or minority). The
analysis measures it independently through the marker composite, so the two can
be compared rather than the expectation being assumed.

Confidence:
  clear        one widely accepted origin, and a GTEx tissue that can represent it
  mixed        several or contested origins — not forced into one answer
  uncertain    origin accepted, but no adequate GTEx tissue exists

Only `clear` rows enter the enrichment test. `mixed` and `uncertain` rows are
scored where possible and reported separately, never silently pooled.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

MIN_HOME_SAMPLES = 10   # D37: the same floor Step 4 applies to any group

# cancer type, cell of origin, CellMarker search name, home tissue,
# expected share of home tissue, confidence, rationale
ROWS = [
    ("Adrenocortical Cancer", "adrenal cortical cell", "Adrenal cortical cell",
     "Adrenal Gland", "majority", "clear",
     "Cortex is the bulk of the adrenal gland, so little dilution expected."),
    ("Brain Lower Grade Glioma", "astrocyte / oligodendrocyte",
     "Astrocyte", "Brain - Cortex", "minority", "mixed",
     "LGG spans astrocytoma and oligodendroglioma; a single origin would be "
     "a forced choice."),
    ("Breast Invasive Carcinoma", "luminal epithelial cell",
     "Luminal epithelial cell", "Breast - Mammary Tissue", "minority", "clear",
     "Breast tissue is largely adipose and stroma; luminal epithelium is a "
     "small fraction."),
    ("Cervical & Endocervical Cancer", "squamous and glandular epithelium",
     "Squamous epithelial cell", "Esophagus - Mucosa", "majority", "mixed",
     "Squamous and adenocarcinoma are both common. GTEx cervix has 4-6 "
     "samples, so oesophageal mucosa stands in as stratified squamous "
     "epithelium, as in D40."),
    ("Cholangiocarcinoma", "cholangiocyte", "Cholangiocyte",
     "Liver", "minority", "clear",
     "Cholangiocytes line bile ducts inside liver but are a small percentage "
     "of liver cells. Strong dilution expected."),
    ("Colon Adenocarcinoma", "colonocyte", "Enterocyte",
     "Colon - Transverse", "majority", "clear",
     "Transverse colon is mucosa-bearing (D40); epithelium is well "
     "represented."),
    ("Diffuse Large B-Cell Lymphoma", "B cell", "B cell",
     "", "minority", "uncertain",
     "DLBCL arises in lymph node, which GTEx does not sample. Spleen was "
     "considered and rejected by me: it is a lymphoid organ but not the "
     "tissue this tumour comes from, so a shift measured against it would not "
     "mean what the column says it means."),
    ("Esophageal Carcinoma", "squamous and glandular epithelium",
     "Squamous epithelial cell", "Esophagus - Mucosa", "majority", "mixed",
     "TCGA ESCA pools squamous carcinoma and Barrett-associated "
     "adenocarcinoma, which have different origins."),
    ("Glioblastoma Multiforme", "astrocyte / OPC / neural stem cell",
     "Astrocyte", "Brain - Cortex", "minority", "mixed",
     "The cell of origin is actively contested between astrocyte, "
     "oligodendrocyte precursor and neural stem cell."),
    ("Head & Neck Squamous Cell Carcinoma", "squamous epithelial cell",
     "Squamous epithelial cell", "Esophagus - Mucosa", "majority", "clear",
     "Stratified squamous mucosa, matched to oesophageal mucosa as in the "
     "Step 4 map."),
    ("Kidney Chromophobe", "collecting-duct intercalated cell",
     "Intercalated cell", "Kidney - Cortex", "minority", "uncertain",
     "Origin in the distal nephron is accepted, but collecting ducts are "
     "largely medullary and GTEx samples cortex only."),
    ("Kidney Clear Cell Carcinoma", "proximal tubule epithelial cell",
     "Proximal tubule cell", "Kidney - Cortex", "majority", "clear",
     "Proximal tubule dominates renal cortex."),
    ("Kidney Papillary Cell Carcinoma", "proximal tubule epithelial cell",
     "Proximal tubule cell", "Kidney - Cortex", "majority", "clear",
     "Same lineage and tissue as clear cell."),
    ("Liver Hepatocellular Carcinoma", "hepatocyte", "Hepatocyte",
     "Liver", "majority", "clear",
     "Hepatocytes are the bulk of liver; minimal dilution expected."),
    ("Lung Adenocarcinoma", "alveolar type II cell", "Alveolar type II cell",
     "Lung", "minority", "clear",
     "AT2 cells are a minority of lung tissue by cell count."),
    ("Lung Squamous Cell Carcinoma", "bronchial basal cell", "Basal cell",
     "Lung", "minority", "clear",
     "Airway basal cells are a small fraction of whole lung."),
    ("Mesothelioma", "mesothelial cell", "Mesothelial cell",
     "", "minority", "uncertain",
     "Mesothelium is a single cell layer on serosal surfaces. No GTEx tissue "
     "represents it; visceral adipose contains some but is not mesothelium."),
    ("Ovarian Serous Cystadenocarcinoma", "fallopian tube secretory cell",
     "Secretory cell", "", "minority", "uncertain",
     "Current consensus places the origin in fallopian tube epithelium. GTEx "
     "Fallopian Tube has 5 samples, below the 10-sample floor (D37)."),
    ("Pancreatic Adenocarcinoma", "ductal epithelial cell", "Ductal cell",
     "Pancreas", "minority", "clear",
     "Pancreas is overwhelmingly acinar; ductal cells are a few percent. "
     "Strong dilution expected."),
    ("Pheochromocytoma & Paraganglioma", "chromaffin cell", "Chromaffin cell",
     "Adrenal Gland", "minority", "clear",
     "Chromaffin cells sit in the medulla, a minority of the adrenal gland. "
     "Strong dilution expected."),
    ("Prostate Adenocarcinoma", "luminal epithelial cell",
     "Luminal epithelial cell", "Prostate", "majority", "clear",
     "Glandular epithelium is well represented in prostate tissue."),
    ("Rectum Adenocarcinoma", "colonocyte", "Enterocyte",
     "Colon - Transverse", "majority", "clear",
     "Same lineage as colon; transverse colon used for the mucosal layer "
     "(D40)."),
    ("Sarcoma", "multiple mesenchymal lineages", "",
     "", "minority", "mixed",
     "TCGA SARC pools leiomyosarcoma, liposarcoma, UPS and others with "
     "different origins. No single lineage applies."),
    ("Skin Cutaneous Melanoma", "melanocyte", "Melanocyte",
     "Skin - Sun Exposed (Lower Leg)", "minority", "clear",
     "Melanocytes are a small fraction of skin, which is mostly "
     "keratinocytes."),
    ("Stomach Adenocarcinoma", "gastric epithelial cell",
     "Gastric epithelial cell", "Stomach", "majority", "clear",
     "Gastric mucosa is well represented in the GTEx stomach sample."),
    ("Testicular Germ Cell Tumor", "germ cell", "Germ cell",
     "Testis", "majority", "clear",
     "Germ cells at various stages dominate adult testis."),
    ("Thyroid Carcinoma", "thyroid follicular cell", "Follicular cell",
     "Thyroid", "majority", "clear",
     "Follicular cells are the bulk of thyroid."),
    ("Uterine Carcinosarcoma", "endometrial epithelium + mesenchymal",
     "Endometrial epithelial cell", "Uterus", "minority", "mixed",
     "Biphasic by definition, with both epithelial and mesenchymal "
     "components."),
    ("Uterine Corpus Endometrioid Carcinoma", "endometrial epithelial cell",
     "Endometrial epithelial cell", "Uterus", "minority", "clear",
     "GTEx uterus is whole myometrium-dominated wall; endometrium is a "
     "minority. Dilution expected."),
    ("Uveal Melanoma", "melanocyte", "Melanocyte",
     "Skin - Sun Exposed (Lower Leg)", "minority", "clear",
     "Uveal melanocytes have no GTEx tissue; cutaneous melanocytes in skin "
     "are the closest available lineage match, as in the Step 4 map."),
]

COLS = ["cancer_type", "cell_of_origin", "cellmarker_name", "home_tissue",
        "lineage_share_of_home", "confidence", "rationale"]


def main() -> int:
    RESULTS.mkdir(exist_ok=True)
    df = pd.DataFrame(ROWS, columns=COLS)

    # every cancer type carried by Step 7 must appear exactly once
    elig = pd.read_csv(RESULTS / "step7_ranking_full.csv")
    expected = set(elig.cancer_type.unique())
    got = set(df.cancer_type)
    assert not expected - got, f"cancer types missing from the table: {sorted(expected - got)}"
    assert not got - expected, f"cancer types not in Step 7: {sorted(got - expected)}"
    assert not df.cancer_type.duplicated().any(), "duplicated cancer type"

    # home tissue must be a real GTEx tissue in our data, with enough samples
    samp = pd.read_parquet(PROC / "step4_samples.parquet")
    g = samp[(samp.study == "GTEX") & (samp.sample_type == "Normal Tissue")]
    n_by_tissue = g.disease.value_counts().to_dict()

    bad = [t for t in df.home_tissue if t and t not in n_by_tissue]
    assert not bad, f"home tissue not present in GTEx data: {bad}"
    df["home_n_samples"] = df.home_tissue.map(lambda t: n_by_tissue.get(t, 0))

    # a tissue below the floor cannot represent a lineage, whatever the biology
    thin = df[(df.home_tissue != "") & (df.home_n_samples < MIN_HOME_SAMPLES)]
    assert thin.empty, (f"home tissue below the {MIN_HOME_SAMPLES}-sample floor: "
                        f"{thin[['cancer_type','home_tissue','home_n_samples']].to_dict('records')}")

    df["scorable"] = (df.home_tissue != "") & (df.cellmarker_name != "")
    df["in_enrichment_test"] = df.scorable & (df.confidence == "clear")
    df = df.sort_values(["confidence", "cancer_type"],
                        key=lambda s: s.map({"clear": 0, "mixed": 1, "uncertain": 2})
                        if s.name == "confidence" else s)
    df.to_csv(RESULTS / "step12_cell_of_origin.csv", index=False)

    rows_by_ct = elig.cancer_type.value_counts().to_dict()
    df["eligible_rows"] = df.cancer_type.map(rows_by_ct)

    print(f"cancer types                   : {len(df)}")
    print(f"  clear                        : {int((df.confidence=='clear').sum())}")
    print(f"  mixed                        : {int((df.confidence=='mixed').sum())}")
    print(f"  uncertain                    : {int((df.confidence=='uncertain').sum())}")
    print(f"scorable (home tissue + lineage): {int(df.scorable.sum())}")
    print(f"in the enrichment test          : {int(df.in_enrichment_test.sum())}")
    print(f"\neligible gene x cancer rows    : {len(elig):,}")
    print(f"  in scorable cancer types      : "
          f"{int(df.loc[df.scorable, 'eligible_rows'].sum()):,}")
    print(f"  in enrichment-test types      : "
          f"{int(df.loc[df.in_enrichment_test, 'eligible_rows'].sum()):,}")

    show = df[["cancer_type", "cell_of_origin", "home_tissue", "home_n_samples",
               "lineage_share_of_home", "confidence", "eligible_rows"]]
    print("\n" + show.to_string(index=False))
    print(f"\nwrote {(RESULTS / 'step12_cell_of_origin.csv').relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
