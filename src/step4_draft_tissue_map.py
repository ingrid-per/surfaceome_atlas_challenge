"""Step 4, prerequisite: draft the TCGA cancer type -> GTEx normal tissue map.

This is a DRAFT for review and correct by hand. Proposal P3 states the
mapping is a scientific judgment and must not be silently generated — so this
script supplies the counts and the literature precedent, and the verdict column
is hers.

Literature the draft is built on:

* Wang et al., "Unifying cancer and normal RNA sequencing data from different
  sources", Scientific Data 5:180061 (2018). Their published pairing file
  (mskcc/RNAseqDB, configuration/tissue-conf.txt) matches 14 GTEx tissues to 19
  TCGA studies. Rows traceable to it are marked source=Wang2018.

* Aran et al., "Comprehensive analysis of normal adjacent to tumor
  transcriptomes", Nature Communications 8:1077 (2017). Shows adjacent normal
  (NAT) is an intermediate state, not healthy tissue: on average 440 genes up
  and 554 down versus healthy, enriched for inflammatory signalling, and 55.3%
  MORE differentially expressed genes are found when healthy tissue is the
  comparator instead of NAT. Their conclusion, verbatim: "using healthy tissues
  is more accurate and provides additional information obscured when using NAT."

Rows not in Wang et al. are marked source=extended and carry my reasoning;
those are the ones most in need of review.

Counts come from the phenotype table, never typed by hand.
"""

from __future__ import annotations

import csv
import gzip
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
DOCS = ROOT / "docs"

# (tcga label) -> (gtex tissues, grade, source, rationale, second arm)
#
# grade records how good the comparator is, and travels with every result:
#   direct       same organ and compartment
#   imperfect    right organ, wrong compartment or a minority cell type
#   lineage_only different site, same cell lineage
#   none         no defensible comparator; skips Step 4, still scored in Step 5
#
# The grade is DATA, not a filter. Whether a weak comparison may support a
# finalist is a Step 7 question and is left open there (see P15).
MAP: dict[str, dict] = {
    # ---- direct: same organ and compartment -------------------------------
    "Breast Invasive Carcinoma": dict(
        gtex=["Breast - Mammary Tissue"], grade="direct", source="Wang2018",
        why="Direct tissue of origin."),
    "Lung Adenocarcinoma": dict(
        gtex=["Lung"], grade="direct", source="Wang2018", why="Direct."),
    "Lung Squamous Cell Carcinoma": dict(
        gtex=["Lung"], grade="direct", source="Wang2018",
        why="Direct; Wang pairs both lung subtypes to lung."),
    "Thyroid Carcinoma": dict(
        gtex=["Thyroid"], grade="direct", source="Wang2018", why="Direct."),
    "Prostate Adenocarcinoma": dict(
        gtex=["Prostate"], grade="direct", source="Wang2018", why="Direct."),
    "Stomach Adenocarcinoma": dict(
        gtex=["Stomach"], grade="direct", source="Wang2018", why="Direct."),
    "Liver Hepatocellular Carcinoma": dict(
        gtex=["Liver"], grade="direct", source="Wang2018", why="Direct."),
    "Colon Adenocarcinoma": dict(
        gtex=["Colon - Transverse", "Colon - Sigmoid"], grade="direct",
        source="Wang2018",
        why="Both colon subsites pooled; TCGA does not record subsite consistently."),
    "Rectum Adenocarcinoma": dict(
        gtex=["Colon - Sigmoid"], grade="direct", source="Wang2018",
        why="Wang pairs rectum to colon; sigmoid is the closest subsite."),
    "Uterine Corpus Endometrioid Carcinoma": dict(
        gtex=["Uterus"], grade="direct", source="Wang2018", why="Direct."),
    "Uterine Carcinosarcoma": dict(
        gtex=["Uterus"], grade="direct", source="Wang2018", why="Direct."),
    "Adrenocortical Cancer": dict(
        gtex=["Adrenal Gland"], grade="direct", source="extended",
        why="Adrenal cortex; direct."),
    "Testicular Germ Cell Tumor": dict(
        gtex=["Testis"], grade="direct", source="extended",
        why="Direct. Testis is immune-privileged with unusually broad expression, "
            "so little will look testis-specific."),
    "Kidney Clear Cell Carcinoma": dict(
        gtex=["Kidney - Cortex"], grade="direct", source="Wang2018",
        why="Direct, but only 28 GTEx samples - the thinnest well-matched normal."),
    "Kidney Papillary Cell Carcinoma": dict(
        gtex=["Kidney - Cortex"], grade="direct", source="Wang2018", why="As KIRC."),
    "Kidney Chromophobe": dict(
        gtex=["Kidney - Cortex"], grade="direct", source="Wang2018", why="As KIRC."),
    "Bladder Urothelial Carcinoma": dict(
        gtex=["Bladder"], grade="direct", source="Wang2018",
        why="Direct, but only 9 GTEx samples - too few for a stable median."),
    "Cervical & Endocervical Cancer": dict(
        gtex=["Cervix - Ectocervix", "Cervix - Endocervix"], grade="direct",
        source="Wang2018",
        why="Direct, but 6 + 4 GTEx samples. Expect this row to fail the minimum."),
    "Esophageal Carcinoma": dict(
        gtex=["Esophagus - Mucosa", "Esophagus - Gastroesophageal Junction"],
        grade="direct", source="Wang2018",
        why="Wang pairs ESCA to all three oesophageal subregions. Muscularis "
            "dropped: smooth muscle, not the epithelium these tumours arise from. "
            "ESCA mixes squamous (mucosa) and adenocarcinoma (GEJ) histologies."),
    "Brain Lower Grade Glioma": dict(
        gtex=["Brain - Cortex", "Brain - Frontal Cortex (Ba9)"], grade="direct",
        source="extended",
        why="Not in Wang2018. LGG arises in the cerebral hemispheres, so cortex is "
            "the right compartment; pooling all 13 GTEx brain regions would mix in "
            "cerebellum and basal ganglia, which differ sharply."),
    "Glioblastoma Multiforme": dict(
        gtex=["Brain - Cortex", "Brain - Frontal Cortex (Ba9)"], grade="direct",
        source="extended", why="As LGG."),
    "Head & Neck Squamous Cell Carcinoma": dict(
        gtex=["Esophagus - Mucosa", "Minor Salivary Gland"], grade="direct",
        source="extended",
        why="DEPARTS FROM Wang2018, which pairs HNSC to salivary gland alone. HNSC "
            "is mostly squamous mucosa; oesophageal mucosa is stratified squamous "
            "epithelium and the closest match. Salivary gland retained as well. "
            "Confirmed by me 2026-10-21."),

    # ---- imperfect: right organ, wrong compartment or minority cell -------
    "Cholangiocarcinoma": dict(
        gtex=["Liver"], grade="imperfect", source="Wang2018",
        why="Bile duct, not hepatocyte. Cholangiocytes are a minority of bulk "
            "liver, so liver-high genes may look falsely tumour-low."),
    "Pancreatic Adenocarcinoma": dict(
        gtex=["Pancreas"], grade="imperfect", source="Wang2018",
        why="GTEx pancreas is acinar-dominant while PDAC is ductal, and PDAC is "
            "heavily stromal. Expect a noisy comparison."),
    "Pheochromocytoma & Paraganglioma": dict(
        gtex=["Adrenal Gland"], grade="imperfect", source="extended",
        why="Arises from adrenal medulla and extra-adrenal paraganglia; GTEx "
            "adrenal gland is cortex-dominant. Right organ, wrong compartment."),
    "Skin Cutaneous Melanoma": dict(
        gtex=["Skin - Sun Exposed (Lower Leg)",
              "Skin - Not Sun Exposed (Suprapubic)"],
        grade="imperfect", source="extended",
        why="Melanocytes are a small minority of bulk skin, so melanocyte genes "
            "will look tumour-high on composition alone."),
    "Ovarian Serous Cystadenocarcinoma": dict(
        gtex=["Ovary"], grade="imperfect", source="extended",
        gtex_alt=["Fallopian Tube"],
        why="Ovary is conventional, but high-grade serous ovarian cancer is widely "
            "thought to arise from fallopian tube epithelium. Both run as parallel "
            "arms per I 2026-10-21; fallopian tube has only 5 samples."),

    # ---- lineage_only: different site, same cell lineage -----------------
    "Mesothelioma": dict(
        gtex=["Adipose - Visceral (Omentum)", "Lung"], grade="lineage_only",
        source="extended",
        why="GTEx has no pleura or peritoneum. The omentum carries a mesothelial "
            "lining, so visceral adipose genuinely contains the cell type MESO "
            "arises from; lung added as the usual anatomical neighbour."),
    "Uveal Melanoma": dict(
        gtex=["Skin - Sun Exposed (Lower Leg)",
              "Skin - Not Sun Exposed (Suprapubic)"],
        grade="lineage_only", source="extended",
        why="No ocular tissue in GTEx. UVM is melanocytic and skin is the only "
            "melanocyte-containing GTEx tissue: same lineage, wrong site."),
    "Sarcoma": dict(
        gtex=["Adipose - Subcutaneous", "Adipose - Visceral (Omentum)",
              "Muscle - Skeletal", "Nerve - Tibial"],
        grade="lineage_only", source="extended",
        why="SARC pools several mesenchymal lineages. These four cover the main "
            "ones, but the comparator is heterogeneous by construction."),
    "Diffuse Large B-Cell Lymphoma": dict(
        gtex=["Whole Blood", "Spleen"], grade="lineage_only", source="extended",
        why="Blood lineage. Lineage-matched comparator, but tagged so it is "
            "reported separately from solid tumours per the runbook."),

    # ---- none ------------------------------------------------------------
    "Thymoma": dict(
        gtex=[], grade="none", source="extended",
        why="GTEx has no thymus and thymic epithelium has no substitute - spleen "
            "would supply lymphocytes, the wrong cell. No differential test. Still "
            "scored in Step 5, whose specificity score needs no matched normal."),
}


def main() -> int:
    DOCS.mkdir(exist_ok=True)

    tumors: Counter = Counter()
    tcga_normal: Counter = Counter()
    gtex: Counter = Counter()
    with gzip.open(RAW / "TcgaTargetGTEX_phenotype.txt.gz", "rt",
                   encoding="latin-1") as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        for r in rd:
            study = (r.get("_study") or "").strip()
            stype = (r.get("_sample_type") or "").strip()
            disease = (r.get("primary disease or tissue") or "").strip()
            if study == "TCGA" and stype == "Primary Tumor":
                tumors[disease] += 1
            elif study == "TCGA" and stype == "Solid Tissue Normal":
                tcga_normal[disease] += 1
            elif study == "GTEX" and stype == "Normal Tissue":
                gtex[disease] += 1

    rows = []
    for cancer in sorted(tumors, key=lambda c: -tumors[c]):
        e = MAP.get(cancer, dict(gtex=[], grade="UNMAPPED", source="",
                                 why="not yet considered"))
        tissues, grade = e["gtex"], e["grade"]
        source, why = e["source"], e["why"]
        alt = e.get("gtex_alt", [])
        missing = [x for x in list(tissues) + list(alt) if x not in gtex]
        n_gtex = sum(gtex.get(t, 0) for t in tissues)
        rows.append({
            "tcga_cancer_type": cancer,
            "n_tumor": tumors[cancer],
            "n_tcga_adjacent_normal": tcga_normal.get(cancer, 0),
            "gtex_tissues": " | ".join(tissues),
            "n_gtex_normal": n_gtex,
            "comparator_grade": grade,
            "gtex_alt_tissues": " | ".join(alt),
            "n_gtex_alt": sum(gtex.get(x, 0) for x in alt),
            "source": source,
            "label_error": " | ".join(missing),
            "my_verdict": "",      # keep / change / drop  <- yours
            "rationale": why,
        })

    out = DOCS / "step4_tissue_map_DRAFT.tsv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # ---- summary to the terminal ----
    usable = [r for r in rows if r["gtex_tissues"]]
    none = [r for r in rows if not r["gtex_tissues"]]
    thin = [r for r in usable if r["n_gtex_normal"] < 30]
    extended = [r for r in usable if r["source"] == "extended"]
    from collections import Counter as _C
    grades = _C(r["comparator_grade"] for r in rows)
    print(f"TCGA cancer types                 : {len(rows)}")
    print(f"  with a proposed GTEx comparator : {len(usable)}")
    print(f"    from Wang et al. 2018         : "
          f"{len([r for r in usable if r['source'] == 'Wang2018'])}")
    print(f"    extended beyond the paper     : {len(extended)}")
    print(f"  no defensible comparator        : {len(none)}")
    print(f"  GTEx normals under 30 samples   : {len(thin)}")
    print(f"  comparator grades               : {dict(grades)}")
    bad = [r for r in rows if r["label_error"]]
    print(f"  label errors (tissue not found) : {len(bad)}")
    for r in bad:
        print(f"    !! {r['tcga_cancer_type']}: {r['label_error']}")
    print(f"\nWrote {out}")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
