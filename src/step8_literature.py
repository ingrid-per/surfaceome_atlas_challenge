"""Step 8a: an evidence check on every ranked candidate.

Two independent passes, neither of which relies on anything remembered:

1. LITERATURE. Query PubMed's E-utilities for each (gene, disease) pair and
   record how many papers exist. A target with 90 papers is well trodden; one
   with 0 is either novel or wrong, and the ranking cannot tell those apart.
   Counts come from NCBI, not from recall.

2. RED FLAGS. Structural problems the ranking cannot see, applied uniformly:
   a short exposed region, no plasma-membrane evidence, heterogeneous
   expression, a weak comparator, a worst-tissue that is also the tumour's own
   organ, and gene classes that are common artefacts.

The disease terms are the literature's vocabulary, not TCGA's labels. Querying
"Kidney Clear Cell Carcinoma" finds far less than "clear cell renal cell
carcinoma", and that mismatch was silently breaking the deliverable's links.

Writes results/step8_literature.csv and results/step8_redflags.csv.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
UA = "surface-prot-challenge/0.1 (academic target discovery)"
PAUSE = 0.36          # NCBI allows ~3 requests/second without an API key

# TCGA label -> the term the literature actually uses. Also used for the
# PubMed and ClinicalTrials.gov links in the deliverable.
DISEASE_TERM = {
    "Adrenocortical Cancer": "adrenocortical carcinoma",
    "Bladder Urothelial Carcinoma": "bladder urothelial carcinoma",
    "Brain Lower Grade Glioma": "low grade glioma",
    "Breast Invasive Carcinoma": "breast cancer",
    "Cervical & Endocervical Cancer": "cervical cancer",
    "Cholangiocarcinoma": "cholangiocarcinoma",
    "Colon Adenocarcinoma": "colorectal cancer",
    "Diffuse Large B-Cell Lymphoma": "diffuse large B-cell lymphoma",
    "Esophageal Carcinoma": "esophageal carcinoma",
    "Glioblastoma Multiforme": "glioblastoma",
    "Head & Neck Squamous Cell Carcinoma": "head and neck squamous cell carcinoma",
    "Kidney Chromophobe": "chromophobe renal cell carcinoma",
    "Kidney Clear Cell Carcinoma": "clear cell renal cell carcinoma",
    "Kidney Papillary Cell Carcinoma": "papillary renal cell carcinoma",
    "Liver Hepatocellular Carcinoma": "hepatocellular carcinoma",
    "Lung Adenocarcinoma": "lung adenocarcinoma",
    "Lung Squamous Cell Carcinoma": "lung squamous cell carcinoma",
    "Mesothelioma": "mesothelioma",
    "Ovarian Serous Cystadenocarcinoma": "ovarian serous carcinoma",
    "Pancreatic Adenocarcinoma": "pancreatic ductal adenocarcinoma",
    "Pheochromocytoma & Paraganglioma": "pheochromocytoma",
    "Prostate Adenocarcinoma": "prostate cancer",
    "Rectum Adenocarcinoma": "rectal cancer",
    "Sarcoma": "soft tissue sarcoma",
    "Skin Cutaneous Melanoma": "melanoma",
    "Stomach Adenocarcinoma": "gastric cancer",
    "Testicular Germ Cell Tumor": "testicular germ cell tumor",
    "Thymoma": "thymoma",
    "Thyroid Carcinoma": "thyroid carcinoma",
    "Uterine Carcinosarcoma": "uterine carcinosarcoma",
    "Uterine Corpus Endometrioid Carcinoma": "endometrial carcinoma",
    "Uveal Melanoma": "uveal melanoma",
}

THERAPY = "(immunotherapy OR antibody OR CAR-T OR target OR therapeutic)"


def pubmed_count(term: str) -> int | None:
    url = f"{EUTILS}?db=pubmed&retmode=json&term={urllib.parse.quote(term)}"
    req = urllib.request.Request(url, headers={"User-Agent": UA})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                d = json.loads(r.read())
            return int(d["esearchresult"]["count"])
        except Exception:
            time.sleep(1.2 * (attempt + 1))
    return None


def main() -> int:
    full = pd.read_csv(RESULTS / "step7_ranking_full.csv")
    per = (full.sort_values("S", ascending=False, kind="mergesort")
           .groupby("cancer_type", sort=True).head(10))
    top10 = pd.read_csv(RESULTS / "step7_top10.csv")
    rows = pd.concat([per, top10]).drop_duplicates(["ensembl_gene", "cancer_type"])
    print(f"checking {len(rows):,} candidate x cancer-type pairs "
          f"across {rows.cancer_type.nunique()} cancer types")

    out = []
    for i, (_, r) in enumerate(rows.iterrows(), 1):
        sym, cancer = str(r.symbol), r.cancer_type
        disease = DISEASE_TERM.get(cancer, cancer)
        n_gene = pubmed_count(f"{sym}[tiab]")
        time.sleep(PAUSE)
        n_pair = pubmed_count(f"{sym}[tiab] AND {disease}")
        time.sleep(PAUSE)
        n_ther = pubmed_count(f"{sym}[tiab] AND {disease} AND {THERAPY}")
        time.sleep(PAUSE)
        out.append({"symbol": sym, "ensembl_gene": r.ensembl_gene,
                    "cancer_type": cancer, "disease_term": disease,
                    "S": r.S, "tier": r.tier,
                    "pubmed_gene": n_gene, "pubmed_gene_disease": n_pair,
                    "pubmed_therapeutic": n_ther})
        if i % 25 == 0:
            print(f"  {i}/{len(rows)}", flush=True)

    lit = pd.DataFrame(out)
    lit.to_csv(RESULTS / "step8_literature.csv", index=False)

    # ---------------- red flags ----------------
    cand = pd.read_csv(PROC / "surface_candidates.csv")
    cand = cand[cand.is_candidate].drop_duplicates("ensembl_gene").set_index("ensembl_gene")
    flags = []
    for _, r in rows.iterrows():
        s = cand.loc[r.ensembl_gene] if r.ensembl_gene in cand.index else None
        f = []
        if s is not None:
            try:
                ecto = float(s.extracellular_reported)
                if ecto < 50:
                    f.append(f"short exposed region ({ecto:.0f} aa)")
            except (ValueError, TypeError):
                pass
            if s.compartment == "unannotated":
                f.append("no plasma-membrane evidence")
            if s.tm_sources == 1 and not s.gpi_any:
                f.append("only one source calls it transmembrane")
        if r.frac_tumours_above_10tpm < 0.7:
            f.append(f"heterogeneous ({r.frac_tumours_above_10tpm:.0%} of tumours)")
        if r.comparator_grade == "lineage_only":
            f.append("comparator is a different organ")
        if r.low_power:
            f.append("few normal samples")
        if str(r.symbol).startswith("OR") and str(r.symbol)[2:3].isdigit():
            f.append("olfactory receptor family - 7TM with small loops")
        flags.append({"symbol": r.symbol, "cancer_type": r.cancer_type,
                      "S": r.S, "n_flags": len(f), "flags": "; ".join(f)})
    fl = pd.DataFrame(flags).sort_values(["n_flags", "S"], ascending=[False, False])
    fl.to_csv(RESULTS / "step8_redflags.csv", index=False)
    print(f"\nwrote results/step8_literature.csv and step8_redflags.csv")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
