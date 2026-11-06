"""EXTRA: export the candidate set for Mantis (mantis.csail.mit.edu).

Not a runbook step and nothing downstream reads it. Mantis embeds a text
column, projects to 2D and draws a map where proximity means similarity, so
the export has to supply something worth embedding.

That is the whole difficulty. Every result file in this project is numbers and
short labels; embedding `symbol` would cluster CLDN6 beside CLDN3 because the
strings look alike, which is a map of gene naming, not biology. So the
descriptive text is joined back in from the raw sources — UniProt protein
names and subcellular location, HPA's process/function/class annotations — and
concatenated into one `description` column.

Two files, because they answer different questions:

  mantis_candidates.csv   one row per eligible gene x cancer type (3,442).
                          The same gene appears once per cancer type it passes
                          in, so CLDN6 shows up as both a testicular and an
                          ovarian candidate and you can see them separately.

  mantis_genes.csv        one row per eligible gene (1,272), through its best
                          cancer type. One point per gene, which is the
                          cleaner map when the question is "where does my top
                          10 sit in the surfaceome".

Both carry coordinate1/coordinate2 (tumour and worst-healthy expression), so
the same file can also drive the explicit scatter the runbook asks for in
Step 11 instead of an embedding layout.

Every value is read from the committed results; nothing is typed.
"""

from __future__ import annotations

import io
import re
import zipfile
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

MAX_DESC = 400   # the upload endpoint rejects bodies over ~1 MB (nginx default)
ROUND = 3        # full float64 text is ~17 chars a value; 3 dp is plenty here

# UniProt location strings carry evidence codes like {ECO:0000250|...} that
# add nothing to an embedding and a lot of tokens.
ECO = re.compile(r"\{[^}]*\}")
NOTE = re.compile(r"\s*Note=.*$", re.S)


def clean_location(s: str) -> str:
    if not isinstance(s, str):
        return ""
    s = s.replace("SUBCELLULAR LOCATION:", "")
    s = ECO.sub("", s)
    s = NOTE.sub("", s)
    return re.sub(r"\s+", " ", s).strip(" .;")


def uniprot_names() -> pd.DataFrame:
    """Protein names live only in the raw download, not the processed table."""
    import gzip
    with gzip.open(RAW / "uniprot_human_reviewed_topology.tsv.gz", "rt",
                   encoding="utf8") as fh:
        u = pd.read_csv(fh, sep="\t", usecols=["Entry", "Protein names"])
    return u.rename(columns={"Entry": "accession",
                             "Protein names": "protein_name"})


def hpa_text() -> pd.DataFrame:
    with zipfile.ZipFile(RAW / "proteinatlas.tsv.zip") as z, \
            z.open(z.namelist()[0]) as fh:
        h = pd.read_csv(io.TextIOWrapper(fh, encoding="utf8"), sep="\t",
                        low_memory=False)
    keep = ["Ensembl", "Gene description", "Protein class",
            "Biological process", "Molecular function", "Disease involvement"]
    h = h[[c for c in keep if c in h.columns]].copy()
    return h.rename(columns={"Ensembl": "ensembl_gene"}).drop_duplicates("ensembl_gene")


def build_description(r) -> str:
    """One readable sentence-ish blob per row — this is what gets embedded.

    Ordered most-specific first: what the protein is, then what it does, then
    where it lives, then the clinical context. Empty fields are dropped rather
    than leaving 'nan' in the text.
    """
    # Four fields only, and the omissions are deliberate rather than purely
    # for size:
    #   HPA "Protein class" is near-identical across rows ("Predicted membrane
    #     proteins"), so it adds length without separating anything.
    #   "Raised in <cancer type>" and the healthy tissue are already categoric
    #     columns; repeating them here would make the embedding cluster by
    #     cancer type instead of by protein biology, which is the one thing the
    #     map can show that the existing figures cannot.
    txt = lambda k: r[k] if isinstance(r.get(k), str) and r[k].strip() else None
    bits = [
        txt("protein_name"),
        f"Family: {r['gene_group']}" if txt("gene_group") else None,
        f"Function: {r['Molecular function']}" if txt("Molecular function") else None,
        f"Process: {r['Biological process']}" if txt("Biological process") else None,
        f"Location: {r['location_clean']}" if txt("location_clean") else None,
    ]
    out = " ".join(str(b).strip().rstrip(".") + "." for b in bits if b)
    return out[:MAX_DESC].rsplit(" ", 1)[0] if len(out) > MAX_DESC else out


def main() -> int:
    elig = pd.read_csv(RESULTS / "step7_ranking_full.csv")
    uni = pd.read_csv(PROC / "step3_uniprot_topology.csv")
    cand = pd.read_csv(PROC / "surface_candidates.csv")
    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)

    # --- UniProt text, reached through every Ensembl id an accession maps to
    #     (the MHC-region problem from D81: one accession, many ids) ---
    c = cand[cand.is_candidate].drop(columns=["ensembl_gene"])
    link = (gmap[gmap.accession.isin(set(c.accession))][["accession", "ensembl_gene"]]
            .dropna().drop_duplicates())
    j = (link.merge(c, on="accession")
         .merge(uni[["accession", "subcellular_location"]], on="accession")
         .merge(uniprot_names(), on="accession", how="left"))
    j["_rank"] = j.tier.map({"S1": 0, "S2": 1, "S3": 2, "Flagged": 3})
    j = (j.sort_values(["ensembl_gene", "_rank"], kind="mergesort")
         .drop_duplicates("ensembl_gene"))
    j["location_clean"] = j.subcellular_location.map(clean_location)

    d = elig.merge(
        j[["ensembl_gene", "accession", "protein_name", "location_clean",
           "protein_existence", "extracellular_reported", "compartment"]],
        on="ensembl_gene", how="left")
    d = d.merge(hpa_text(), on="ensembl_gene", how="left")

    # --- stability from the robustness grid, where it exists ---
    try:
        st = pd.read_csv(RESULTS / "step9_stability.csv").set_index("symbol").stability
        d["threshold_stability"] = d.symbol.map(st).fillna(0.0).round(3)
    except FileNotFoundError:
        d["threshold_stability"] = pd.NA

    # --- where each row sits in the published output ---
    top = pd.read_csv(RESULTS / "step7_top10.csv")[["ensembl_gene", "cancer_type"]]
    res = pd.read_csv(RESULTS / "step7_reserve.csv")[["ensembl_gene", "cancer_type"]]
    key = lambda f: set(zip(f.ensembl_gene, f.cancer_type))
    tk, rk = key(top), key(res)
    d["selection"] = [
        "top 10" if k in tk else "reserve" if k in rk else "eligible"
        for k in zip(d.ensembl_gene, d.cancer_type)]

    # Mantis's Color / Size / Opacity controls only accept NUMERIC columns, so
    # the categoricals worth colouring by are also exported as numbers. Without
    # these you cannot find the finalists on the map at all.
    d["is_finalist"] = d.selection.map(
        {"top 10": 2, "reserve": 1, "eligible": 0}).astype(int)
    d["global_rank"] = pd.NA
    order = (pd.read_csv(RESULTS / "step7_top10.csv").symbol.tolist()
             + pd.read_csv(RESULTS / "step7_reserve.csv").symbol.tolist())
    rank = {sym: i + 1 for i, sym in enumerate(order)}
    d.loc[d.selection != "eligible", "global_rank"] = (
        d.loc[d.selection != "eligible", "symbol"].map(rank))
    d["tier_num"] = d.tier.map({"S1": 1, "S2": 2, "S3": 3})
    d["comparator_num"] = d.comparator_grade.map(
        {"direct": 2, "imperfect": 1, "lineage_only": 0})
    d["is_low_power"] = d.low_power.astype(bool).astype(int)
    d["protein_unconfirmed"] = (
        d.protein_existence != "Evidence at protein level").astype(int)

    d["description"] = d.apply(build_description, axis=1)
    d["healthy_max_tpm"] = (2 ** d.healthy_max_log2 - 1).round(2)
    d["hpa_link"] = ("https://www.proteinatlas.org/" + d.ensembl_gene
                     + "-" + d.symbol.astype(str))
    d["uniprot_link"] = "https://www.uniprot.org/uniprotkb/" + d.accession.astype(str)
    d["label"] = d.symbol.astype(str) + " — " + d.cancer_type.astype(str)

    cols = [
        # title / semantic / links
        "label", "symbol", "description", "hpa_link", "uniprot_link",
        # coordinates for the explicit scatter
        "median_tumour_log2tpm1", "healthy_max_log2",
        # numerics
        "S", "S_relaxed", "log2fc", "fdr", "tumour_tpm", "healthy_max_tpm",
        "frac_tumours_above_10tpm", "tumour_tpm_q25",
        "healthy_breadth_above_10tpm", "gtex_tau", "threshold_stability",
        # numeric encodings of the categoricals, for Color / Size / Opacity
        "is_finalist", "global_rank", "tier_num", "comparator_num",
        "is_low_power", "protein_unconfirmed",
        # categoricals
        "selection", "cancer_type", "tier", "comparator_grade", "compartment",
        "protein_existence", "agreement", "trajectory", "gene_group",
        "healthy_max_tissue", "healthy_max_source", "extracellular_reported",
        "low_power", "gpi", "cspa", "surfy",
        # keys
        "ensembl_gene", "accession",
    ]
    out = d[[c for c in cols if c in d.columns]].rename(columns={
        "median_tumour_log2tpm1": "coordinate1_tumour_log2tpm",
        "healthy_max_log2": "coordinate2_healthy_log2tpm"})

    num = out.select_dtypes("number").columns
    out[num] = out[num].round(ROUND)

    p1 = RESULTS / "mantis_candidates.csv"
    out.to_csv(p1, index=False)

    # one row per gene, through its best cancer type
    best = out.loc[out.groupby("ensembl_gene").S.idxmax()].sort_values(
        "S", ascending=False)
    p2 = RESULTS / "mantis_genes.csv"
    best.to_csv(p2, index=False)

    # The full gene x cancer-type file is ~2.3 MB and the upload endpoint
    # refuses bodies over roughly 1 MB, so ship the per-cancer-type shortlist
    # as well: the top 10 by S within each of the 30 cancer types. Small, and
    # it is the more interesting map anyway — the best candidate in every
    # indication rather than a long tail of genes that merely cleared the bar.
    per = (out.sort_values("S", ascending=False, kind="mergesort")
           .groupby("cancer_type", sort=True).head(10))
    p3 = RESULTS / "mantis_per_cancer.csv"
    per.to_csv(p3, index=False)

    leaked = out.description.str.contains(r"\bnan\b", case=False, na=False)
    assert not leaked.any(), \
        f"{int(leaked.sum())} descriptions contain a literal 'nan'"
    assert out.description.notna().all() and (out.description.str.len() > 20).all()

    n_txt = int((out.description.str.len() > 80).sum())
    mb = lambda p: p.stat().st_size / 1048576
    print(f"wrote {p1.relative_to(ROOT)}   {len(out):,} rows x "
          f"{out.shape[1]} cols   {mb(p1):.2f} MB")
    print(f"wrote {p2.relative_to(ROOT)}        {len(best):,} rows"
          f"{'':>14}{mb(p2):.2f} MB")
    print(f"wrote {p3.relative_to(ROOT)}  {len(per):,} rows"
          f"{'':>16}{mb(p3):.2f} MB")
    for p in (p2, p3):
        if mb(p) > 0.95:
            print(f"  WARNING: {p.name} is over the ~1 MB upload limit")
    print(f"\n  {p1.name} is {mb(p1):.1f} MB and will be refused on upload; "
          f"use {p2.name} or {p3.name}")
    print(f"\ndescription column: median {int(out.description.str.len().median())} chars, "
          f"{n_txt:,} rows over 80 chars ({n_txt/len(out):.0%})")
    print(f"selection breakdown: {out.selection.value_counts().to_dict()}")
    print("\nexample description:")
    print("  " + out.loc[out.selection == "top 10", "description"].iloc[0][:300])
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
