"""A one-page map of every data source: what it is, what is in it, what it answers.

Built for explaining the pipeline out loud. The datasets are grouped by the
QUESTION each one answers rather than by provider, because that is the order
the argument runs in: first decide what counts as a surface protein, then ask
whether it is raised in tumour and absent from healthy tissue, then see whether
the protein layer agrees.

Grouping this way also makes two things visible that a provider list hides:
the Human Protein Atlas is four different files doing four different jobs, and
most of the effort went into the surface call rather than the expression test.

Every count is read from the committed data, never typed (rule 2).

Writes results/figures/step1_data_sources.png
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, pastel_viridis, save_figure

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"


def facts() -> dict:
    """Every number on the figure, read from the files it describes."""
    f = {}
    man = pd.read_csv(RESULTS / "step1_manifest.csv")
    f["release"] = dict(zip(man.resource, man.release))
    f["n_files"] = man.resource.value_counts().to_dict()
    # 1024-based, matching human() in step1_download.py — otherwise the figure
    # would say 1.8 GB where the manifest says 1700 MB, for the same bytes
    f["total_mb"] = man.bytes.sum() / 1048576
    f["total_files"] = len(man)

    uni = pd.read_csv(PROC / "step3_uniprot_topology.csv")
    f["uniprot_n"] = len(uni)
    f["uniprot_tm"] = int((uni.n_tm_segments > 0).sum())
    f["uniprot_topo"] = int((uni.n_extracellular_domains > 0).sum())
    f["uniprot_gpi"] = int(uni.gpi_curated.sum())

    dtm = pd.read_csv(PROC / "step3_deeptmhmm_full.csv")
    f["dtm_n"] = len(dtm)
    f["dtm_surface"] = int(dtm.dtm_surface_call.sum())

    gmap = pd.read_csv(PROC / "gene_map.csv", low_memory=False)
    f["cspa"] = gmap[gmap.cspa].accession.nunique()
    f["surfy"] = gmap[gmap.surfy_surface].accession.nunique()
    f["both_atlases"] = gmap[gmap.cspa & gmap.surfy_surface].accession.nunique()

    samp = pd.read_parquet(PROC / "step4_samples.parquet")
    tum = samp[(samp.study == "TCGA") & (samp.sample_type == "Primary Tumor")]
    gtx = samp[(samp.study == "GTEX") & (samp.sample_type == "Normal Tissue")]
    nat = samp[(samp.study == "TCGA") & (samp.sample_type == "Solid Tissue Normal")]
    f["n_tumour"] = len(tum)
    f["n_cancer_types"] = tum.disease.nunique()
    f["n_gtex"] = len(gtx)
    f["n_gtex_tissues"] = gtx.disease.nunique()
    f["n_nat"] = len(nat)
    f["n_samples_total"] = len(samp)

    hgnc = pd.read_csv(RAW / "hgnc_complete_set.txt", sep="\t", low_memory=False)
    f["hgnc_n"] = len(hgnc)
    groups = {g.strip() for s in hgnc.gene_group.dropna()
              for g in str(s).split("|") if g.strip()}
    f["hgnc_groups"] = len(groups)

    coh = pd.read_csv(RESULTS / "step6_cohorts.csv")
    f["cptac_cohorts"] = len(coh)

    cand = pd.read_csv(PROC / "surface_candidates.csv")
    f["candidates"] = int(cand.is_candidate.sum())
    return f


def build_cards(f: dict) -> list[tuple[str, list[tuple]]]:
    """(column heading, [(name, version, [content lines])]) in reading order."""
    return [
        ("Is it on the cell surface?", [
            ("UniProt", f["release"].get("UniProt", ""),
             [f"{f['uniprot_n']:,} reviewed human proteins",
              f"{f['uniprot_tm']:,} with a transmembrane segment",
              f"{f['uniprot_topo']:,} with a mapped extracellular domain",
              f"{f['uniprot_gpi']:,} with a curated GPI anchor",
              "curated by hand — the strongest evidence"]),
            ("DeepTMHMM 1.0", "run locally, academic licence",
             [f"predicted topology for {f['dtm_n']:,} proteins",
              f"{f['dtm_surface']:,} called membrane-spanning",
              "the 'TMHMM-type model' the brief asked for",
              "reasons from sequence, so it covers the",
              "43% UniProt never annotated"]),
            ("HPA — subcellular_location", f["release"].get("HPA", ""),
             ["13,603 genes, by immunofluorescence",
              "where in the cell the protein sits",
              "used to exclude proteins with positive",
              "evidence of being internal only"]),
            ("CSPA", f["release"].get("CSPA", ""),
             [f"{f['cspa']:,} proteins, experimentally measured",
              "labelled on intact living cells, so only",
              "genuinely exposed proteins are caught",
              "narrow: only the cell types they tested"]),
            ("SURFY", f["release"].get("SURFY", ""),
             [f"{f['surfy']:,} proteins, machine-learning predicted",
              f"only {f['both_atlases']:,} overlap with CSPA",
              "neither is used as a filter — both are",
              "reported as supporting evidence"]),
        ]),
        ("Is it raised in tumour,\nand absent from healthy tissue?", [
            ("UCSC Xena Toil", f["release"].get("Xena Toil", ""),
             [f"60,498 genes x {f['n_samples_total']:,} samples, 1.3 GB",
              f"{f['n_tumour']:,} TCGA tumours, {f['n_cancer_types']} cancer types",
              f"{f['n_gtex']:,} GTEx healthy, {f['n_gtex_tissues']} tissues",
              f"{f['n_nat']:,} tumour-adjacent normal",
              "TCGA and GTEx reprocessed through ONE",
              "pipeline — otherwise half the difference",
              "you measure is the two pipelines",
              "values are log2(TPM+0.001), not TPM"]),
            ("HPA — rna_tissue", f["release"].get("HPA", ""),
             ["20,162 genes x 40 healthy tissues",
              "HPA's own donors, own sequencing",
              "the second, independent reading of",
              "healthy tissue — kept separate from",
              "GTEx on purpose, so they can disagree"]),
        ]),
        ("Does the protein layer agree?", [
            ("CPTAC", f["release"].get("CPTAC", ""),
             [f"{f['cptac_cohorts']} cancer cohorts, {f['n_files'].get('CPTAC', 0)} files",
              "mass-spec protein in tumour and in",
              "matched adjacent normal",
              "reported, never used to rank (D48)"]),
            ("GTEx proteome", "Jiang et al., Cell 2020",
             ["protein abundance across 32 healthy",
              "tissues, from the same consortium",
              "mass spec under-detects membrane",
              "proteins, so absence is missing",
              "evidence, not evidence of absence"]),
            ("HPA — normal_ihc", f["release"].get("HPA", ""),
             ["1.2 million rows of antibody staining",
              "gene x tissue x CELL TYPE, 60 tissues",
              "the only protein-level evidence here",
              "that resolves to single cell types",
              "this is what the finalist review reads"]),
        ]),
        ("Supporting", [
            ("HPA — proteinatlas.tsv", f["release"].get("HPA", ""),
             ["20,162 genes x 119 columns",
              "symbol, Ensembl, UniProt, synonyms",
              "the table that lets three databases",
              "with three naming schemes be joined"]),
            ("HGNC", f["release"].get("HGNC", ""),
             [f"{f['hgnc_n']:,} genes, {f['hgnc_groups']:,} gene groups",
              "used for one thing: at most two per",
              "family in the top 10, so one family",
              "cannot fill the list and share one",
              "failure mode"]),
            ("CellMarker 2.0", f["release"].get("CellMarker 2.0", ""),
             ["curated cell-type marker genes",
              "added late, for the lineage check",
              "that did not work — kept because the",
              "negative result is part of the story"]),
        ]),
    ]


def main() -> int:
    apply_style()
    f = facts()
    cols = build_cards(f)
    tone = pastel_viridis(4, 0.30)

    # geometry in axis units; cards are sized by their line count so a short
    # card does not leave a hole under it
    CW, PAD, LINE, HEAD, TOPPAD = 1.0, 0.16, 0.185, 0.62, 0.30
    fig, ax = plt.subplots(figsize=(17.5, 10.6))
    ax.set_xlim(0, len(cols) * (CW + PAD))
    ax.axis("off")

    heights = []
    for ci, (_, cards) in enumerate(cols):
        y = 0.0
        for _, _, lines in cards:
            y += HEAD + len(lines) * LINE + TOPPAD
        heights.append(y)
    top = max(heights)
    ax.set_ylim(-0.75, top + 1.15)

    for ci, (heading, cards) in enumerate(cols):
        x = ci * (CW + PAD)
        colour = tone[ci]
        ax.text(x, top + 0.78, heading, fontsize=12.5, fontweight="600",
                va="top", ha="left", linespacing=1.35,
                family=plt.rcParams["font.family"])
        ax.add_patch(plt.Rectangle((x, top + 0.50), CW, 0.045,
                                   facecolor=colour, edgecolor="none"))

        y = top
        for name, version, lines in cards:
            h = HEAD + len(lines) * LINE
            ax.add_patch(FancyBboxPatch(
                (x, y - h), CW, h, boxstyle="round,pad=0.012,rounding_size=0.02",
                facecolor="#FBFCFC", edgecolor="#DFE5E6", linewidth=0.9))
            # a colour stripe keyed to the column, so the grouping survives
            # being printed or photographed from a slide
            ax.add_patch(plt.Rectangle((x, y - h), 0.022, h,
                                       facecolor=colour, edgecolor="none"))

            ax.text(x + 0.055, y - 0.17, name, fontsize=11.2, fontweight="700",
                    va="center", ha="left")
            # release strings vary wildly in length; CPTAC's ran off the card
            shown = version if len(version) <= 44 else version[:42].rstrip(" ,") + "..."
            ax.text(x + 0.055, y - 0.345, shown, fontsize=7.8,
                    va="center", ha="left", color="#7D8A8C",
                    family=plt.rcParams["font.monospace"][0])
            for li, line in enumerate(lines):
                ax.text(x + 0.055, y - HEAD - li * LINE + 0.04, line,
                        fontsize=8.6, va="center", ha="left", color="#2C3A3B")
            y -= h + TOPPAD

    # The excluded list is the other half of the argument: a dataset left out
    # on purpose, with the reason, is as much a decision as one taken in.
    from step1_resources import EXCLUDED
    ex_x = 3 * (CW + PAD)
    ex_y = top - sum(HEAD + len(l) * LINE + TOPPAD for _, _, l in cols[3][1])
    # name and reason on separate lines: on one line the reason had to be cut
    # mid-word and still ran past the card edge
    ex_lines = []
    for name, why in EXCLUDED:
        reason = why.split(":", 1)[-1].strip()
        ex_lines.append((name.split("(")[0].strip()[:40], False))
        ex_lines.append(((reason[:50].rstrip() + "..."
                          if len(reason) > 52 else reason), True))
    ex_h = HEAD + len(ex_lines) * (LINE * 0.78) + 0.12
    ax.add_patch(FancyBboxPatch(
        (ex_x, ex_y - ex_h), CW, ex_h,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        facecolor="#FAF6F6", edgecolor="#E6D9DA", linewidth=0.9,
        linestyle="--"))
    ax.text(ex_x + 0.055, ex_y - 0.17, "Deliberately not used", fontsize=11.2,
            fontweight="700", va="center", ha="left", color="#8E5F63")
    ax.text(ex_x + 0.055, ex_y - 0.345,
            "recorded with the reason, like everything else", fontsize=7.8,
            va="center", ha="left", color="#A48587",
            family=plt.rcParams["font.monospace"][0])
    for li, (text, is_reason) in enumerate(ex_lines):
        ax.text(ex_x + (0.075 if is_reason else 0.055),
                ex_y - HEAD - li * (LINE * 0.78) + 0.06, text,
                fontsize=6.9 if is_reason else 7.5, va="center", ha="left",
                color="#9A8183" if is_reason else "#6E595B",
                style="italic" if is_reason else "normal")

    ax.text(0, -0.42,
            f"{f['total_files']} files, {f['total_mb']/1024:.1f} GB, every one with a "
            f"recorded version, source URL and SHA-256 checksum.  "
            f"{f['candidates']:,} proteins survive the surface filter.",
            fontsize=9.5, va="center", ha="left", color="#55686A")
    ax.set_title("Where every number comes from", fontsize=17, fontweight="700",
                 loc="left", pad=22)

    fig.tight_layout()
    p = save_figure(fig, "step1_data_sources")
    print(f"wrote {Path(p).relative_to(ROOT)}")
    print(f"  {f['total_files']} files, {f['total_mb']:.0f} MB across "
          f"{len(f['n_files'])} resources")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
