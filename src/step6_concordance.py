"""Step 6b: do RNA and protein agree, for the same gene in the same cancer type?

D48: this is INFORMATION, not a filter. Nothing is dropped for disagreeing.
A gene where both layers agree is a stronger bet; a gene where only RNA moves
is a flag for whoever designs the next experiment, because RNA is a proxy and
the protein is what a drug actually binds.

D49: CPTAC's normal is tumour-ADJACENT tissue, so the matched RNA comparison is
the `tcga_nat` arm. Comparing protein-vs-adjacent-normal against
RNA-vs-GTEx-healthy would blame the molecular layer for a difference that is
really the comparator. The GTEx-healthy arm is carried alongside for reference.

Writes results/step6_rna_protein_concordance.csv and three figures.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import (ACCENT, GREY, apply_style, pastel_viridis,
                        pastel_viridis_cmap, save_figure)

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

RNA_FDR, RNA_LFC = 0.01, 2.0        # D39
PROT_FDR, PROT_LFC = 0.01, 1.0      # D50, scale-matched


def build() -> pd.DataFrame:
    p = pd.read_csv(RESULTS / "step6_protein_tumor_vs_normal.csv")
    r = pd.read_csv(RESULTS / "step4_tumor_vs_normal.csv")

    nat = r[r.arm == "tcga_nat"][
        ["ensembl_gene", "cancer_type", "log2fc", "fdr", "passes_default"]
    ].rename(columns={"log2fc": "rna_log2fc_nat", "fdr": "rna_fdr_nat",
                      "passes_default": "rna_passes_nat"})
    gte = r[r.arm == "gtex"][
        ["ensembl_gene", "cancer_type", "log2fc", "fdr", "passes_default", "S"]
        if "S" in r.columns else
        ["ensembl_gene", "cancer_type", "log2fc", "fdr", "passes_default"]
    ].rename(columns={"log2fc": "rna_log2fc_gtex", "fdr": "rna_fdr_gtex",
                      "passes_default": "rna_passes_gtex"})

    m = p.merge(nat, left_on=["ensembl_gene", "tcga_cancer_type"],
                right_on=["ensembl_gene", "cancer_type"], how="left")
    m = m.merge(gte, left_on=["ensembl_gene", "tcga_cancer_type"],
                right_on=["ensembl_gene", "cancer_type"], how="left",
                suffixes=("", "_g"))

    rna_up = (m.rna_fdr_nat < RNA_FDR) & (m.rna_log2fc_nat >= RNA_LFC)
    rna_dn = (m.rna_fdr_nat < RNA_FDR) & (m.rna_log2fc_nat <= -RNA_LFC)
    pro_up, pro_dn = m.protein_up.fillna(False), m.protein_down.fillna(False)
    no_rna = m.rna_log2fc_nat.isna()

    m["agreement"] = np.select(
        [no_rna,
         rna_up & pro_up, rna_dn & pro_dn,
         rna_up & pro_dn, rna_dn & pro_up,
         rna_up & ~pro_up & ~pro_dn, rna_dn & ~pro_up & ~pro_dn,
         ~rna_up & ~rna_dn & (pro_up | pro_dn)],
        ["rna_not_tested",
         "both_up", "both_down",
         "opposite", "opposite",
         "rna_only_up", "rna_only_down",
         "protein_only"],
        default="neither")
    return m


def fig_volcano(p: pd.DataFrame) -> Path:
    cohorts = sorted(p.cohort.unique())
    fig, axes = plt.subplots(2, 4, figsize=(14, 7.2))
    axes = axes.ravel()
    col = pastel_viridis(3, lightness=0.30)
    for ax, c in zip(axes, cohorts):
        s = p[p.cohort == c].copy()
        s["nlp"] = -np.log10(np.clip(s.protein_fdr, 1e-300, None))
        ax.scatter(s.protein_log2fc, s.nlp, s=4, color=GREY, alpha=0.22,
                   edgecolor="none", rasterized=True)
        sc = s[s.is_surface_candidate]
        ax.scatter(sc.protein_log2fc, sc.nlp, s=8, color=col[0], alpha=0.75,
                   edgecolor="none", rasterized=True, label="surface candidate")
        hit = sc[sc.protein_up]
        ax.scatter(hit.protein_log2fc, hit.nlp, s=16, color=ACCENT, alpha=0.95,
                   edgecolor="white", linewidth=0.3, label="up in tumour")
        ax.axvline(PROT_LFC, color=ACCENT, lw=0.7, ls=":")
        ax.axvline(-PROT_LFC, color=ACCENT, lw=0.7, ls=":")
        ax.axhline(-np.log10(PROT_FDR), color=ACCENT, lw=0.7, ls=":")
        ax.set_title(f"{c}  ({len(hit)} surface up)", fontsize=9)
        ax.set_xlim(-4.5, 4.5)
    axes[0].legend(fontsize=7, loc="upper left")
    for ax in axes[len(cohorts):]:
        ax.axis("off")
    fig.suptitle("Step 6 — protein, tumour vs adjacent normal (CPTAC)\n"
                 "all measured proteins in grey; surface candidates coloured",
                 fontsize=11, y=1.003)
    fig.supxlabel("protein log2 fold change", fontsize=9)
    fig.supylabel("-log10 FDR", fontsize=9)
    fig.tight_layout()
    return save_figure(fig, "step6_protein_volcano")


def fig_concordance(m: pd.DataFrame) -> Path:
    # Only surface candidates have an RNA result - Step 4 ran on the candidate
    # set - so every point here is a candidate. No "all proteins" layer exists
    # to compare against, and pretending otherwise would be misleading.
    s = m.dropna(subset=["rna_log2fc_nat", "protein_log2fc"])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.2, 5.0))
    ax1.scatter(s.rna_log2fc_nat, s.protein_log2fc, s=7, alpha=0.30,
                color=pastel_viridis(3, 0.28)[1], edgecolor="none",
                rasterized=True)
    rho = spearmanr(s.rna_log2fc_nat, s.protein_log2fc).statistic
    agree = float(((s.rna_log2fc_nat > 0) == (s.protein_log2fc > 0)).mean())
    ax1.axhline(0, color=GREY, lw=0.8, ls="--")
    ax1.axvline(0, color=GREY, lw=0.8, ls="--")
    ax1.set(xlim=(-7, 7), ylim=(-3.5, 3.5),
            xlabel="RNA log2FC  (tumour vs adjacent normal)",
            ylabel="protein log2FC  (tumour vs adjacent normal)",
            title=f"Same contrast, both layers — surface candidates\n"
                  f"Spearman ρ = {rho:.2f}   sign agreement = {agree:.0%}   "
                  f"n = {len(s):,} gene × cancer-type pairs")

    order = ["both_up", "rna_only_up", "protein_only", "neither",
             "opposite", "rna_only_down", "both_down"]
    sc_all = m[m.is_surface_candidate]
    counts = sc_all.agreement.value_counts()
    keep = [o for o in order if counts.get(o, 0) > 0]
    vals = [int(counts[o]) for o in keep]
    ax2.barh(range(len(keep)), vals, color=pastel_viridis(len(keep), 0.40),
             edgecolor="white", linewidth=0.6)
    ax2.set_yticks(range(len(keep)))
    ax2.set_yticklabels([k.replace("_", " ") for k in keep], fontsize=8.5)
    ax2.invert_yaxis()
    ax2.set_xlabel("surface candidate × cancer-type pairs")
    ax2.set_title("Where the surface candidates fall\n"
                  "(reported, never used to filter — D48)")
    for i, v in enumerate(vals):
        ax2.text(v + max(vals) * 0.01, i, f"{v:,}", va="center", fontsize=8)
    ax2.set_xlim(0, max(vals) * 1.14)
    fig.tight_layout()
    return save_figure(fig, "step6_rna_protein_concordance")


def main() -> int:
    apply_style()
    p = pd.read_csv(RESULTS / "step6_protein_tumor_vs_normal.csv")
    m = build()
    m.to_csv(RESULTS / "step6_rna_protein_concordance.csv", index=False)

    sc = m[m.is_surface_candidate]
    tested = sc[sc.agreement != "rna_not_tested"]
    both = sc.dropna(subset=["rna_log2fc_nat", "protein_log2fc"])
    rho = spearmanr(both.rna_log2fc_nat, both.protein_log2fc).statistic
    agree = float(((both.rna_log2fc_nat > 0) == (both.protein_log2fc > 0)).mean())

    lines = [
        "STEP 6b — RNA and protein concordance",
        "Generated by src/step6_concordance.py — do not edit by hand.",
        "Reported as information (D48). Nothing is dropped for disagreeing.",
        "Matched contrast: both layers are tumour vs ADJACENT NORMAL (D49).",
        "",
        f"rows where both layers were measured     : {len(both):,}",
        f"  distinct genes                         : {both.ensembl_gene.nunique():,}",
        f"  cancer types                           : {both.tcga_cancer_type.nunique()}",
        "",
        f"Spearman rho, RNA vs protein log2FC      : {rho:+.2f}",
        f"sign agreement                           : {agree:.1%}",
        "",
        "surface candidate x cancer-type pairs, by agreement:",
    ]
    for k, v in sc.agreement.value_counts().items():
        lines.append(f"  {k:<18} : {int(v):>6,}  ({v / len(sc):>5.1%})")
    lines += [
        "",
        "the practical question - of the RNA hits, how many does protein back?",
    ]
    rna_hits = sc[sc.agreement.isin(["both_up", "rna_only_up"])]
    if len(rna_hits):
        nb = int((rna_hits.agreement == "both_up").sum())
        lines += [
            f"  RNA says up (vs adjacent normal)       : {len(rna_hits):,}",
            f"    protein agrees                       : {nb:,}  "
            f"({nb / len(rna_hits):.0%})",
            f"    protein does not move                : {len(rna_hits) - nb:,}  "
            f"({1 - nb / len(rna_hits):.0%})",
            "",
            "  Neither group is discarded. The second is a flag for whoever",
            "  designs the next experiment: RNA is a proxy, protein is what a",
            "  drug binds, and mass spectrometry misses hydrophobic membrane",
            "  proteins, so 'protein does not move' is weak evidence either way.",
        ]
    report = "\n".join(lines)
    print(report)
    (RESULTS / "step6_concordance.txt").write_text(report + "\n")

    for fn in (fig_volcano, fig_concordance):
        arg = p if fn is fig_volcano else m
        print(f"\n  wrote {fn(arg).relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
