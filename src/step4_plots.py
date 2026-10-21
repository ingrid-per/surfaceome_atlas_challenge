"""Step 4c: diagnostic figures. Style comes from plot_style (D36).

Three figures:

  1. p-value histograms, one panel per comparison. The runbook asks for these
     as the main diagnostic. Read them with the caveat that most genes here are
     NOT null: we pre-selected 3,912 surface proteins and compare whole
     different tissues, so a flat tail is not the expectation.
  2. power vs significance — the plot that shows the heavy p<0.01 mass is
     sample size, not a batch artefact, and that the effect-size threshold
     removes the dependence.
  3. arm concordance — GTEx log2FC against TCGA-adjacent-normal log2FC, which
     D34 makes a result rather than a nuisance.
  4. volcano per cancer type, with the D39 thresholds drawn.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import ACCENT, GREY, PALETTE, apply_style, pastel_viridis, save_figure

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

FDR_MAX, LOG2FC_MIN = 0.01, 2.0


def fig_pvalue_histograms(d: pd.DataFrame) -> Path:
    groups = sorted(d.groupby(["cancer_type", "arm"]).groups)
    ncol = 6
    nrow = int(np.ceil(len(groups) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 2.3, nrow * 1.9),
                             sharex=True, sharey=False)
    axes = np.atleast_1d(axes).ravel()
    colours = {"gtex": PALETTE[1], "tcga_nat": PALETTE[4]}
    for ax, key in zip(axes, groups):
        g = d[(d.cancer_type == key[0]) & (d.arm == key[1])]
        ax.hist(g.p_value, bins=20, range=(0, 1),
                color=colours.get(key[1], PALETTE[2]), edgecolor="white",
                linewidth=0.4)
        ax.axhline(len(g) / 20, color=GREY, lw=0.8, ls=":")
        short = key[0].replace(" Carcinoma", "").replace(" Adenocarcinoma", "")
        ax.set_title(f"{short[:22]}\n{key[1]}  n={g.n_normal.iloc[0]}",
                     fontsize=7, pad=4)
        ax.set_yticks([])
    for ax in axes[len(groups):]:
        ax.axis("off")
    fig.suptitle("Step 4 — p-value distributions per comparison\n"
                 "dotted line = uniform expectation if every gene were null",
                 fontsize=11, y=1.005)
    fig.supxlabel("Wilcoxon p-value", fontsize=9)
    fig.tight_layout()
    return save_figure(fig, "step4_pvalue_histograms")


def fig_power(d: pd.DataFrame) -> Path:
    rows = []
    for (c, a), g in d.groupby(["cancer_type", "arm"]):
        rows.append(dict(arm=a, n_eff=min(g.n_tumour.iloc[0], g.n_normal.iloc[0]),
                         frac_sig=(g.p_value < 0.01).mean(),
                         frac_pass=g.passes_default.mean()))
    r = pd.DataFrame(rows)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.5, 3.9))
    for arm, col in (("gtex", PALETTE[1]), ("tcga_nat", PALETTE[4])):
        s = r[r.arm == arm]
        ax1.scatter(s.n_eff, s.frac_sig * 100, s=42, color=col, alpha=0.85,
                    edgecolor="white", linewidth=0.6, label=arm)
        ax2.scatter(s.n_eff, s.frac_pass * 100, s=42, color=col, alpha=0.85,
                    edgecolor="white", linewidth=0.6, label=arm)
    rho1 = spearmanr(r.n_eff, r.frac_sig).statistic
    rho2 = spearmanr(r.n_eff, r.frac_pass).statistic
    ax1.set(xscale="log", xlabel="samples in the smaller group",
            ylabel="genes with p < 0.01  (%)",
            title=f"Significance alone tracks power\nSpearman ρ = {rho1:.2f}")
    ax2.set(xscale="log", xlabel="samples in the smaller group",
            ylabel="genes passing FDR < 0.01 and log2FC ≥ 2  (%)",
            title=f"Adding the effect-size threshold removes it\n"
                  f"Spearman ρ = {rho2:.2f}")
    ax1.legend(title="arm")
    fig.tight_layout()
    return save_figure(fig, "step4_power_vs_significance")


def fig_arm_concordance(d: pd.DataFrame) -> Path:
    g = d[d.arm.isin(["gtex", "tcga_nat"])]
    wide = g.pivot_table(index=["ensembl_gene", "cancer_type"], columns="arm",
                         values="log2fc")
    wide = wide.dropna()
    if wide.empty:
        return Path()
    fig, ax = plt.subplots(figsize=(5.4, 5.0))
    ax.scatter(wide.tcga_nat, wide.gtex, s=4, alpha=0.18,
               color=PALETTE[2], edgecolor="none", rasterized=True)
    lim = float(np.nanpercentile(np.abs(wide.values), 99.8))
    ax.plot([-lim, lim], [-lim, lim], color=GREY, lw=0.9, ls="--",
            label="equal in both arms")
    for v in (LOG2FC_MIN, -LOG2FC_MIN):
        ax.axhline(v, color=ACCENT, lw=0.7, ls=":")
        ax.axvline(v, color=ACCENT, lw=0.7, ls=":")
    rho = spearmanr(wide.tcga_nat, wide.gtex).statistic
    agree = float(((wide.gtex > 0) == (wide.tcga_nat > 0)).mean())
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim),
           xlabel="log2FC vs TCGA adjacent normal",
           ylabel="log2FC vs GTEx healthy tissue",
           title=f"Arm concordance (D34)\nSpearman ρ = {rho:.2f}   "
                 f"sign agreement = {agree:.1%}   n = {len(wide):,}")
    ax.legend(loc="upper left")
    fig.tight_layout()
    return save_figure(fig, "step4_arm_concordance")


def fig_volcano(d: pd.DataFrame) -> Path:
    g = d[d.arm == "gtex"].copy()
    top = (g.groupby("cancer_type").passes_default.sum()
           .sort_values(ascending=False).head(12).index.tolist())
    fig, axes = plt.subplots(3, 4, figsize=(13, 8.4))
    axes = axes.ravel()
    cols = pastel_viridis(3, lightness=0.35)
    for ax, c in zip(axes, top):
        s = g[g.cancer_type == c].copy()
        s["nlp"] = -np.log10(np.clip(s.fdr, 1e-300, None))
        ax.scatter(s.log2fc, s.nlp, s=5, color=GREY, alpha=0.25,
                   edgecolor="none", rasterized=True)
        hit = s[s.passes_default]
        ax.scatter(hit.log2fc, hit.nlp, s=11, color=cols[0], alpha=0.9,
                   edgecolor="white", linewidth=0.3)
        ax.axvline(LOG2FC_MIN, color=ACCENT, lw=0.7, ls=":")
        ax.axhline(-np.log10(FDR_MAX), color=ACCENT, lw=0.7, ls=":")
        ax.set_title(f"{c[:30]}  ({len(hit)} pass)", fontsize=8.5)
        ax.set_xlim(-8, 12)
    for ax in axes[len(top):]:
        ax.axis("off")
    fig.suptitle("Step 4 — tumour vs GTEx healthy tissue, 12 cancer types with "
                 "the most hits", fontsize=11, y=1.002)
    fig.supxlabel("log2 fold change (tumour − normal)", fontsize=9)
    fig.supylabel("−log10 FDR", fontsize=9)
    fig.tight_layout()
    return save_figure(fig, "step4_volcano")


def main() -> int:
    apply_style()
    d = pd.read_csv(RESULTS / "step4_tumor_vs_normal.csv")
    for fn in (fig_pvalue_histograms, fig_power, fig_arm_concordance, fig_volcano):
        p = fn(d)
        if p:
            print(f"  wrote {p.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
