"""Step 5 figures. Style from plot_style (D36).

The first figure is the one the runbook asks for in the report: every surface
candidate as a dot, tumour expression against worst-case healthy expression.
The ideal region is upper-left — high in tumour, low in every healthy tissue.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import (ACCENT, GREY, apply_style, pastel_viridis,
                        pastel_viridis_cmap, save_figure)

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

LABEL = ["EGFR", "ERBB2", "MSLN", "FOLR1", "CEACAM5", "TACSTD2", "NECTIN4",
         "CLDN6", "CA9", "GPC3", "CD70", "UPK3B"]


def fig_landscape(d: pd.DataFrame) -> Path:
    cf = d[d.carried_forward].copy()
    best = cf.loc[cf.groupby("symbol").S.idxmax()]

    fig, ax = plt.subplots(figsize=(7.6, 6.6))
    sc = ax.scatter(best.healthy_max_log2, best.median_tumour_log2tpm1,
                    c=best.S, cmap=pastel_viridis_cmap(0.18), s=16,
                    alpha=0.85, edgecolor="white", linewidth=0.25,
                    vmin=-4, vmax=4, rasterized=True)
    lim = [0, max(best.healthy_max_log2.max(),
                  best.median_tumour_log2tpm1.max()) * 1.04]
    ax.plot(lim, lim, color=GREY, lw=1.0, ls="--", zorder=1,
            label="equal in tumour and healthy tissue  (S = 0)")

    named = best[best.symbol.isin(LABEL)]
    ax.scatter(named.healthy_max_log2, named.median_tumour_log2tpm1,
               s=46, facecolor="none", edgecolor=ACCENT, linewidth=1.3, zorder=3)
    for _, r in named.iterrows():
        ax.annotate(r.symbol, (r.healthy_max_log2, r.median_tumour_log2tpm1),
                    textcoords="offset points", xytext=(6, 4), fontsize=8,
                    color="#2b2f36", zorder=4)

    ax.set(xlim=lim, ylim=lim,
           xlabel="worst healthy tissue, log2(TPM+1)   (less safe to the right)",
           ylabel="tumour median, log2(TPM+1)   (more target to hit, upward)",
           title="Step 5 — the target landscape\n"
                 f"{len(best):,} surface candidates, best cancer type each; "
                 "ideal region is upper-left")
    ax.legend(loc="upper left")
    cb = fig.colorbar(sc, ax=ax, shrink=0.72, pad=0.02)
    cb.set_label("specificity score S", fontsize=9)
    fig.tight_layout()
    return save_figure(fig, "step5_target_landscape")


def fig_variants(d: pd.DataFrame) -> Path:
    cf = d[d.carried_forward]
    best = cf.loc[cf.groupby("symbol").S.idxmax()]
    variants = [("S", "all tissues"), ("S_nobrain", "brain excluded"),
                ("S_no_reproductive_fetal", "testis + placenta excluded"),
                ("S_relaxed", "all three excluded")]
    cols = pastel_viridis(4, lightness=0.32)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.6, 4.3))
    for (c, lab), col in zip(variants, cols):
        ax1.hist(best[c], bins=60, range=(-8, 8), histtype="step",
                 linewidth=1.8, color=col, label=lab)
    ax1.axvline(0, color=GREY, lw=0.9, ls="--")
    ax1.set(xlabel="specificity score S", ylabel="surface candidates",
            title="Effect of excluding special tissues (D46)")
    ax1.legend(fontsize=8)

    known = best[best.symbol.isin(LABEL)].sort_values("S")
    y = np.arange(len(known))
    ax2.hlines(y, known.S, known.S_relaxed, color=GREY, lw=1.1, zorder=1)
    ax2.scatter(known.S, y, s=38, color=cols[0], zorder=2, label="all tissues")
    ax2.scatter(known.S_relaxed, y, s=38, color=cols[3], zorder=2,
                label="brain+testis+placenta excluded")
    ax2.axvline(0, color=ACCENT, lw=0.9, ls=":")
    ax2.set_yticks(y)
    ax2.set_yticklabels(known.symbol, fontsize=8.5)
    ax2.set(xlabel="specificity score S",
            title="Known targets move most when brain\nand placenta are excluded")
    ax2.legend(fontsize=8, loc="lower right")
    fig.tight_layout()
    return save_figure(fig, "step5_score_variants")


def main() -> int:
    apply_style()
    d = pd.read_csv(RESULTS / "step5_specificity.csv")
    for fn in (fig_landscape, fig_variants):
        print(f"  wrote {fn(d).relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
