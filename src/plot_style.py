"""Shared figure style for the whole project (D36).

seaborn aesthetic, pastel viridis, clean. Defined once here and imported by
every plotting script, so the figures in the final report look like one set
rather than each script inventing its own look.

Usage:
    from plot_style import apply_style, save_figure, PALETTE, pastel_viridis
    apply_style()
    ...
    save_figure(fig, "step4_pvalue_histograms")
"""

from __future__ import annotations

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
from matplotlib.colors import LinearSegmentedColormap

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

DPI = 200


def pastel_viridis(n: int = 6, lightness: float = 0.45):
    """Viridis pulled toward white, so it reads as pastel but keeps viridis order.

    Plain viridis is heavily saturated and dark at one end. Blending each colour
    with white by `lightness` keeps the perceptually-uniform progression — which
    is the reason to use viridis at all — while softening it.
    """
    base = mpl.colormaps["viridis"](np.linspace(0.08, 0.92, n))[:, :3]
    return [tuple(c * (1 - lightness) + lightness) for c in base]


def pastel_viridis_cmap(lightness: float = 0.35) -> LinearSegmentedColormap:
    """Continuous version, for heatmaps and density."""
    base = mpl.colormaps["viridis"](np.linspace(0, 1, 256))[:, :3]
    soft = base * (1 - lightness) + lightness
    return LinearSegmentedColormap.from_list("pastel_viridis", soft)


# Six-step categorical palette; most figures need far fewer.
PALETTE = pastel_viridis(6)

# A single accent for thresholds and reference lines. Deliberately outside the
# viridis range so a threshold never reads as another data series.
ACCENT = "#c1666b"
GREY = "#8a8f98"


def apply_style() -> None:
    sns.set_theme(
        context="notebook",
        style="whitegrid",
        palette=PALETTE,
        font_scale=0.95,
    )
    mpl.rcParams.update({
        "figure.dpi": 110,
        "savefig.dpi": DPI,
        "savefig.bbox": "tight",
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "axes.edgecolor": GREY,
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.titlesize": 11,
        "axes.titleweight": "semibold",
        "axes.titlepad": 10,
        "axes.labelcolor": "#2b2f36",
        "axes.labelsize": 9.5,
        "grid.color": "#e6e8ec",
        "grid.linewidth": 0.7,
        "xtick.color": "#2b2f36",
        "ytick.color": "#2b2f36",
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "legend.frameon": False,
        "legend.fontsize": 8.5,
        "patch.edgecolor": "white",
        "patch.linewidth": 0.4,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica Neue", "Helvetica", "Arial", "DejaVu Sans"],
    })


def save_figure(fig, name: str, subdir: str = "figures") -> Path:
    """Write a figure to results/figures/<name>.png and return the path."""
    out_dir = RESULTS / subdir
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{name}.png"
    fig.savefig(path, dpi=DPI, bbox_inches="tight")
    plt.close(fig)
    return path
