"""Step 4d: the three-state trajectory — healthy -> adjacent normal -> tumour.

Aran et al. (Nat Commun 2017) showed adjacent normal (NAT) is not healthy
tissue but an intermediate state. D34 already uses that to justify two separate
arms. This goes one step further and asks, per gene and per cancer type, WHICH
shape the three medians make.

Why it matters for a target search, which is the point of the whole project:

  tumour_specific_up   healthy ~= NAT  <  tumour
      The clean profile. The gene is raised in the tumour and NOT in the
      tissue right next to it, so a drug has somewhere safe to miss.

  progressive_up       healthy  <  NAT  <  tumour
  field_effect_up      healthy  <  NAT ~= tumour
      A field effect: the surrounding tissue is already elevated. Same tumour
      fold-change as the clean case if you only compare tumour to healthy, but
      a worse safety profile, because the drug will also meet the target in
      histologically normal tissue around the tumour.

That distinction is invisible to a two-group test. It only appears with three
states, and it is exactly the kind of judgment the runbook asks a human to
bring. The classification is DATA — it is reported, not used to filter.

All three medians already exist in results/step4_tumor_vs_normal.csv (the
tumour median is shared between arms), so nothing is recomputed.

Writes results/step4_trajectory.csv and two figures.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import (ACCENT, GREY, PALETTE, apply_style, pastel_viridis,
                        save_figure)

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"

# A step counts as a real change at 2-fold. Deliberately softer than the D39
# headline threshold of log2FC >= 2: this describes the SHAPE of a trajectory,
# it does not decide whether a gene is a hit.
DELTA = 1.0

ORDER = ["tumour_specific_up", "progressive_up", "field_effect_up",
         "nat_transient_up", "flat", "nat_transient_down",
         "field_effect_down", "progressive_down", "tumour_specific_down"]


def classify(d1: np.ndarray, d2: np.ndarray, total: np.ndarray) -> np.ndarray:
    """d1 = NAT - healthy (field step); d2 = tumour - NAT (tumour step)."""
    up1, dn1 = d1 >= DELTA, d1 <= -DELTA
    up2, dn2 = d2 >= DELTA, d2 <= -DELTA
    flat1, flat2 = ~up1 & ~dn1, ~up2 & ~dn2
    return np.select(
        [
            flat1 & up2,            # healthy ~= NAT < tumour
            up1 & up2,              # healthy < NAT < tumour
            up1 & flat2,            # healthy < NAT ~= tumour
            up1 & dn2,              # up in NAT, back down in tumour
            flat1 & dn2,            # healthy ~= NAT > tumour
            dn1 & dn2,              # healthy > NAT > tumour
            dn1 & flat2,            # healthy > NAT ~= tumour
            dn1 & up2,              # down in NAT, back up in tumour
        ],
        [
            "tumour_specific_up", "progressive_up", "field_effect_up",
            "nat_transient_up", "tumour_specific_down", "progressive_down",
            "field_effect_down", "nat_transient_down",
        ],
        default="flat",
    )


def build() -> pd.DataFrame:
    d = pd.read_csv(RESULTS / "step4_tumor_vs_normal.csv")
    g = d[d.arm == "gtex"][["ensembl_gene", "symbol", "tier", "cancer_type",
                            "median_tumour_log2tpm1", "median_normal_log2tpm1",
                            "log2fc", "fdr", "passes_default"]].rename(
        columns={"median_normal_log2tpm1": "healthy",
                 "median_tumour_log2tpm1": "tumour",
                 "log2fc": "log2fc_vs_healthy", "fdr": "fdr_vs_healthy",
                 "passes_default": "passes_vs_healthy"})
    n = d[d.arm == "tcga_nat"][["ensembl_gene", "cancer_type",
                                "median_normal_log2tpm1", "log2fc", "fdr",
                                "passes_default"]].rename(
        columns={"median_normal_log2tpm1": "nat",
                 "log2fc": "log2fc_vs_nat", "fdr": "fdr_vs_nat",
                 "passes_default": "passes_vs_nat"})
    m = g.merge(n, on=["ensembl_gene", "cancer_type"], how="inner")

    m["field_step"] = m.nat - m.healthy          # d1
    m["tumour_step"] = m.tumour - m.nat          # d2
    m["total_step"] = m.tumour - m.healthy
    m["trajectory"] = classify(m.field_step.values, m.tumour_step.values,
                               m.total_step.values)
    return m


def fig_trajectory_mix(m: pd.DataFrame) -> Path:
    counts = (m.trajectory.value_counts()
              .reindex(ORDER).fillna(0).astype(int))
    cols = pastel_viridis(len(ORDER), lightness=0.40)
    fig, ax = plt.subplots(figsize=(8.2, 4.0))
    bars = ax.barh(range(len(ORDER)), counts.values, color=cols,
                   edgecolor="white", linewidth=0.6)
    ax.set_yticks(range(len(ORDER)))
    ax.set_yticklabels([o.replace("_", " ") for o in ORDER], fontsize=8.5)
    ax.invert_yaxis()
    ax.set_xlabel("gene × cancer-type pairs")
    ax.set_title("Three-state trajectory: healthy to adjacent normal to tumour\n"
                 f"{len(m):,} pairs across 15 cancer types, step threshold "
                 f"{DELTA:g} log2 (2-fold)")
    total = counts.sum()
    for i, (b, v) in enumerate(zip(bars, counts.values)):
        if v:
            ax.text(b.get_width() + total * 0.006, b.get_y() + b.get_height() / 2,
                    f"{v:,}  ({v / total:.1%})", va="center", fontsize=8,
                    color="#2b2f36")
    ax.set_xlim(0, counts.max() * 1.22)
    fig.tight_layout()
    return save_figure(fig, "step4_trajectory_mix")


def fig_trajectory_profiles(m: pd.DataFrame) -> Path:
    """Mean profile of each up-trajectory, plus where the real hits sit."""
    show = ["tumour_specific_up", "progressive_up", "field_effect_up"]
    cols = pastel_viridis(3, lightness=0.30)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10.4, 4.2))

    x = [0, 1, 2]
    for c, lab in zip(cols, show):
        s = m[m.trajectory == lab]
        if s.empty:
            continue
        prof = [s.healthy.median(), s.nat.median(), s.tumour.median()]
        ax1.plot(x, prof, "-o", color=c, lw=2.2, ms=7, label=f"{lab.replace('_', ' ')}  (n={len(s):,})")
    ax1.set_xticks(x)
    ax1.set_xticklabels(["healthy\n(GTEx)", "adjacent\nnormal", "tumour"])
    ax1.set_ylabel("median expression, log2(TPM+1)")
    ax1.set_title("Median profile of each rising trajectory")
    ax1.legend(fontsize=8)

    # Among genes that pass the D39 defaults against healthy tissue, which
    # trajectory are they actually on? This is the safety-relevant split.
    hits = m[m.passes_vs_healthy]
    counts = hits.trajectory.value_counts()
    keep = [o for o in ORDER if counts.get(o, 0) > 0]
    vals = [counts[o] for o in keep]
    bar_cols = pastel_viridis(len(keep), lightness=0.40)
    ax2.bar(range(len(keep)), vals, color=bar_cols, edgecolor="white", linewidth=0.6)
    ax2.set_xticks(range(len(keep)))
    ax2.set_xticklabels([k.replace("_", "\n") for k in keep], fontsize=7.5)
    ax2.set_ylabel("gene × cancer-type pairs")
    ax2.set_title(f"Trajectory of the {len(hits):,} hits that pass\n"
                  "FDR < 0.01 and log2FC ≥ 2 against healthy tissue")
    for i, v in enumerate(vals):
        ax2.text(i, v + max(vals) * 0.015, f"{v:,}", ha="center", fontsize=8)
    fig.tight_layout()
    return save_figure(fig, "step4_trajectory_profiles")


def main() -> int:
    apply_style()
    m = build()
    m.to_csv(RESULTS / "step4_trajectory.csv", index=False)

    counts = m.trajectory.value_counts()
    hits = m[m.passes_vs_healthy]
    hit_counts = hits.trajectory.value_counts()

    lines = [
        "STEP 4d — three-state trajectory (healthy -> adjacent normal -> tumour)",
        "Generated by src/step4_trajectory.py — do not edit by hand.",
        f"Step threshold {DELTA:g} log2 (2-fold). Classification is reported, "
        "never used to filter.",
        "",
        f"cancer types with all three states       : "
        f"{m.cancer_type.nunique()}",
        f"gene x cancer-type pairs                 : {len(m):,}",
        "",
        "trajectory mix (all pairs):",
    ]
    for o in ORDER:
        v = int(counts.get(o, 0))
        lines.append(f"  {o:<22} : {v:>7,}  ({v / len(m):>5.1%})")

    lines += [
        "",
        f"of the {len(hits):,} pairs that pass FDR<0.01 and log2FC>=2 against",
        "healthy tissue, the trajectory they are actually on:",
    ]
    for o in ORDER:
        v = int(hit_counts.get(o, 0))
        if v:
            lines.append(f"  {o:<22} : {v:>7,}  ({v / len(hits):>5.1%})")

    clean = int(hit_counts.get("tumour_specific_up", 0))
    field = int(hit_counts.get("progressive_up", 0) +
                hit_counts.get("field_effect_up", 0))
    lines += [
        "",
        f"  tumour-specific (adjacent tissue NOT raised) : {clean:,}",
        f"  field effect    (adjacent tissue IS raised)  : {field:,}",
        "",
        "The second group would look identical to the first in a two-group",
        "tumour-vs-healthy test, but carries a worse safety profile: the target",
        "is already present in histologically normal tissue around the tumour.",
    ]
    report = "\n".join(lines)
    print(report)
    (RESULTS / "step4_trajectory.txt").write_text(report + "\n")

    for p in (fig_trajectory_mix(m), fig_trajectory_profiles(m)):
        print(f"\n  wrote {p.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
