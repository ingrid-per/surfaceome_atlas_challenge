"""Step 4e: is the "neighbouring tissue is already high" pattern global, or a
few patients?

Everything so far has compared GROUP medians: all tumours against all normals.
That cannot tell whether every adjacent-normal sample is somewhat elevated, or
whether a handful of samples carry the whole effect and the rest look healthy.

94% of TCGA adjacent-normal samples have a tumour from the same patient in this
data, so the samples can be scored individually.

This is DESCRIPTIVE. It measures the spread and asks whether the group-level
signal survives when the most elevated samples are removed. It deliberately
draws no conclusion about WHY a sample is elevated — contamination,
inflammation, field cancerization and tissue-composition differences all remain
on the table and are not distinguished here.

Writes results/step4_nat_heterogeneity.csv, ..._samples.csv and two figures.
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import ACCENT, GREY, apply_style, pastel_viridis, save_figure

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

MIN_PAIRED = 10
MIN_GENES = 5
# A sample counts as "elevated" at 2-fold above the healthy median, the same
# step threshold the trajectory classification uses (D41).
ELEVATED = 1.0


def main() -> int:
    apply_style()
    mat = pd.read_parquet(PROC / "step4_candidate_tpm.parquet")
    samp = pd.read_parquet(PROC / "step4_samples.parquet")
    traj = pd.read_csv(RESULTS / "step4_trajectory.csv")
    tmap = pd.read_csv(ROOT / "docs" / "step4_tissue_map_DRAFT.tsv", sep="\t")

    pos = {c: i for i, c in enumerate(mat.columns)}
    V = mat.to_numpy(np.float32)
    gidx = {g: i for i, g in enumerate(mat.index)}

    tum = samp[(samp.study == "TCGA") & (samp.sample_type == "Primary Tumor")].copy()
    nat = samp[(samp.study == "TCGA")
               & (samp.sample_type == "Solid Tissue Normal")].copy()
    for d in (tum, nat):
        d["patient"] = d["sample"].str.rsplit("-", n=1).str[0]
    paired_patients = set(tum.patient) & set(nat.patient)

    gtex = samp[(samp.study == "GTEX") & (samp.sample_type == "Normal Tissue")]
    tissue_of = {r.tcga_cancer_type: [x.strip()
                                      for x in str(r.gtex_tissues).split("|")
                                      if x.strip()]
                 for _, r in tmap.iterrows()}

    rows, sample_rows = [], []
    for cancer, grp in nat.groupby("disease"):
        grp = grp[grp.patient.isin(paired_patients)]
        if len(grp) < MIN_PAIRED:
            continue
        # the genes that showed the pattern in THIS cancer type
        t = traj[(traj.cancer_type == cancer)
                 & traj.passes_vs_healthy
                 & traj.trajectory.isin(["field_effect_up", "progressive_up"])]
        genes = [g for g in t.ensembl_gene if g in gidx]
        if len(genes) < MIN_GENES:
            continue
        gi = [gidx[g] for g in genes]

        hids = [pos[s] for s in gtex.loc[gtex.disease.isin(tissue_of.get(cancer, [])),
                                         "sample"] if s in pos]
        nids = [pos[s] for s in grp["sample"] if s in pos]
        tids = [pos[s] for s in tum.loc[tum.disease == cancer, "sample"] if s in pos]
        if not hids or not nids:
            continue

        healthy_med = np.median(V[np.ix_(gi, hids)], axis=1)      # per gene
        tumour_med = np.median(V[np.ix_(gi, tids)], axis=1)

        # per-sample score: median across the pattern genes of
        # (this sample - healthy median for that gene)
        per_sample = np.median(V[np.ix_(gi, nids)] - healthy_med[:, None], axis=0)
        order = np.argsort(per_sample)

        for s, sc in zip(grp["sample"].tolist(), per_sample):
            sample_rows.append({"cancer_type": cancer, "sample": s,
                                "elevation_log2": float(sc),
                                "elevated": bool(sc >= ELEVATED)})

        # Does the group-level pattern survive without the most elevated quarter?
        keep = np.array(nids)[order[:max(1, int(len(nids) * 0.75))]]
        trimmed_med = np.median(V[np.ix_(gi, list(keep))], axis=1)
        full_med = np.median(V[np.ix_(gi, nids)], axis=1)

        rows.append({
            "cancer_type": cancer,
            "n_paired_normals": len(nids),
            "n_pattern_genes": len(genes),
            "elevation_median": float(np.median(per_sample)),
            "elevation_mean": float(np.mean(per_sample)),
            "elevation_iqr": float(np.percentile(per_sample, 75)
                                   - np.percentile(per_sample, 25)),
            "elevation_min": float(per_sample.min()),
            "elevation_max": float(per_sample.max()),
            "pct_samples_elevated": float((per_sample >= ELEVATED).mean()),
            # group-level gap, all samples vs bottom 75%
            "group_gap_full": float(np.median(full_med - healthy_med)),
            "group_gap_trimmed": float(np.median(trimmed_med - healthy_med)),
            "tumour_gap": float(np.median(tumour_med - healthy_med)),
        })

    summ = pd.DataFrame(rows).sort_values("n_paired_normals", ascending=False)
    smp = pd.DataFrame(sample_rows)
    summ["retained_after_trim"] = summ.group_gap_trimmed / summ.group_gap_full
    summ.to_csv(RESULTS / "step4_nat_heterogeneity.csv", index=False)
    smp.to_csv(RESULTS / "step4_nat_heterogeneity_samples.csv", index=False)

    # ---------------------------------------------------------------- figures
    order = summ.sort_values("elevation_median").cancer_type.tolist()
    cols = pastel_viridis(len(order), lightness=0.42)
    fig, ax = plt.subplots(figsize=(9.2, 0.42 * len(order) + 2.2))
    for i, (c, col) in enumerate(zip(order, cols)):
        v = smp.loc[smp.cancer_type == c, "elevation_log2"].values
        jitter = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.30
        ax.scatter(v, np.full(len(v), i) + jitter, s=16, color=col, alpha=0.75,
                   edgecolor="white", linewidth=0.3)
        ax.plot([np.median(v)], [i], marker="|", ms=20, color="#2b2f36", mew=2)
    ax.axvline(0, color=GREY, lw=0.9, ls="--")
    ax.axvline(ELEVATED, color=ACCENT, lw=0.9, ls=":")
    ax.set_yticks(range(len(order)))
    ax.set_yticklabels([c[:38] for c in order], fontsize=8.5)
    ax.set_xlabel("per-sample elevation above the healthy median, log2  "
                  "(one dot = one patient's adjacent-normal sample)")
    ax.set_title("Is the neighbouring-tissue signal global, or a few patients?\n"
                 "dashed = healthy level, dotted = 2-fold above; "
                 "bar = cancer-type median")
    fig.tight_layout()
    p1 = save_figure(fig, "step4_nat_heterogeneity")

    fig, ax = plt.subplots(figsize=(7.6, 4.4))
    s = summ.sort_values("pct_samples_elevated")
    y = np.arange(len(s))
    ax.barh(y, s.pct_samples_elevated * 100, color=pastel_viridis(len(s), 0.42),
            edgecolor="white", linewidth=0.6)
    ax.set_yticks(y)
    ax.set_yticklabels([c[:34] for c in s.cancer_type], fontsize=8)
    ax.set_xlabel("% of that cancer type's adjacent-normal samples "
                  "sitting ≥2-fold above healthy")
    ax.set_title("How widely shared is the elevation within each cancer type?")
    for i, (v, n) in enumerate(zip(s.pct_samples_elevated * 100, s.n_paired_normals)):
        ax.text(v + 1.2, i, f"{v:.0f}%  (n={n})", va="center", fontsize=7.5)
    ax.set_xlim(0, 112)
    fig.tight_layout()
    p2 = save_figure(fig, "step4_nat_elevation_prevalence")

    # ---------------------------------------------------------------- report
    lines = [
        "STEP 4e — how the neighbouring-tissue signal is distributed across patients",
        "Generated by src/step4_nat_heterogeneity.py — do not edit by hand.",
        "DESCRIPTIVE ONLY. This does not attempt to say WHY a sample is elevated.",
        "",
        f"cancer types examined                    : {len(summ)}",
        f"adjacent-normal samples scored           : {len(smp):,}",
        "",
        f"{'cancer type':<36} {'n':>4} {'genes':>6} {'median':>7} {'IQR':>6} "
        f"{'range':>14} {'%elev':>6} {'kept':>6}",
    ]
    for _, r in summ.iterrows():
        lines.append(
            f"  {r.cancer_type[:34]:<34} {r.n_paired_normals:>4} "
            f"{r.n_pattern_genes:>6} {r.elevation_median:>7.2f} "
            f"{r.elevation_iqr:>6.2f} "
            f"{r.elevation_min:>6.1f}..{r.elevation_max:<6.1f} "
            f"{r.pct_samples_elevated:>5.0%} {r.retained_after_trim:>5.0%}")
    lines += [
        "",
        "columns:",
        "  median / IQR  per-sample elevation above the healthy median, log2",
        "  %elev         share of samples at least 2-fold above healthy",
        "  kept          how much of the group-level gap remains after dropping",
        "                the most elevated 25% of samples. Near 100% means the",
        "                signal is broadly shared; a low value means a minority",
        "                of samples carried it.",
    ]
    report = "\n".join(lines)
    print(report)
    (RESULTS / "step4_nat_heterogeneity.txt").write_text(report + "\n")
    for p in (p1, p2):
        print(f"\n  wrote {p.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
