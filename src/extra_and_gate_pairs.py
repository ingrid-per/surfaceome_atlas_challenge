"""EXTRA (not a runbook step): AND-gate pairs — two surface markers a conditionally-active drug could require.

Asked for by me, framed around what Orionis Biosciences' A-Kine platform
needs: cytokines engineered to act only in the presence of cancer-relevant
biomarkers. For that kind of drug a surface protein is an ANCHOR, not something
to destroy, and if activation requires two markers the unit of analysis is the
PAIR, not the gene.

A good pair, for each cancer type:

  1. both members eligible in that cancer type (Step 7 rules)
  2. CO-EXPRESSED across that cancer's tumours — a patient high in one is high
     in the other, so the pair is present together on the same tumours
  3. both present in most tumours of that type — an AND gate fires only where
     both are, so patient coverage is the product, not the max
  4. DISJOINT in healthy tissue — no single normal tissue carries both, which
     is the entire point: either marker alone may be unsafe, the combination
     need not be

Point 4 is what a single-gene ranking cannot express. The score here is the
worst healthy tissue for the PAIR: for each normal tissue take the LOWER of the
two genes (an AND gate is limited by its weaker input), then take the worst
tissue. A pair is safe where every tissue fails at least one of the two.

Reports only. No gate is applied to the main ranking.

Writes results/extra_and_gates.csv and results/extra_counts.txt.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

MIN_TUMOURS = 30          # a correlation needs enough tumours to mean anything
MAX_CANDIDATES = 60       # per cancer type, by S — keeps pairs tractable
MIN_COEXPR = 0.30         # Spearman across that cancer's tumours
MIN_BOTH_FRAC = 0.50      # both markers above 10 TPM in >= half the tumours
EXPRESSED_TPM = 10.0
TOP_PAIRS = 15            # reported per cancer type


def main() -> int:
    mat = pd.read_parquet(PROC / "step4_candidate_tpm.parquet")
    samp = pd.read_parquet(PROC / "step4_samples.parquet")
    full = pd.read_csv(RESULTS / "step7_ranking_full.csv")
    prof = pd.read_csv(RESULTS / "step5_healthy_profile.csv", index_col=0)

    tum = samp[(samp.study == "TCGA") & (samp.sample_type == "Primary Tumor")]
    gtex = samp[(samp.study == "GTEX") & (samp.sample_type == "Normal Tissue")]

    # healthy tissue profile, genes x tissues, in log2(TPM+1)
    cols = {}
    for tissue, grp in gtex.groupby("disease"):
        ids = [s for s in grp["sample"] if s in mat.columns]
        if ids:
            cols[tissue] = mat[ids].median(axis=1)
    healthy = pd.DataFrame(cols)

    rows = []
    for cancer, grp in full.groupby("cancer_type"):
        ids = [s for s in tum.loc[tum.disease == cancer, "sample"]
               if s in mat.columns]
        if len(ids) < MIN_TUMOURS:
            continue
        sel = (grp.sort_values("S", ascending=False)
               .drop_duplicates("ensembl_gene").head(MAX_CANDIDATES))
        genes = [g for g in sel.ensembl_gene if g in mat.index]
        if len(genes) < 2:
            continue

        E = mat.loc[genes, ids]
        tpm = np.power(2.0, E.to_numpy()) - 1.0
        # Spearman via ranks, done once for the whole block
        R = pd.DataFrame(E.to_numpy().T).rank().to_numpy()
        Rc = R - R.mean(axis=0)
        Rn = np.sqrt((Rc ** 2).sum(axis=0))
        Rn[Rn == 0] = np.nan
        C = (Rc.T @ Rc) / np.outer(Rn, Rn)

        present = tpm >= EXPRESSED_TPM          # genes x tumours, boolean
        sym = dict(zip(sel.ensembl_gene, sel.symbol))
        Sval = dict(zip(sel.ensembl_gene, sel.S))
        H = healthy.reindex(genes)

        n = len(genes)
        for i in range(n):
            for j in range(i + 1, n):
                r = C[i, j]
                if not np.isfinite(r) or r < MIN_COEXPR:
                    continue
                both = float((present[i] & present[j]).mean())
                if both < MIN_BOTH_FRAC:
                    continue
                # an AND gate is limited by its weaker input in each tissue
                pair_min = np.minimum(H.iloc[i].to_numpy(), H.iloc[j].to_numpy())
                k = int(np.nanargmax(pair_min))
                worst_pair = float(pair_min[k])
                worst_tissue = H.columns[k]
                # what each gene alone would have scored
                solo = float(max(np.nanmax(H.iloc[i].to_numpy()),
                                 np.nanmax(H.iloc[j].to_numpy())))
                med_t = float(np.median(np.minimum(tpm[i], tpm[j])))
                rows.append({
                    "cancer_type": cancer,
                    "gene_a": sym[genes[i]], "gene_b": sym[genes[j]],
                    "ensembl_a": genes[i], "ensembl_b": genes[j],
                    "coexpression_rho": round(float(r), 3),
                    "frac_tumours_both": round(both, 3),
                    "pair_tumour_tpm": round(med_t, 1),
                    "S_a": round(float(Sval[genes[i]]), 2),
                    "S_b": round(float(Sval[genes[j]]), 2),
                    "worst_healthy_pair_log2": round(worst_pair, 2),
                    "worst_healthy_pair_tissue": worst_tissue,
                    "worst_healthy_single_log2": round(solo, 2),
                    "safety_gain_log2": round(solo - worst_pair, 2),
                    "n_tumours": len(ids),
                })

    d = pd.DataFrame(rows)
    if d.empty:
        print("no pairs met the criteria")
        return 1
    # the pair's own specificity: tumour level minus the pair's worst tissue
    d["S_pair"] = np.log2(d.pair_tumour_tpm + 1) - d.worst_healthy_pair_log2
    d = d.sort_values(["cancer_type", "S_pair"], ascending=[True, False])
    d.to_csv(RESULTS / "step9_and_gates.csv", index=False)

    lines = [
        "STEP 9 — AND-gate pairs",
        "Generated by src/step9_and_gates.py — do not edit by hand.",
        "Two surface markers a conditionally-active drug could require together.",
        "Reported only; no gate is applied to the main ranking.",
        "",
        f"criteria: co-expression rho >= {MIN_COEXPR}; both markers above "
        f"{EXPRESSED_TPM:g} TPM in >= {MIN_BOTH_FRAC:.0%} of tumours;",
        f"          top {MAX_CANDIDATES} candidates per cancer type by S; "
        f">= {MIN_TUMOURS} tumours.",
        "",
        f"pairs found                              : {len(d):,}",
        f"cancer types with at least one pair      : {d.cancer_type.nunique()}",
        f"pairs where the pair is SAFER than either",
        f"  gene alone (safety gain > 1 log2)      : "
        f"{int((d.safety_gain_log2 > 1).sum()):,}",
        "",
        "The safety gain is the point: it is how much lower the worst healthy",
        "tissue becomes when a drug requires BOTH markers instead of either.",
        "",
        "best pair per cancer type, by the pair's own specificity S_pair:",
        "",
        f"  {'cancer type':<34} {'pair':<24} {'rho':>5} {'both%':>6} "
        f"{'S_pair':>7} {'gain':>6}  worst tissue for the pair",
    ]
    for cancer, g in d.groupby("cancer_type"):
        r = g.iloc[0]
        lines.append(
            f"  {cancer[:34]:<34} {r.gene_a + '+' + r.gene_b:<24} "
            f"{r.coexpression_rho:>5.2f} {r.frac_tumours_both:>6.0%} "
            f"{r.S_pair:>7.2f} {r.safety_gain_log2:>6.2f}  "
            f"{str(r.worst_healthy_pair_tissue)[:26]}")

    lines += ["", "largest safety gains overall "
              "(pair far safer than either gene alone):", ""]
    for _, r in d.nlargest(12, "safety_gain_log2").iterrows():
        lines.append(
            f"  {r.gene_a + ' + ' + r.gene_b:<26} {r.cancer_type[:30]:<30} "
            f"gain {r.safety_gain_log2:>5.2f} log2   "
            f"single worst {r.worst_healthy_single_log2:>5.2f} -> "
            f"pair {r.worst_healthy_pair_log2:>5.2f}")

    report = "\n".join(lines)
    print(report)
    (RESULTS / "step9_counts.txt").write_text(report + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
