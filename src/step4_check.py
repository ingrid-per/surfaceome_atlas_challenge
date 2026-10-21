"""Step 4 checks (register rows 4.1-4.7).

A note on 4.3, because it is the check most likely to be misread.

The runbook says a healthy p-value histogram shows "a spike near zero and
roughly flat elsewhere", and that another shape means a bug or a batch problem.
That heuristic assumes most genes are null. Here they are not: we pre-selected
3,912 surface proteins and compare entire different tissues, and several
comparisons have hundreds of samples per group. Between 50% and 89% of genes
fall below p < 0.01, and the tail near 1 is thin.

That is power, not an artefact, and 4.3 tests it as such:
  * the significant fraction should track sample size (it does, rho ~0.9)
  * that dependence should largely vanish once the effect-size threshold is
    applied (it does, rho ~0.2) — which is D39 doing its job
  * there should be no global directional shift (median log2FC ~ 0)
  * both directions should be represented among significant genes
A real batch effect would push everything one way and would NOT scale with n.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

MIN_N = 10
XENA_FLOOR = -9.9658
POSITIVE = ["EGFR", "ERBB2", "TACSTD2", "NECTIN4", "MSLN", "CLDN6", "CD19",
            "FOLR1", "CEACAM5"]


def main() -> int:
    d = pd.read_csv(RESULTS / "step4_tumor_vs_normal.csv")
    log = pd.read_csv(RESULTS / "step4_comparisons.csv")
    failures, checks = [], 0

    # 4.1 extraction integrity
    checks += 3
    sys.path.insert(0, str(ROOT / "src"))
    from step4_extract_matrix import to_log2_tpm1
    floor = float(to_log2_tpm1(np.array([XENA_FLOOR]))[0])
    mat = pd.read_parquet(PROC / "step4_candidate_tpm.parquet")
    vmin = float(mat.to_numpy().min())
    print(f"4.1  back-transform: floor {XENA_FLOOR} -> log2(TPM+1) = {floor:.2e}")
    print(f"     matrix {mat.shape[0]:,} genes x {mat.shape[1]:,} samples, "
          f"min value {vmin:.3f}")
    if abs(floor) > 1e-6:
        failures.append(f"the matrix floor should map to 0, got {floor:.3e} (D11)")
    if vmin < 0:
        failures.append(f"negative expression after back-transform: {vmin}")
    if mat.shape[1] != 19_131:
        failures.append(f"expected 19,131 sample columns, got {mat.shape[1]:,}")

    # 4.2 no comparison ran under the minimum
    checks += 1
    ran = log[log.status != "SKIPPED"]
    bad = ran[(ran.n_tumour < MIN_N) | (ran.n_normal < MIN_N)]
    print(f"4.2  comparisons run {len(ran)}, skipped "
          f"{int((log.status == 'SKIPPED').sum())}, "
          f"low_power {int((log.status == 'low_power').sum())}; "
          f"violating the n>={MIN_N} rule: {len(bad)}")
    if len(bad):
        failures.append(f"{len(bad)} comparisons ran below the D37 minimum")

    # 4.3 the p-value shape is power, not a batch effect
    checks += 4
    rows = []
    for (c, a), g in d.groupby(["cancer_type", "arm"]):
        rows.append(dict(n_eff=min(g.n_tumour.iloc[0], g.n_normal.iloc[0]),
                         frac_sig=(g.p_value < 0.01).mean(),
                         frac_pass=g.passes_default.mean()))
    r = pd.DataFrame(rows)
    rho_sig = spearmanr(r.n_eff, r.frac_sig).statistic
    rho_pass = spearmanr(r.n_eff, r.frac_pass).statistic
    med_lfc = float(d.log2fc.median())
    sig = d[d.fdr < 0.01]
    frac_up = float((sig.log2fc > 0).mean())
    print(f"4.3  significance vs sample size  rho = {rho_sig:.2f}  (power)")
    print(f"     after effect-size filter     rho = {rho_pass:.2f}  (dependence gone)")
    print(f"     global median log2FC = {med_lfc:+.3f}   "
          f"significant genes {frac_up:.0%} up / {1 - frac_up:.0%} down")
    if rho_sig < 0.5:
        failures.append(
            f"significant fraction does not track sample size (rho={rho_sig:.2f}); "
            "the heavy p<0.01 mass may be an artefact rather than power")
    if abs(rho_pass) > 0.6:
        failures.append(
            f"the effect-size threshold did not remove the power dependence "
            f"(rho={rho_pass:.2f})")
    if abs(med_lfc) > 0.25:
        failures.append(
            f"global median log2FC is {med_lfc:+.2f}, not ~0 — a one-directional "
            "shift of this size suggests a batch effect")
    if not 0.2 < frac_up < 0.8:
        failures.append(
            f"significant genes are {frac_up:.0%} up — a sound two-sided test on "
            "real signal gives both directions")

    # 4.4 FDR computed within each (cancer type, arm) group
    checks += 1
    grp = d[(d.cancer_type == d.cancer_type.iloc[0]) & (d.arm == d.arm.iloc[0])]
    from statsmodels.stats.multitest import multipletests
    recomputed = multipletests(grp.p_value.values, method="fdr_bh")[1]
    ok = np.allclose(np.sort(recomputed), np.sort(grp.fdr.values), atol=1e-9)
    print(f"4.4  FDR recomputed within one group matches: {ok}")
    if not ok:
        failures.append("FDR was not computed within each (cancer type, arm) group")

    # 4.5 arm concordance reported
    checks += 1
    wide = d.pivot_table(index=["ensembl_gene", "cancer_type"], columns="arm",
                         values="log2fc").dropna(subset=["gtex", "tcga_nat"])
    agree = float(((wide.gtex > 0) == (wide.tcga_nat > 0)).mean())
    rho_arm = spearmanr(wide.gtex, wide.tcga_nat).statistic
    print(f"4.5  arm concordance n={len(wide):,}  rho={rho_arm:.2f}  "
          f"sign agreement={agree:.1%}")
    if rho_arm < 0.2:
        failures.append(f"the two arms barely correlate (rho={rho_arm:.2f}) — "
                        "one of them is probably wrong")

    # 4.6 positive controls behave
    checks += len(POSITIVE)
    print("4.6  positive controls (GTEx arm):")
    for s in POSITIVE:
        g = d[(d.symbol == s) & (d.arm == "gtex")]
        if g.empty:
            failures.append(f"{s} absent from the results")
            print(f"       {s:<9} MISSING")
            continue
        n = int(g.passes_default.sum())
        best = g.loc[g.log2fc.idxmax()]
        print(f"       {s:<9} passes in {n:>2} cancer types; best "
              f"{best.cancer_type[:28]:<28} log2FC={best.log2fc:>5.2f}")
        if n == 0:
            failures.append(f"{s} passes in no cancer type")

    # 4.7 deterministic
    checks += 1
    before = (RESULTS / "step4_tumor_vs_normal.csv").read_bytes()
    subprocess.run([sys.executable, str(ROOT / "src" / "step4_differential.py")],
                   cwd=ROOT, capture_output=True, timeout=1800)
    same = (RESULTS / "step4_tumor_vs_normal.csv").read_bytes() == before
    print(f"4.7  rerun produces identical output: {same}")
    if not same:
        failures.append("step4_differential.py is not deterministic")

    print()
    if failures:
        print(f"FAILED: {len(failures)} problem(s) across {checks} checks\n")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"OK: {checks} checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
