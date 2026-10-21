"""Step 5 checks (register rows 5.1-5.6).

Check 5.3 is the one worth reading. The runbook asks that a ubiquitous surface
protein score near zero. It originally named B2M, which our surface filter
correctly excludes (no anchor of its own), so P6b replaced it with a control
chosen by rule. That control must score low, and the known targets must NOT be
required to pass — D45 removed the pass/fail precisely because S > 0 rejects
every approved drug target.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

MIN_TUMOUR_TPM = 5.0

# Known targets. Used to CHARACTERISE the score, never to tune it (D45): the
# threshold was deliberately not calibrated on these, so Step 8 stays honest.
KNOWN = ["EGFR", "ERBB2", "MSLN", "FOLR1", "CEACAM5", "TACSTD2", "NECTIN4",
         "CLDN6", "CD19", "CA9", "GPC3"]


def main() -> int:
    d = pd.read_csv(RESULTS / "step5_specificity.csv")
    h = pd.read_csv(RESULTS / "step5_healthy_profile.csv", index_col=0)
    failures, checks = [], 0

    # 5.1 every carried-forward row has a healthy maximum and a named tissue
    checks += 3
    cf = d[d.carried_forward]
    missing = int(cf.healthy_max_log2.isna().sum())
    no_tissue = int(cf.healthy_max_tissue.isna().sum())
    print(f"5.1  carried forward {len(cf):,} rows; "
          f"missing healthy max {missing}; missing tissue name {no_tissue}")
    if missing or no_tissue:
        failures.append("some carried-forward rows have no healthy maximum")
    if not (cf.tumour_tpm >= MIN_TUMOUR_TPM).all():
        failures.append(f"a carried-forward row is below {MIN_TUMOUR_TPM} TPM (D43)")
    if (d.healthy_max_log2 < 0).any():
        failures.append("negative healthy expression")

    # 5.2 the maximum is never below the 90th percentile, by construction
    checks += 1
    bad = int((d.S > d.S_p90_crosscheck + 1e-9).sum())
    print(f"5.2  rows where S exceeds S_p90 (impossible): {bad}")
    if bad:
        failures.append(f"{bad} rows have S > S_p90; the maximum is not the maximum")

    # 5.3 the ubiquitous control scores low (P6b)
    checks += 1
    cand = pd.read_csv(PROC / "surface_candidates.csv")
    pm = cand[cand.is_candidate].groupby("ensembl_gene").agg(
        pm=("plasma_membrane_evidence", "any"))
    ctrl = h.join(d.drop_duplicates("ensembl_gene")
                  .set_index("ensembl_gene")[["symbol", "tier"]], how="inner")
    ctrl = ctrl.join(pm, how="left")
    n_tis = int(ctrl.gtex_n_tissues.max())
    rule = ctrl[(ctrl.tier == "S1") & ctrl.pm.fillna(False)
                & (ctrl.gtex_n_tissues_above_10tpm >= n_tis)].sort_values("gtex_tau")
    print(f"5.3  P6b control set by rule: {len(rule)} genes, "
          f"top = {rule.symbol.iloc[0] if len(rule) else 'NONE'}")
    if rule.empty:
        failures.append("the P6b rule selected no control gene")
    else:
        top = rule.index[:5]
        scores = d[d.ensembl_gene.isin(top)].groupby("symbol").S.max()
        for sym, s in scores.items():
            print(f"       {sym:<10} best S across cancer types = {s:+.2f}")
        worst = float(scores.max()) if len(scores) else np.nan
        if worst > 2.0:
            failures.append(
                f"a ubiquitous control scores S = {worst:+.2f}; a protein "
                "present in every tissue should not look tumour-specific")

    # 5.4 known targets are characterised, NOT required to pass (D45)
    checks += 1
    print("5.4  known targets (characterisation only, no pass/fail — D45):")
    found = 0
    for sym in KNOWN:
        g = d[d.symbol == sym]
        if g.empty:
            continue
        found += 1
        b = g.loc[g.S.idxmax()]
        print(f"       {sym:<9} S={b.S:+6.2f}  S_relaxed={b.S_relaxed:+6.2f}  "
              f"worst={str(b.healthy_max_tissue)[:22]:<22} "
              f"({b.cancer_type[:24]})")
    if found < len(KNOWN) - 2:
        failures.append(f"only {found} of {len(KNOWN)} known targets found")

    # 5.5 the tissue-exclusion variants are ordered as they must be
    checks += 1
    viol = int(((d.S_nobrain < d.S - 1e-9)
                | (d.S_relaxed < d.S_nobrain - 1e-9)).sum())
    print(f"5.5  variant ordering S <= S_nobrain <= S_relaxed violated in {viol} rows")
    if viol:
        failures.append(f"{viol} rows violate the S variant ordering (D46)")

    # 5.6 deterministic
    checks += 1
    before = (RESULTS / "step5_specificity.csv").read_bytes()
    subprocess.run([sys.executable, str(ROOT / "src" / "step5_specificity.py")],
                   cwd=ROOT, capture_output=True, timeout=1800)
    same = (RESULTS / "step5_specificity.csv").read_bytes() == before
    print(f"5.6  rerun produces identical output: {same}")
    if not same:
        failures.append("step5_specificity.py is not deterministic")

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
