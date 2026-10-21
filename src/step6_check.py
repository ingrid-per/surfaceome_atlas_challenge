"""Step 6 checks (register rows 6.1-6.6)."""
from __future__ import annotations
import subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.stats import spearmanr

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
MIN_N = 10


def main() -> int:
    p = pd.read_csv(RESULTS / "step6_protein_tumor_vs_normal.csv")
    m = pd.read_csv(RESULTS / "step6_rna_protein_concordance.csv")
    log = pd.read_csv(RESULTS / "step6_cohorts.csv")
    fails, checks = [], 0

    # 6.1 all 8 cohorts with matched normals ran (D12)
    checks += 1
    ok = log[log.status == "ok"]
    print(f"6.1  cohorts run {len(ok)} of 8 (D12: BRCA and GBM have no normals)")
    if len(ok) != 8:
        fails.append(f"only {len(ok)} of 8 CPTAC cohorts ran")

    # 6.2 no cohort below the sample minimum
    checks += 1
    bad = ok[(ok.n_tumour < MIN_N) | (ok.n_normal < MIN_N)]
    print(f"6.2  cohorts violating n>={MIN_N}: {len(bad)}")
    if len(bad):
        fails.append(f"{len(bad)} cohorts ran below the D37 minimum")

    # 6.3 run unbiased: more proteins tested than surface candidates (D48)
    checks += 1
    tot, cand = p.ensembl_gene.nunique(), p[p.is_surface_candidate].ensembl_gene.nunique()
    print(f"6.3  proteins tested {tot:,}, of which surface candidates {cand:,}")
    if tot <= cand:
        fails.append("Step 6 appears to have run only on surface candidates (D48)")

    # 6.4 no global directional shift (D51)
    checks += 2
    med = float(p.protein_log2fc.median())
    sig = p[p.protein_fdr < 0.01]
    up = float((sig.protein_log2fc > 0).mean())
    print(f"6.4  median protein log2FC {med:+.3f}; FDR<0.01 split {up:.0%} up")
    if abs(med) > 0.25:
        fails.append(f"global protein log2FC shift of {med:+.2f} suggests a batch effect")
    if not 0.3 < up < 0.7:
        fails.append(f"significant proteins are {up:.0%} up - not a balanced two-sided result")

    # 6.5 concordance measured against the matched arm (D49)
    checks += 1
    both = m.dropna(subset=["rna_log2fc_nat", "protein_log2fc"])
    rho = spearmanr(both.rna_log2fc_nat, both.protein_log2fc).statistic
    print(f"6.5  RNA/protein concordance on the matched arm: rho={rho:+.2f}, n={len(both):,}")
    if rho < 0.2:
        fails.append(f"RNA and protein barely correlate (rho={rho:.2f})")
    if "rna_log2fc_gtex" not in m.columns:
        fails.append("the GTEx arm was not carried alongside for reference (D49)")

    # 6.6 nothing was filtered on protein evidence (D48)
    checks += 1
    n_cand_rows = int(m.is_surface_candidate.sum())
    print(f"6.6  surface-candidate rows retained: {n_cand_rows:,}; "
          f"agreement categories: {m.loc[m.is_surface_candidate,'agreement'].nunique()}")
    if "rna_not_tested" not in set(m.loc[m.is_surface_candidate, "agreement"]):
        fails.append("rows without an RNA result were dropped rather than labelled")

    print()
    if fails:
        print(f"FAILED: {len(fails)} problem(s) across {checks} checks\n")
        for f in fails:
            print(f"  - {f}")
        return 1
    print(f"OK: {checks} checks passed")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
