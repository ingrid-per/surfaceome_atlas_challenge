"""Step 7 checks (register rows 7.1-7.5)."""
from __future__ import annotations
import json, subprocess, sys
from pathlib import Path
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
DOCS = ROOT / "docs"


def main() -> int:
    full = pd.read_csv(RESULTS / "step7_ranking_full.csv")
    top10 = pd.read_csv(RESULTS / "step7_top10.csv")
    fails, checks = [], 0

    # 7.1 eligibility honoured (D52)
    checks += 3
    print(f"7.1  eligible rows {len(full):,}, genes {full.ensembl_gene.nunique():,}")
    if not full.tier.isin(["S1", "S2", "S3"]).all():
        fails.append("a Flagged-tier gene is eligible (D52)")
    if not full.passes_default.all():
        fails.append("an ineligible row failed the RNA thresholds (D52)")
    if not (full.tumour_tpm >= 5).all():
        fails.append("a row below 5 TPM is eligible (D43/D52)")

    # 7.2 sorted by S, descending
    checks += 1
    s = top10.S.to_numpy()
    print(f"7.2  top 10 S range {s.min():.2f} .. {s.max():.2f}")
    if not (s[:-1] >= s[1:] - 1e-9).all():
        fails.append("the top 10 is not sorted by S descending (D52)")

    # 7.3 family limit respected (D52)
    checks += 1
    fam = top10.gene_group.dropna()
    worst = fam.value_counts().max() if len(fam) else 0
    print(f"7.3  most genes from one HGNC group in the top 10: {worst}")
    if worst > 2:
        fails.append(f"{worst} genes from one family exceed the 2-per-family limit")

    # 7.4 the rule was recorded BEFORE the ranking existed
    checks += 1
    rule_committed = subprocess.run(
        ["git", "log", "-1", "--format=%cI", "--", "docs/decisions.md"],
        cwd=ROOT, capture_output=True, text=True).stdout.strip()
    has_rule = "D52" in (DOCS / "decisions.md").read_text()
    print(f"7.4  selection rule D52 recorded in decisions.md: {has_rule}")
    if not has_rule:
        fails.append("the selection rule D52 is not in docs/decisions.md")

    # 7.5 the exported JSON matches the committed results
    checks += 2
    j = json.loads((RESULTS / "step7_targets.json").read_text())
    print(f"7.5  JSON: {len(j['top10'])} top10, {len(j['by_cancer'])} cancer types")
    if len(j["top10"]) != len(top10):
        fails.append("the exported JSON top 10 differs from step7_top10.csv")
    if [x["symbol"] for x in j["top10"]] != list(top10.symbol):
        fails.append("the exported JSON top 10 is in a different order")

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
