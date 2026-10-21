"""Checks for the metric-sensitivity extra.

Same shape as src/step9_check.py: re-derive the claims from the written files,
independently of the script that produced them.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from extra_metric_sensitivity import GATES, REFS, TOP_N


def main() -> int:
    runs = pd.read_csv(RESULTS / "extra_metric_runs.csv")
    st = pd.read_csv(RESULTS / "extra_metric_stability.csv")
    top7 = pd.read_csv(RESULTS / "step7_top10.csv").symbol.tolist()
    fails = []

    def check(n, label, ok, detail=""):
        print(f"  {n:<6} {label:<54} {'PASS' if ok else 'FAIL'}"
              f"{'  ' + detail if detail else ''}")
        if not ok:
            fails.append(n)

    expected = (len(REFS)                      # M1
                + len(REFS) * len(GATES)       # M2
                + len(REFS)                    # M3
                + 2                            # M4, GTEx and HPA only
                + len(REFS)                    # M5
                + 1)                           # M6
    check("E.1", "ranking count matches the specification",
          len(runs) == expected, f"{len(runs)} vs {expected}")

    check("E.2", "no duplicated (metric, reference) pair",
          not runs.duplicated(["metric", "reference"]).any())

    # the baseline metric must reproduce Step 7 — the same guarantee as 9.1
    base = runs[(runs.metric == "M1 S") & (runs.reference == "both")]
    got = base.iloc[0].top20.split(";")[:10] if len(base) == 1 else []
    check("E.3", "M1 S (both) reproduces the Step 7 top 10", got == top7,
          "" if got == top7 else f"{got[:3]}... != {top7[:3]}...")

    check("E.4", "stability values lie in [0, 1]",
          bool(st.stability.between(0, 1).all()))

    counts = {}
    for t in runs.top20:
        for s in str(t).split(";"):
            counts[s] = counts.get(s, 0) + 1
    rec = pd.Series(counts).sort_index()
    stored = st.set_index("symbol").rankings_in_top20.sort_index()
    check("E.5", "stability recomputes from extra_metric_runs.csv",
          rec.equals(stored.astype(int)))

    sizes = runs.top20.map(lambda t: len(str(t).split(";")))
    check("E.6", f"every ranking returned 1..{TOP_N} genes",
          bool((sizes > 0).all() and (sizes <= TOP_N).all()),
          f"min {sizes.min()} max {sizes.max()}")

    # a gate that removes nothing would make M2 a duplicate of M1
    for r in REFS:
        m1 = runs[(runs.metric == "M1 S") & (runs.reference == r)].iloc[0].top20
        strictest = runs[(runs.metric == f"M2 S, prevalence>={max(GATES):.2f}")
                         & (runs.reference == r)].iloc[0].top20
        check(f"E.7{r[0]}", f"the strictest prevalence gate changes the {r} list",
              m1 != strictest)

    print(f"\n  {len(runs)} rankings, {len(st)} genes ever in the top {TOP_N}")
    if fails:
        print(f"  FAILED: {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
