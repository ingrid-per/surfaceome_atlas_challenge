"""Checks 9.1-9.6 for the robustness grid.

These re-derive the claims from the written outputs, independently of the
script that produced them.
"""

from __future__ import annotations

import itertools
import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
RESULTS = ROOT / "results"
sys.path.insert(0, str(Path(__file__).resolve().parent))
from step9_robustness import FDR, LOG2FC, MIN_EXTRA, SOURCE, SUMMARY, TOP_N

AXES = ["min_extra", "fdr", "log2fc", "source", "summary"]


def main() -> int:
    runs = pd.read_csv(RESULTS / "step9_runs.csv")
    st = pd.read_csv(RESULTS / "step9_stability.csv")
    top7 = pd.read_csv(RESULTS / "step7_top10.csv").symbol.tolist()
    fails = []

    def check(n: int, label: str, ok: bool, detail: str = "") -> None:
        print(f"  {n:<5} {label:<52} {'PASS' if ok else 'FAIL'}"
              f"{'  ' + detail if detail else ''}")
        if not ok:
            fails.append(n)

    grid_size = len(MIN_EXTRA) * len(FDR) * len(LOG2FC) * len(SOURCE) * len(SUMMARY)
    check(9.1, "run count equals the grid size", len(runs) == grid_size,
          f"{len(runs)} vs {grid_size}")

    # every cell of the grid is present exactly once
    want = {(a, float(b), float(c), d, e) for a, b, c, d, e
            in itertools.product(MIN_EXTRA, FDR, LOG2FC, SOURCE, SUMMARY)}
    have = {(r.min_extra, float(r.fdr), float(r.log2fc), r.source, r.summary)
            for r in runs.itertuples()}
    check(9.2, "every setting appears exactly once",
          have == want and not runs.duplicated(AXES).any())

    # the default cell is in the grid and its top 10 is Step 7's
    d = runs[runs.is_default]
    got = d.iloc[0].top20.split(";")[:10] if len(d) == 1 else []
    check(9.3, "default cell reproduces the Step 7 top 10", got == top7,
          "" if got == top7 else f"{got} != {top7}")

    check(9.4, "stability values lie in [0, 1]",
          bool(st.stability.between(0, 1).all()))

    # stability is recomputable from the run table
    counts: dict[str, int] = {}
    for t in runs.top20:
        for s in str(t).split(";"):
            counts[s] = counts.get(s, 0) + 1
    recomputed = pd.Series(counts).sort_index()
    stored = st.set_index("symbol").runs_in_top20.sort_index()
    check(9.5, "stability recomputes from step9_runs.csv",
          recomputed.equals(stored.astype(int)))

    # no run returned more than TOP_N, and none is empty
    sizes = runs.top20.map(lambda t: len(str(t).split(";")))
    check(9.6, f"every run returned 1..{TOP_N} genes",
          bool((sizes > 0).all() and (sizes <= TOP_N).all()),
          f"min {sizes.min()} max {sizes.max()}")

    print(f"\n  {len(runs)} runs, {len(st)} genes ever in the top {TOP_N}, "
          f"{int((st.stability >= .75).sum())} robust at >= 75%")
    if fails:
        print(f"  FAILED: {fails}")
    return 1 if fails else 0


if __name__ == "__main__":
    raise SystemExit(main())
