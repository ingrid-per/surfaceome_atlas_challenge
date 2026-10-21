"""Step 3 checks (register rows 3.1-3.9).

The controls are the point. A method that cannot put EGFR in S1 and keep GAPDH
out has no business ranking novel targets.

Use --allow-partial while the DeepTMHMM run is still going: it skips the checks
that need the whole proteome but still runs every control assertion.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
PROC = ROOT / "data" / "processed"
RESULTS = ROOT / "results"

SPINE = 20_431
# The runbook's expected scale: SURFY predicts 2,886 and CSPA measured 1,492.
CAND_MIN, CAND_MAX = 1_500, 4_500

# (accession, symbol, required tier). None = any surface tier.
POSITIVE = [("P00533", "EGFR", "S1"), ("P04626", "ERBB2", "S1"),
            ("P09758", "TACSTD2", "S1"), ("Q96NY8", "NECTIN4", "S1"),
            ("P56747", "CLDN6", "S1"), ("P15391", "CD19", "S1")]
GPI_CONTROLS = [("Q13421", "MSLN", "S3"), ("P15328", "FOLR1", "S3")]
NEGATIVE = [("P04406", "GAPDH"), ("P60709", "ACTB"), ("P68431", "H3C1"),
            ("P62805", "H4C1"), ("P0DP23", "CALM1")]
# Not a candidate by design: surface-localised but with no anchor of its own.
# See proposal P6b — B2M is documented here rather than used as a Step 5 control.
EXCLUDED_BY_DESIGN = [("P61769", "B2M", "no TM segment and no GPI anchor; "
                       "reaches the surface via MHC class I")]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--allow-partial", action="store_true",
                    help="skip whole-proteome checks while the run is in flight")
    args = ap.parse_args()

    path = PROC / "surface_candidates.csv"
    if not path.exists():
        print("surface_candidates.csv missing — run step3_build_candidates.py")
        return 1
    df = pd.read_csv(path)
    by_acc = df.set_index("accession")
    failures, checks = [], 0

    # 3.1 every input protein accounted for
    checks += 1
    print(f"3.1  proteins in table: {len(df):,} (spine {SPINE:,})")
    if len(df) != SPINE:
        failures.append(f"table has {len(df):,} rows, expected {SPINE:,}")

    cand = df[df.is_candidate]
    counts = df.tier.value_counts()
    print(f"3.2  candidates: {len(cand):,}   " +
          "  ".join(f"{t}={int(counts.get(t, 0)):,}"
                    for t in ["S1", "S2", "S3", "Flagged"]))

    # 3.2 candidate count in the expected range (full run only)
    checks += 1
    if args.allow_partial:
        print("       (range check skipped: --allow-partial)")
    elif not (CAND_MIN <= len(cand) <= CAND_MAX):
        failures.append(
            f"candidate count {len(cand):,} outside {CAND_MIN:,}-{CAND_MAX:,}; "
            "investigate before Step 4")

    # 3.3 positive controls
    print("3.3  positive controls:")
    for acc, sym, want in POSITIVE:
        checks += 1
        if acc not in by_acc.index:
            failures.append(f"{sym} absent from the table")
            print(f"       {sym:<9} MISSING")
            continue
        r = by_acc.loc[acc]
        ok = r.tier == want
        print(f"       {sym:<9} tier={r.tier:<8} tm_sources={r.tm_sources} "
              f"extracellular={r.extracellular_reported:<6} "
              f"{'ok' if ok else 'EXPECTED ' + want}")
        if args.allow_partial and r.predictor_missing:
            continue   # cannot reach S1 until the predictor has seen it
        if not ok:
            failures.append(f"{sym}: tier {r.tier}, expected {want}")

    # 3.4 GPI controls reach S3 and carry no numeric length (D20)
    print("3.4  GPI controls:")
    for acc, sym, want in GPI_CONTROLS:
        checks += 2
        r = by_acc.loc[acc]
        print(f"       {sym:<9} tier={r.tier:<8} gpi_source={r.gpi_source:<9} "
              f"extracellular={r.extracellular_reported}")
        if r.tier != want:
            failures.append(f"{sym}: tier {r.tier}, expected {want}")
        if str(r.extracellular_reported) != "not_applicable":
            failures.append(
                f"{sym}: extracellular length is "
                f"{r.extracellular_reported!r}, must be not_applicable for a GPI protein (D20)")

    # 3.5 negative controls absent
    print("3.5  negative controls:")
    for acc, sym in NEGATIVE:
        checks += 1
        r = by_acc.loc[acc]
        print(f"       {sym:<9} tier={r.tier:<8} candidate={bool(r.is_candidate)}")
        if r.is_candidate:
            failures.append(f"{sym} is a surface candidate but must not be")

    # 3.6 excluded-by-design proteins, documented rather than silently absent
    print("3.6  excluded by design:")
    for acc, sym, why in EXCLUDED_BY_DESIGN:
        checks += 1
        r = by_acc.loc[acc]
        print(f"       {sym:<9} tier={r.tier:<8} — {why}")
        if r.is_candidate:
            failures.append(f"{sym} became a candidate; the surface filter changed")

    # 3.7 no S1 without two sources
    checks += 1
    bad = df[(df.tier == "S1") & ((df.tm_sources < 2) | df.predictor_missing)]
    print(f"3.7  S1 rows lacking two independent sources: {len(bad)}")
    if len(bad):
        failures.append(f"{len(bad)} S1 rows have fewer than two sources")

    # 3.8 GPI source recorded for every GPI call (D22)
    checks += 1
    gpi = df[df.gpi_any]
    missing_src = gpi[gpi.gpi_source.isna() | (gpi.gpi_source == "")]
    pred_only = int((df.gpi_source == "predicted").sum())
    print(f"3.8  GPI calls: {len(gpi):,}, predictor-only: {pred_only:,}, "
          f"without a recorded source: {len(missing_src)}")
    if len(missing_src):
        failures.append(f"{len(missing_src)} GPI calls carry no source (D22)")

    # 3.10 D32: a curated GPI anchor must not be overridden by a predicted TM
    # helix. This is the CEACAM5 regression — it sat in S2 on a predicted helix
    # UniProt does not annotate, which is the classic GPI-signal-as-TM error.
    checks += 1
    gpi_pred_only = df[df.gpi_any & (df.n_tm_segments == 0)]
    wrong = gpi_pred_only[gpi_pred_only.tier != "S3"]
    print(f"3.10 curated-GPI, no UniProt TM -> must be S3: "
          f"{len(gpi_pred_only)} proteins, {len(wrong)} misrouted")
    if len(wrong):
        failures.append(
            f"{len(wrong)} curated-GPI proteins not in S3 "
            f"(e.g. {list(wrong.symbol.head(5))}) — a predicted TM helix is "
            "overriding a curated anchor (D32)")
    for sym in ["CEACAM5", "CD48", "FCGR3B", "TNFRSF10C", "ART1"]:
        checks += 1
        r = df[df.symbol == sym]
        if r.empty or r.iloc[0].tier != "S3":
            failures.append(
                f"{sym} must be S3 (curated GPI, no UniProt TM) — got "
                f"{r.iloc[0].tier if not r.empty else 'absent'}")

    # 3.9 predictor coverage, or every gap flagged
    checks += 1
    miss = int(df.predictor_missing.sum())
    print(f"3.9  proteins without a predictor result: {miss:,}")
    if not args.allow_partial and miss:
        unflagged = df[df.predictor_missing & (df.tier == "S1")]
        if len(unflagged):
            failures.append(f"{len(unflagged)} proteins reached S1 without a prediction")
        print(f"       all {miss:,} are flagged and capped below S1")

    print()
    if failures:
        print(f"FAILED: {len(failures)} problem(s) across {checks} checks\n")
        for f in failures:
            print(f"  - {f}")
        return 1
    print(f"OK: {checks} checks passed"
          + ("  (partial: whole-proteome checks skipped)" if args.allow_partial else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
