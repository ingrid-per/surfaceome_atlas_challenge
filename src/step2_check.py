"""Step 2 checks (register rows 2.1-2.6).

The spot-check table and the NECTIN4 regression are the two that matter most:
the first proves well-known genes resolve, the second pins the specific bug a
symbol-based join would cause.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
MAP = ROOT / "data" / "processed" / "gene_map.csv"
RESULTS = ROOT / "results"

SPINE_SIZE = 20_431

# Verified by hand against UniProt, HPA and the GENCODE v23 probemap.
# (accession, current symbol, ensembl gene)
SPOT_CHECKS = [
    ("P00533", "EGFR", "ENSG00000146648"),
    ("P04626", "ERBB2", "ENSG00000141736"),
    ("Q13421", "MSLN", "ENSG00000102854"),
    ("P09758", "TACSTD2", "ENSG00000184292"),
    ("P04406", "GAPDH", "ENSG00000111640"),
    ("P60709", "ACTB", "ENSG00000075624"),
    ("P61769", "B2M", "ENSG00000166710"),
    ("P15328", "FOLR1", "ENSG00000110195"),
    ("Q96NY8", "NECTIN4", "ENSG00000143217"),
    ("P56747", "CLDN6", "ENSG00000184697"),
    ("P15391", "CD19", "ENSG00000177455"),
]

# Genes renamed since GENCODE v23. A symbol join silently drops these; an
# Ensembl-gene join keeps them. NECTIN4 is the one with an approved drug.
RENAMED = [("ENSG00000143217", "NECTIN4", "PVRL4"),
           ("ENSG00000110400", "NECTIN1", "PVRL1"),
           ("ENSG00000130202", "NECTIN2", "PVRL2"),
           ("ENSG00000177707", "NECTIN3", "PVRL3")]


def main() -> int:
    if not MAP.exists():
        print("data/processed/gene_map.csv missing — run step2_build_gene_map.py")
        return 1

    df = pd.read_csv(MAP, low_memory=False)
    failures: list[str] = []
    checks = 0

    # 2.1 coverage is reported and not catastrophic
    checks += 1
    n_acc = df["accession"].nunique()
    mapped = df[~df["no_ensembl"]]["accession"].nunique()
    pct = mapped / n_acc
    print(f"2.1  mapped {mapped:,}/{n_acc:,} accessions ({pct:.1%})")
    if n_acc != SPINE_SIZE:
        failures.append(f"spine is {n_acc:,}, expected {SPINE_SIZE:,}")
    if pct < 0.90:
        failures.append(f"only {pct:.1%} of accessions mapped — investigate before Step 3")

    # 2.2 every ambiguous row is logged
    checks += 1
    multi = pd.read_csv(RESULTS / "step2_multimapped.csv")
    logged = set(multi["accession"].astype(str))
    amb = df[df["ambiguous_mapping"]]
    unlogged = {a for a in amb["accession"] if a not in logged} - \
               {g for g in amb["ensembl_gene"].dropna() if g in logged}
    # an ambiguous row is explained either by its accession or by its gene
    unexplained = [a for a in amb["accession"].unique()
                   if a not in logged
                   and not set(amb[amb.accession == a]["ensembl_gene"]) & logged]
    print(f"2.2  {len(amb):,} ambiguous rows, {len(multi):,} logged cases, "
          f"{len(unexplained)} unexplained")
    if unexplained:
        failures.append(f"{len(unexplained)} ambiguous accessions absent from "
                        f"step2_multimapped.csv, e.g. {unexplained[:5]}")

    # 2.3 spot checks
    checks += len(SPOT_CHECKS)
    print("2.3  spot checks:")
    for acc, sym, gene in SPOT_CHECKS:
        hit = df[(df.accession == acc) & (df.ensembl_gene == gene)]
        if hit.empty:
            failures.append(f"spot check failed: {sym} {acc} -> {gene} not in map")
            print(f"       {sym:<9} MISSING")
            continue
        r = hit.iloc[0]
        ok = r.symbol == sym and bool(r.in_hpa) and bool(r.in_xena)
        print(f"       {sym:<9} {acc}  {gene}  hpa={r.hpa_symbol} "
              f"xena={r.xena_symbol}  {'ok' if ok else 'PROBLEM'}")
        if r.symbol != sym:
            failures.append(f"{acc}: symbol is {r.symbol}, expected {sym}")
        if not r.in_hpa:
            failures.append(f"{sym}: not matched to HPA")
        if not r.in_xena:
            failures.append(f"{sym}: not matched to the Xena matrix")

    # 2.4 the renamed-gene regression — the reason D17 exists
    checks += len(RENAMED)
    print("2.4  renamed-gene regression (symbol join would lose these):")
    for gene, current, old in RENAMED:
        hit = df[df.ensembl_gene == gene]
        if hit.empty:
            failures.append(f"{current}: {gene} absent from the map entirely")
            print(f"       {current:<9} MISSING")
            continue
        r = hit.iloc[0]
        if not r.in_xena:
            failures.append(f"{current} ({gene}) did not match the Xena matrix")
        if r.xena_symbol != old:
            failures.append(
                f"{current}: expected Xena symbol {old}, got {r.xena_symbol!r}")
        print(f"       {current:<9} {gene}  xena_symbol={r.xena_symbol!r}  "
              f"in_xena={bool(r.in_xena)}")

    # 2.5 nothing lost silently
    checks += 1
    unmapped = pd.read_csv(RESULTS / "step2_unmapped.csv")
    total = df[~df["no_ensembl"]]["accession"].nunique() + len(unmapped)
    print(f"2.5  mapped {mapped:,} + unmapped {len(unmapped):,} = {total:,} "
          f"(spine {SPINE_SIZE:,})")
    if total != SPINE_SIZE:
        failures.append(f"accessions unaccounted for: {total:,} != {SPINE_SIZE:,}")

    # 2.6 deterministic
    checks += 1
    before = MAP.read_bytes()
    subprocess.run([sys.executable, str(ROOT / "src" / "step2_build_gene_map.py")],
                   cwd=ROOT, capture_output=True, timeout=900)
    same = MAP.read_bytes() == before
    print(f"2.6  rerun produces identical output: {same}")
    if not same:
        failures.append("gene_map.csv is not deterministic across runs")

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
