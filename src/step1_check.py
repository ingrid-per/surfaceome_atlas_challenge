"""Step 1 checks (register rows 1.1-1.4).

Fails loudly rather than letting a missing file or an impossible count through
into Step 2.
"""

from __future__ import annotations

import csv
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from step1_resources import CPTAC_NO_NORMAL, RESOURCES

ROOT = Path(__file__).resolve().parent.parent
RAW = ROOT / "data" / "raw"
RESULTS = ROOT / "results"

# From the Xena Toil matrix header and phenotype table, verified before download.
EXPECT_SAMPLES = 19_131
EXPECT_GENES = 60_498
MIN_TCGA_TUMORS = 9_000
MIN_GTEX_NORMALS = 7_000


def main() -> int:
    failures: list[str] = []
    checks = 0

    # 1.1 every pinned resource is on disk and in the manifest with a checksum
    manifest = RESULTS / "step1_manifest.csv"
    if not manifest.exists():
        print("FAILED: results/step1_manifest.csv missing — run the downloader")
        return 1
    with manifest.open() as fh:
        rows = {r["filename"]: r for r in csv.DictReader(fh)}

    for res in RESOURCES:
        checks += 1
        path = RAW / res.filename
        if not path.exists():
            failures.append(f"missing file: {res.filename}")
            continue
        row = rows.get(res.filename)
        if row is None:
            failures.append(f"no manifest row: {res.filename}")
        elif len(row.get("sha256", "")) != 64:
            failures.append(f"no valid SHA-256 in manifest: {res.filename}")
        elif int(row["bytes"]) != path.stat().st_size:
            failures.append(
                f"size on disk differs from manifest: {res.filename}")

    # files in raw that nobody declared
    checks += 1
    declared = {r.filename for r in RESOURCES} | {".gitkeep"}
    for p in RAW.iterdir():
        if p.name not in declared and not p.name.endswith(".part"):
            failures.append(f"undeclared file in data/raw: {p.name}")

    # 1.2 counts are plausible
    counts = RESULTS / "step1_counts.txt"
    checks += 1
    if not counts.exists():
        failures.append("results/step1_counts.txt missing — run the inventory")
    else:
        text = counts.read_text()
        checks += 3
        if f"{EXPECT_GENES:,} genes x {EXPECT_SAMPLES:,} samples" not in text:
            failures.append(
                f"Toil matrix is not {EXPECT_GENES:,} x {EXPECT_SAMPLES:,} "
                "— the release changed, so re-pin it in the manifest")
        for label, floor in (("TCGA tumor", MIN_TCGA_TUMORS),
                             ("GTEx normal", MIN_GTEX_NORMALS)):
            line = next((l for l in text.splitlines() if "included by D10" in l), "")
            if label not in line:
                failures.append(f"no '{label}' count in step1_counts.txt")
        # no group may be zero
        if " 0  " in text and "no matched normal" not in text:
            failures.append("a sample-count group is zero")

    # 1.3 CPTAC normals present exactly where expected
    prot = "proteomics_gene_abundance_log2_reference_intensity_normalized"
    for res in RESOURCES:
        if res.group != "CPTAC" or not res.filename.endswith("_Normal.txt"):
            continue
        checks += 1
        cohort = res.filename.split("_")[0]
        if cohort in CPTAC_NO_NORMAL:
            failures.append(f"{cohort} is listed as having no normal but one was pinned")

    # 1.5 data/MANIFEST.md itself lists every file on disk.
    #
    # Check 1.1 verifies the CSV. The rule in CLAUDE.md names the MARKDOWN
    # file, and the two are written by the same function from the same list —
    # so a file recorded by hand in the markdown only is invisible to that
    # function and disappears the next time any download runs. That happened to
    # aran2015_tumour_purity.xlsx on 2026-10-22: it sat in data/raw with no
    # manifest row until this check was written.
    checks += 1
    md = (ROOT / "data" / "MANIFEST.md").read_text()
    for p_raw in sorted(RAW.iterdir()):
        if p_raw.name == ".gitkeep" or p_raw.name.endswith(".part"):
            continue
        if f"`{p_raw.name}`" not in md:
            failures.append(f"file on disk but not in data/MANIFEST.md: {p_raw.name}")

    # 1.4 data/raw is read-only
    checks += 1
    if os.access(RAW, os.W_OK):
        failures.append("data/raw is writable — raw data must be read-only")
    for p in RAW.iterdir():
        if p.is_file() and p.name != ".gitkeep" and os.access(p, os.W_OK):
            failures.append(f"writable raw file: {p.name}")
            break

    if failures:
        print(f"FAILED: {len(failures)} problem(s) across {checks} checks\n")
        for f in failures:
            print(f"  - {f}")
        return 1

    total = sum(int(r["bytes"]) for r in rows.values())
    print(f"OK: {checks} checks passed")
    print(f"  {len(rows)} files, {total / 1048576:.1f} MB, all with SHA-256")
    print(f"  Toil matrix {EXPECT_GENES:,} genes x {EXPECT_SAMPLES:,} samples")
    print("  data/raw is read-only")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
