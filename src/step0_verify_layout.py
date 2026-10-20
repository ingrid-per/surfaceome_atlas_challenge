"""Step 0 check: is the project scaffold complete and are the rules actually present?

Re-runnable in any later session as a regression test — if someone deletes a rule
from CLAUDE.md or starts committing raw data, this fails.

Exits non-zero on the first set of failures, listing all of them.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

REQUIRED_DIRS = ["data/raw", "data/processed", "src", "results", "docs"]

REQUIRED_FILES = [
    "CLAUDE.md",
    ".gitignore",
    "environment.yml",
    "docs/decisions.md",
    "docs/checks.md",
    "data/MANIFEST.md",
    "Surface Target Pipeline Claude Code Runbook.md",
]

# Distinctive fragments of the eight rules from the runbook's "File 1" block. Matching
# on a fragment rather than the whole line survives trivial rewording but still
# catches a rule being dropped.
REQUIRED_RULE_FRAGMENTS = [
    "must come from a script in src/",
    "Never pick a threshold, gene list or tissue list silently",
    "version or date, source URL and checksum",
    "rows in, rows out",
    "automated check for every step",
    "data/raw is read-only",
    "too good or too empty",
    "plain-English comments",
]


def main() -> int:
    failures: list[str] = []

    for d in REQUIRED_DIRS:
        if not (ROOT / d).is_dir():
            failures.append(f"missing directory: {d}")

    for f in REQUIRED_FILES:
        p = ROOT / f
        if not p.is_file():
            failures.append(f"missing file: {f}")
        elif p.stat().st_size == 0:
            failures.append(f"empty file: {f}")

    claude = ROOT / "CLAUDE.md"
    if claude.is_file():
        text = claude.read_text()
        missing = [r for r in REQUIRED_RULE_FRAGMENTS if r not in text]
        if missing:
            failures.append(
                f"CLAUDE.md is missing {len(missing)} of "
                f"{len(REQUIRED_RULE_FRAGMENTS)} rules: {missing}"
            )

    # data/raw must be excluded from git, or the first download bloats the repo.
    probe = "data/raw/probe.csv"
    try:
        res = subprocess.run(
            ["git", "check-ignore", "-q", probe],
            cwd=ROOT, capture_output=True, timeout=30,
        )
        if res.returncode != 0:
            failures.append(f".gitignore does not exclude {probe}")
    except FileNotFoundError:
        failures.append("git not available, cannot verify .gitignore")

    checked = len(REQUIRED_DIRS) + len(REQUIRED_FILES) + len(REQUIRED_RULE_FRAGMENTS) + 1
    if failures:
        print(f"FAILED: {len(failures)} problem(s) across {checked} checks\n")
        for f in failures:
            print(f"  - {f}")
        return 1

    print(f"OK: project layout complete, {checked} checks passed")
    print(f"  {len(REQUIRED_DIRS)} directories, {len(REQUIRED_FILES)} files, "
          f"{len(REQUIRED_RULE_FRAGMENTS)} rules present in CLAUDE.md, "
          "data/raw excluded from git")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
