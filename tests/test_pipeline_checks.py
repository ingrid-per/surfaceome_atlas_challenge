"""Integration tests: run the pipeline's own check scripts.

The project keeps one `*_check.py` per step, each asserting the claims that
step makes. Four of them read only files committed to the repository, so they
run anywhere — including CI — without the 1.7 GB of raw data. They are the
tests that would catch a result silently changing.

The remaining checks (steps 1-5) need `data/raw` or `data/processed`, which are
reproduced from `data/MANIFEST.md` by checksum rather than committed. They are
listed here and skipped when the data is absent, so the gap is visible rather
than quietly missing.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "src"

# checks whose inputs are all committed
COMMITTED = ["step6_check", "step7_check", "step9_check", "extra_metric_check"]

# checks that need data reproduced from the manifest first
NEEDS_DATA = {
    "step1_check": ROOT / "data" / "raw",
    "step2_check": ROOT / "data" / "processed" / "gene_map.csv",
    "step3_check": ROOT / "data" / "processed" / "surface_candidates.csv",
    "step4_check": ROOT / "data" / "processed" / "step4_candidate_tpm.parquet",
    "step5_check": ROOT / "data" / "processed" / "step4_candidate_tpm.parquet",
}


def run(script: str) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SRC / f"{script}.py")],
                          cwd=ROOT, capture_output=True, text=True, timeout=900)


@pytest.mark.parametrize("script", COMMITTED)
def test_check_passes(script):
    """Every assertion in the step's own check script holds."""
    r = run(script)
    assert r.returncode == 0, (
        f"{script} failed\n--- stdout ---\n{r.stdout[-2000:]}"
        f"\n--- stderr ---\n{r.stderr[-1000:]}")


@pytest.mark.parametrize("script,needed", sorted(NEEDS_DATA.items()))
def test_check_passes_when_data_present(script, needed):
    """Same, for checks that need data reproduced from the manifest."""
    if not needed.exists():
        pytest.skip(f"needs {needed.relative_to(ROOT)} — see data/MANIFEST.md")
    r = run(script)
    assert r.returncode == 0, f"{script} failed\n{r.stdout[-2000:]}"


def test_layout_is_intact():
    """The project structure and the rules file the pipeline depends on."""
    r = run("step0_verify_layout")
    assert r.returncode == 0, r.stdout[-1500:]


def test_top10_matches_the_published_payload():
    """The ranked result and the deliverable cannot drift apart.

    These are written by different scripts from the same source, so a change
    to one that does not reach the other would otherwise go unnoticed.
    """
    import json

    import pandas as pd
    top = pd.read_csv(ROOT / "results" / "step7_top10.csv")
    payload = json.loads((ROOT / "results" / "step7_targets.json").read_text())
    assert list(top.symbol) == [c["symbol"] for c in payload["top10"]]
    for row, card in zip(top.itertuples(), payload["top10"]):
        assert round(float(row.S), 2) == card["S"]
        assert row.cancer_type == card["cancer"]


def test_readme_figures_match_the_results():
    """Every finalist named in the README carries the S value the data gives."""
    import json

    readme = (ROOT / "README.md").read_text()
    payload = json.loads((ROOT / "results" / "step7_targets.json").read_text())
    for c in payload["top10"]:
        assert c["symbol"] in readme, f"{c['symbol']} missing from the README"
        assert f"{c['S']:.2f}" in readme, f"{c['symbol']} S value stale"
