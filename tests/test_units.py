"""Unit tests for the pure functions the pipeline depends on.

These need no data files, so they run anywhere in under a second. Each one
covers a rule that, if it silently broke, would corrupt everything downstream
while still producing plausible-looking numbers — which is the dangerous kind
of bug in an analysis pipeline.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

from extra_mantis_export import clean_location
from step2_build_gene_map import strip_version
from step4_extract_matrix import XENA_FLOOR, to_log2_tpm1
from step5_specificity import log2p1, tau
from step7_rank import apply_family_limit


class TestStripVersion:
    """D17: every Ensembl id entering any join goes through this one function.

    Versions differ between GENCODE releases, so an unstripped id silently
    fails to join and the gene disappears rather than erroring.
    """

    @pytest.mark.parametrize("raw,want", [
        ("ENSG00000146648.22", "ENSG00000146648"),
        ("ENSG00000146648", "ENSG00000146648"),
        ("  ENSG00000146648.1  ", "ENSG00000146648"),
        ("ENSG00000223972.5", "ENSG00000223972"),
    ])
    def test_strips(self, raw, want):
        assert strip_version(raw) == want

    def test_is_idempotent(self):
        once = strip_version("ENSG00000146648.22")
        assert strip_version(once) == once


class TestToLog2Tpm1:
    """D11: Toil stores log2(TPM + 0.001), not TPM.

    Reading those values as TPM would make every downstream number wrong while
    still looking plausible, which is why this conversion lives in one tested
    function rather than being repeated at each call site.
    """

    def test_floor_maps_to_zero(self):
        # -9.9658 is the matrix floor and means "not expressed"
        assert to_log2_tpm1(np.array([XENA_FLOOR], dtype=np.float32))[0] == 0.0

    def test_never_negative(self):
        v = np.array([XENA_FLOOR, -5.0, -1.0, 0.0, 5.0], dtype=np.float32)
        assert (to_log2_tpm1(v) >= 0).all()

    def test_round_trip(self):
        # 100 TPM stored as log2(100.001) must come back as log2(101)
        stored = np.array([np.log2(100 + 0.001)], dtype=np.float32)
        assert to_log2_tpm1(stored)[0] == pytest.approx(np.log2(101), abs=1e-4)

    def test_monotonic(self):
        v = np.array([-9.9658, -2.0, 0.0, 3.0, 10.0], dtype=np.float32)
        out = to_log2_tpm1(v)
        assert (np.diff(out) >= 0).all()


class TestTau:
    """Tissue-specificity index: 0 = everywhere equally, 1 = one tissue only."""

    def test_single_tissue_is_one(self):
        assert tau(np.array([100.0, 0.0, 0.0, 0.0])) == pytest.approx(1.0)

    def test_uniform_is_zero(self):
        assert tau(np.array([50.0, 50.0, 50.0, 50.0])) == pytest.approx(0.0)

    def test_bounded(self):
        rng = np.random.default_rng(0)
        for _ in range(20):
            t = tau(rng.random(30) * 100)
            assert 0.0 <= t <= 1.0

    def test_all_zero_is_undefined(self):
        # a gene expressed nowhere has no meaningful specificity
        assert np.isnan(tau(np.array([0.0, 0.0, 0.0])))

    def test_ignores_missing(self):
        assert tau(np.array([100.0, np.nan, 0.0, 0.0])) == pytest.approx(1.0)


class TestLog2P1:
    def test_zero_maps_to_zero(self):
        assert log2p1(np.array([0.0]))[0] == 0.0

    def test_one_maps_to_one(self):
        assert log2p1(np.array([1.0]))[0] == pytest.approx(1.0)


class TestFamilyLimit:
    """D52: at most two genes per HGNC gene group, so one family cannot fill
    the list and hand over ten targets sharing a single failure mode."""

    @staticmethod
    def _frame(groups):
        return pd.DataFrame({"symbol": [f"G{i}" for i in range(len(groups))],
                             "gene_group": groups})

    def test_caps_at_two_per_family(self):
        out = apply_family_limit(self._frame(["Claudins"] * 5), 5)
        assert len(out) == 2

    def test_ungrouped_genes_are_never_capped(self):
        out = apply_family_limit(self._frame([np.nan] * 5), 5)
        assert len(out) == 5

    def test_preserves_input_order(self):
        out = apply_family_limit(self._frame(["A", "B", "A", "B", "C"]), 5)
        assert list(out.symbol) == ["G0", "G1", "G2", "G3", "G4"]

    def test_records_what_it_dropped(self):
        out = apply_family_limit(self._frame(["A", "A", "A", "B"]), 4)
        assert [s for s, _ in out.attrs["dropped_by_family"]] == ["G2"]

    def test_stops_at_n(self):
        assert len(apply_family_limit(self._frame(["A", "B", "C", "D"]), 2)) == 2


class TestCleanLocation:
    """UniProt location strings carry evidence codes that add tokens and no
    meaning to an embedding."""

    def test_strips_evidence_codes(self):
        s = clean_location(
            "SUBCELLULAR LOCATION: Cell membrane {ECO:0000250|UniProtKB:Q9Z262}")
        assert "ECO" not in s and "Cell membrane" in s

    def test_drops_the_note_block(self):
        assert "Note=" not in clean_location("Cell membrane. Note=Long prose here.")

    def test_handles_missing(self):
        assert clean_location(None) == ""
        assert clean_location(float("nan")) == ""
