"""Tests for WGS pipeline utility functions.

Tests the pure-logic functions without requiring external tools
(minimap2, samtools, bcftools) or reference genome files.
"""

import json
import pytest
from pathlib import Path
from unittest.mock import patch

from genetic_health.wgs_pipeline import (
    cpu_count,
    _collect_all_rsids,
)


class TestCpuCount:
    def test_returns_positive_int(self):
        result = cpu_count()
        assert isinstance(result, int)
        assert result >= 1

    def test_fallback_on_failure(self):
        with patch('os.cpu_count', return_value=None):
            assert cpu_count() == 4


class TestCollectAllRsids:
    def test_returns_sorted_list(self):
        rsids = _collect_all_rsids()
        assert isinstance(rsids, list)
        assert rsids == sorted(rsids)

    def test_contains_known_snps(self):
        rsids = _collect_all_rsids()
        # From SNP database
        assert "rs1801133" in rsids  # MTHFR
        assert "rs4680" in rsids     # COMT
        # From PRS models
        assert "rs7903146" in rsids  # TCF7L2
        # From ancestry markers
        assert "rs1426654" in rsids  # SLC24A5
        # From star alleles
        assert "rs4244285" in rsids  # CYP2C19
        # From APOE
        assert "rs429358" in rsids
        assert "rs7412" in rsids
        # From blood type
        assert "rs505922" in rsids

    def test_has_hundreds_of_rsids(self):
        rsids = _collect_all_rsids()
        assert len(rsids) > 200

    def test_no_duplicates(self):
        rsids = _collect_all_rsids()
        assert len(rsids) == len(set(rsids))

    def test_all_start_with_rs(self):
        rsids = _collect_all_rsids()
        for rsid in rsids:
            assert rsid.startswith("rs"), f"Unexpected rsID format: {rsid}"
