"""Tests for quality metrics module."""

import tempfile
from pathlib import Path

from genetic_health.quality_metrics import compute_quality_metrics


def _make_genome(entries):
    """Helper to create genome dict from list of (rsid, chrom, pos, gt)."""
    return {
        rsid: {"chromosome": chrom, "position": pos, "genotype": gt}
        for rsid, chrom, pos, gt in entries
    }


class TestComputeQualityMetrics:
    def test_empty_genome(self):
        result = compute_quality_metrics({})
        assert result["total_snps"] == 0
        assert result["autosomal_count"] == 0
        assert result["het_rate"] == 0.0
        assert result["has_mt"] is False
        assert result["has_y"] is False

    def test_basic_counts(self):
        genome = _make_genome([
            ("rs1", "1", "100", "AA"),
            ("rs2", "1", "200", "AG"),
            ("rs3", "2", "300", "CC"),
        ])
        result = compute_quality_metrics(genome)
        assert result["total_snps"] == 3
        assert result["autosomal_count"] == 3
        assert result["chromosomes"]["1"] == 2
        assert result["chromosomes"]["2"] == 1

    def test_heterozygosity_rate(self):
        genome = _make_genome([
            ("rs1", "1", "100", "AG"),
            ("rs2", "1", "200", "AG"),
            ("rs3", "2", "300", "AA"),
            ("rs4", "2", "400", "CC"),
        ])
        result = compute_quality_metrics(genome)
        assert result["het_rate"] == 0.5  # 2 het out of 4 autosomal

    def test_mt_detection(self):
        genome = _make_genome([
            ("rs1", "1", "100", "AA"),
            ("rs2", "MT", "200", "A"),
        ])
        result = compute_quality_metrics(genome)
        assert result["has_mt"] is True
        assert result["mt_snp_count"] == 1

    def test_y_detection(self):
        genome = _make_genome([
            ("rs1", "1", "100", "AA"),
            ("rs2", "Y", "200", "G"),
        ])
        result = compute_quality_metrics(genome)
        assert result["has_y"] is True

    def test_call_rate_with_file(self):
        genome = _make_genome([
            ("rs1", "1", "100", "AA"),
            ("rs2", "1", "200", "AG"),
        ])
        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as f:
            f.write("# header\n")
            f.write("rs1\t1\t100\tAA\n")
            f.write("rs2\t1\t200\tAG\n")
            f.write("rs3\t1\t300\t--\n")
            f.write("rs4\t1\t400\t--\n")
            tmp_path = Path(f.name)

        result = compute_quality_metrics(genome, tmp_path)
        assert result["no_call_count"] == 2
        assert result["call_rate"] == 2 / (2 + 2)
        tmp_path.unlink()

    def test_call_rate_no_file(self):
        genome = _make_genome([("rs1", "1", "100", "AA")])
        result = compute_quality_metrics(genome)
        assert result["no_call_count"] == 0
        assert result["call_rate"] == 1.0

    def test_chromosome_counting(self):
        genome = _make_genome([
            ("rs1", "1", "100", "AA"),
            ("rs2", "X", "200", "AG"),
            ("rs3", "22", "300", "CC"),
            ("rs4", "MT", "400", "T"),
        ])
        result = compute_quality_metrics(genome)
        assert result["chromosomes"]["1"] == 1
        assert result["chromosomes"]["X"] == 1
        assert result["chromosomes"]["22"] == 1
        assert result["chromosomes"]["MT"] == 1
        assert result["autosomal_count"] == 2  # chr 1 and 22

    def test_none_chromosome_skipped(self):
        genome = _make_genome([("rs1", "1", "100", "AA")])
        genome["rs2"] = {"chromosome": None, "position": "200", "genotype": "AG"}
        result = compute_quality_metrics(genome)
        assert "NONE" not in result["chromosomes"]
        assert result["chromosomes"]["1"] == 1
        assert result["total_snps"] == 2  # still counted in total

    def test_empty_chromosome_skipped(self):
        genome = _make_genome([("rs1", "1", "100", "AA")])
        genome["rs2"] = {"chromosome": "", "position": "200", "genotype": "AG"}
        result = compute_quality_metrics(genome)
        assert "" not in result["chromosomes"]
        assert result["chromosomes"]["1"] == 1

    def test_result_keys(self):
        result = compute_quality_metrics({})
        expected_keys = {
            "total_snps", "no_call_count", "call_rate", "chromosomes",
            "has_mt", "has_y", "mt_snp_count", "autosomal_count", "het_rate",
        }
        assert set(result.keys()) == expected_keys
