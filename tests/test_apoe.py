"""Tests for APOE haplotype calling module."""

from genetic_health.apoe import call_apoe_haplotype


def _make_genome(rsid_genotypes):
    """Helper to create genome dict."""
    return {
        rsid: {"chromosome": "19", "position": "100", "genotype": gt}
        for rsid, gt in rsid_genotypes.items()
    }


class TestCallApoeHaplotype:
    def test_e3_e3(self):
        genome = _make_genome({"rs429358": "TT", "rs7412": "CC"})
        result = call_apoe_haplotype(genome)
        assert result["apoe_type"] == "e3/e3"
        assert result["risk_level"] == "average"
        assert result["alzheimer_or"] == 1.0
        assert result["confidence"] == "high"

    def test_e3_e4(self):
        genome = _make_genome({"rs429358": "CT", "rs7412": "CC"})
        result = call_apoe_haplotype(genome)
        assert result["apoe_type"] == "e3/e4"
        assert result["risk_level"] == "elevated"
        assert result["alzheimer_or"] == 2.8

    def test_e4_e4(self):
        genome = _make_genome({"rs429358": "CC", "rs7412": "CC"})
        result = call_apoe_haplotype(genome)
        assert result["apoe_type"] == "e4/e4"
        assert result["risk_level"] == "high"
        assert result["alzheimer_or"] == 12.0

    def test_e2_e3(self):
        genome = _make_genome({"rs429358": "TT", "rs7412": "CT"})
        result = call_apoe_haplotype(genome)
        assert result["apoe_type"] == "e2/e3"
        assert result["risk_level"] == "reduced"

    def test_e2_e2(self):
        genome = _make_genome({"rs429358": "TT", "rs7412": "TT"})
        result = call_apoe_haplotype(genome)
        assert result["apoe_type"] == "e2/e2"
        assert result["risk_level"] == "reduced"

    def test_e2_e4(self):
        genome = _make_genome({"rs429358": "CT", "rs7412": "CT"})
        result = call_apoe_haplotype(genome)
        assert result["apoe_type"] == "e2/e4"
        assert result["risk_level"] == "moderate"

    def test_missing_snps_low_confidence(self):
        result = call_apoe_haplotype({})
        assert result["apoe_type"] == "Unknown"
        assert result["confidence"] == "low"
        assert result["alzheimer_or"] is None

    def test_partial_data_one_snp(self):
        genome = _make_genome({"rs429358": "TT"})
        result = call_apoe_haplotype(genome)
        assert result["apoe_type"] == "Unknown"
        assert result["confidence"] == "low"

    def test_result_keys(self):
        result = call_apoe_haplotype({})
        assert "apoe_type" in result
        assert "risk_level" in result
        assert "alzheimer_or" in result
        assert "description" in result
        assert "confidence" in result
        assert "details" in result

    def test_reversed_allele_order(self):
        """TC should work same as CT for rs429358."""
        genome = _make_genome({"rs429358": "TC", "rs7412": "CC"})
        result = call_apoe_haplotype(genome)
        assert result["apoe_type"] == "e3/e4"
