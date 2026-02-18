"""Tests for blood type prediction module."""

from genetic_health.blood_type import predict_blood_type


def _make_genome(rsid_genotypes):
    """Helper to create genome dict."""
    return {
        rsid: {"chromosome": "9", "position": "100", "genotype": gt}
        for rsid, gt in rsid_genotypes.items()
    }


class TestPredictBloodType:
    def test_type_o(self):
        genome = _make_genome({
            "rs505922": "TT",
            "rs8176746": "CC",
            "rs590787": "CC",
        })
        result = predict_blood_type(genome)
        assert result["abo"] == "O"
        assert result["rh"] == "+"
        assert result["blood_type"] == "O+"

    def test_type_a(self):
        genome = _make_genome({
            "rs505922": "CT",
            "rs8176746": "CC",
            "rs590787": "CC",
        })
        result = predict_blood_type(genome)
        assert result["abo"] == "A"
        assert result["rh"] == "+"

    def test_type_b(self):
        genome = _make_genome({
            "rs505922": "CT",
            "rs8176746": "TT",
            "rs590787": "CC",
        })
        result = predict_blood_type(genome)
        assert result["abo"] == "B"
        assert result["rh"] == "+"

    def test_type_ab(self):
        genome = _make_genome({
            "rs505922": "CC",
            "rs8176746": "CT",
            "rs590787": "CC",
        })
        result = predict_blood_type(genome)
        assert result["abo"] == "AB"
        assert result["rh"] == "+"
        assert result["blood_type"] == "AB+"

    def test_rh_negative(self):
        genome = _make_genome({
            "rs505922": "TT",
            "rs8176746": "CC",
            "rs590787": "TT",
        })
        result = predict_blood_type(genome)
        assert result["rh"] == "-"
        assert result["blood_type"] == "O-"

    def test_missing_snps_low_confidence(self):
        result = predict_blood_type({})
        assert result["confidence"] == "low"
        assert result["blood_type"] == "Unknown"

    def test_partial_data(self):
        genome = _make_genome({"rs505922": "TT"})
        result = predict_blood_type(genome)
        assert result["abo"] == "O"
        assert result["rh"] == "Unknown"
        assert result["confidence"] == "moderate"

    def test_result_keys(self):
        result = predict_blood_type({})
        assert "blood_type" in result
        assert "abo" in result
        assert "rh" in result
        assert "confidence" in result
        assert "details" in result

    def test_cc_no_b_gives_a(self):
        genome = _make_genome({
            "rs505922": "CC",
            "rs8176746": "CC",
        })
        result = predict_blood_type(genome)
        assert result["abo"] == "A"
