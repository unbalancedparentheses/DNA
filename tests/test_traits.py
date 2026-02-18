"""Tests for trait prediction module."""

from genetic_health.traits import predict_traits


def _make_genome(rsid_genotypes):
    """Helper to create genome dict."""
    return {
        rsid: {"chromosome": "15", "position": "100", "genotype": gt}
        for rsid, gt in rsid_genotypes.items()
    }


class TestEyeColor:
    def test_blue_eyes_aa(self):
        genome = _make_genome({"rs12913832": "AA"})
        result = predict_traits(genome)
        assert "blue" in result["eye_color"]["prediction"].lower()
        assert result["eye_color"]["confidence"] == "high"

    def test_brown_eyes_gg(self):
        genome = _make_genome({"rs12913832": "GG"})
        result = predict_traits(genome)
        assert "brown" in result["eye_color"]["prediction"].lower()
        assert result["eye_color"]["confidence"] == "high"

    def test_green_hazel_ag(self):
        genome = _make_genome({"rs12913832": "AG"})
        result = predict_traits(genome)
        pred = result["eye_color"]["prediction"].lower()
        assert "green" in pred or "hazel" in pred

    def test_oca2_modifier(self):
        genome = _make_genome({"rs12913832": "AA", "rs1800407": "GA"})
        result = predict_traits(genome)
        pred = result["eye_color"]["prediction"].lower()
        assert "green" in pred or "blue" in pred

    def test_missing_snps(self):
        result = predict_traits({})
        assert result["eye_color"]["confidence"] == "low"


class TestHairColor:
    def test_red_hair(self):
        genome = _make_genome({"rs1805007": "TT"})
        result = predict_traits(genome)
        assert "red" in result["hair_color"]["prediction"].lower()
        assert result["hair_color"]["confidence"] == "high"

    def test_carrier_one_allele(self):
        genome = _make_genome({"rs1805007": "CT"})
        result = predict_traits(genome)
        pred = result["hair_color"]["prediction"].lower()
        assert "red" in pred or "auburn" in pred

    def test_no_red_variants(self):
        genome = _make_genome({"rs1805007": "CC", "rs1805008": "CC"})
        result = predict_traits(genome)
        assert "non-red" in result["hair_color"]["prediction"].lower()

    def test_missing_snps(self):
        result = predict_traits({})
        assert result["hair_color"]["confidence"] == "low"


class TestEarwax:
    def test_dry_earwax(self):
        genome = _make_genome({"rs17822931": "TT"})
        result = predict_traits(genome)
        assert "dry" in result["earwax_type"]["prediction"].lower()
        assert result["earwax_type"]["confidence"] == "high"

    def test_wet_earwax_cc(self):
        genome = _make_genome({"rs17822931": "CC"})
        result = predict_traits(genome)
        assert "wet" in result["earwax_type"]["prediction"].lower()

    def test_wet_earwax_ct(self):
        genome = _make_genome({"rs17822931": "CT"})
        result = predict_traits(genome)
        assert "wet" in result["earwax_type"]["prediction"].lower()

    def test_missing_snps(self):
        result = predict_traits({})
        assert result["earwax_type"]["confidence"] == "low"


class TestFreckling:
    def test_high_freckling(self):
        genome = _make_genome({"rs1805007": "TT"})
        result = predict_traits(genome)
        assert "high" in result["freckling"]["prediction"].lower()

    def test_moderate_freckling(self):
        genome = _make_genome({"rs1805007": "CT", "rs1805008": "CC"})
        result = predict_traits(genome)
        assert "moderate" in result["freckling"]["prediction"].lower()

    def test_typical(self):
        genome = _make_genome({"rs1805007": "CC", "rs1805008": "CC"})
        result = predict_traits(genome)
        assert "typical" in result["freckling"]["prediction"].lower()


class TestPredictTraits:
    def test_result_keys(self):
        result = predict_traits({})
        assert "eye_color" in result
        assert "hair_color" in result
        assert "earwax_type" in result
        assert "freckling" in result

    def test_each_trait_has_required_keys(self):
        result = predict_traits({})
        for trait_name, trait in result.items():
            assert "prediction" in trait, f"{trait_name} missing prediction"
            assert "confidence" in trait, f"{trait_name} missing confidence"
            assert "snps_used" in trait, f"{trait_name} missing snps_used"
            assert "description" in trait, f"{trait_name} missing description"
