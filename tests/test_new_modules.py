"""Tests for the 6 new analysis modules: pain, histamine, thyroid, hormone, eye, alcohol."""

from genetic_health.pain_sensitivity import profile_pain_sensitivity
from genetic_health.histamine import profile_histamine
from genetic_health.thyroid import profile_thyroid
from genetic_health.hormone_metabolism import profile_hormone_metabolism
from genetic_health.eye_health import profile_eye_health
from genetic_health.alcohol_profile import profile_alcohol


def _genome(**snps):
    """Build a genome dict from rsid=genotype kwargs."""
    return {rsid: {"genotype": gt, "chromosome": "1", "position": "1"}
            for rsid, gt in snps.items()}


# =========================================================================
# PAIN SENSITIVITY
# =========================================================================

class TestPainSensitivity:
    def test_empty_genome(self):
        r = profile_pain_sensitivity({})
        assert r["markers_found"] == 0
        assert r["pain_sensitivity_score"] == 50
        assert "summary" in r

    def test_comt_warrior(self):
        r = profile_pain_sensitivity(_genome(rs4680="GG"))
        assert r["markers_found"] >= 1
        assert r["pain_sensitivity_score"] < 50  # warrior = lower sensitivity

    def test_comt_worrier(self):
        r = profile_pain_sensitivity(_genome(rs4680="AA"))
        assert r["markers_found"] >= 1
        assert r["pain_sensitivity_score"] > 50  # worrier = higher sensitivity

    def test_oprm1_reduced(self):
        r = profile_pain_sensitivity(_genome(rs1799971="AG"))
        assert r["markers_found"] >= 1
        assert "domains" in r
        assert "opioid_response" in r["domains"]

    def test_full_profile(self):
        r = profile_pain_sensitivity(_genome(
            rs1799971="AA", rs4680="AG", rs6746030="GG", rs8065080="CC"))
        assert r["markers_found"] == 4
        assert 0 <= r["pain_sensitivity_score"] <= 100
        assert len(r["recommendations"]) > 0

    def test_result_structure(self):
        r = profile_pain_sensitivity(_genome(rs4680="GG"))
        assert "pain_sensitivity_score" in r
        assert "domains" in r
        assert "markers_found" in r
        assert "summary" in r
        assert "recommendations" in r


# =========================================================================
# HISTAMINE
# =========================================================================

class TestHistamine:
    def test_empty_genome(self):
        r = profile_histamine({})
        assert r["markers_found"] == 0
        assert r["markers_tested"] == 3
        assert "summary" in r

    def test_low_risk(self):
        r = profile_histamine(_genome(rs10156191="CC", rs1049793="GG", rs11558538="CC"))
        assert r["risk_level"] == "low"
        assert r["markers_found"] == 3

    def test_elevated_risk(self):
        r = profile_histamine(_genome(rs10156191="TT", rs1049793="CC", rs11558538="TT"))
        assert r["risk_level"] == "elevated"

    def test_moderate_risk(self):
        r = profile_histamine(_genome(rs10156191="CT", rs11558538="CT"))
        assert r["risk_level"] in ("moderate", "elevated")
        assert r["markers_found"] == 2

    def test_foods_to_watch(self):
        r = profile_histamine(_genome(rs10156191="TT"))
        assert isinstance(r["foods_to_watch"], list)

    def test_result_structure(self):
        r = profile_histamine(_genome(rs10156191="CC"))
        assert "risk_level" in r
        assert "markers_found" in r
        assert "markers_tested" in r
        assert "gene_results" in r
        assert "summary" in r
        assert "recommendations" in r
        assert "foods_to_watch" in r


# =========================================================================
# THYROID
# =========================================================================

class TestThyroid:
    def test_empty_genome(self):
        r = profile_thyroid({})
        assert r["markers_found"] == 0
        assert "summary" in r

    def test_foxe1_risk(self):
        r = profile_thyroid(_genome(rs965513="AA"))
        assert r["markers_found"] >= 1
        rp = r["risk_profile"]
        assert "cancer_risk" in rp

    def test_dio2_conversion(self):
        r = profile_thyroid(_genome(rs225014="TT"))
        assert r["markers_found"] >= 1
        rp = r["risk_profile"]
        assert "conversion_efficiency" in rp

    def test_autoimmune_risk(self):
        r = profile_thyroid(_genome(rs2071403="CC", rs179247="AA"))
        assert r["markers_found"] >= 2

    def test_full_profile(self):
        r = profile_thyroid(_genome(
            rs965513="AG", rs2071403="CT", rs179247="AG",
            rs11206244="CC", rs225014="CT"))
        assert r["markers_found"] == 5
        assert len(r["recommendations"]) > 0

    def test_result_structure(self):
        r = profile_thyroid(_genome(rs965513="GG"))
        assert "risk_profile" in r
        assert "markers_found" in r
        assert "summary" in r
        assert "recommendations" in r


# =========================================================================
# HORMONE METABOLISM
# =========================================================================

class TestHormoneMetabolism:
    def test_empty_genome(self):
        r = profile_hormone_metabolism({})
        assert r["markers_found"] == 0
        assert "summary" in r

    def test_aromatase(self):
        r = profile_hormone_metabolism(_genome(rs4646="TT"))
        assert r["markers_found"] >= 1
        assert "domains" in r

    def test_androgen(self):
        r = profile_hormone_metabolism(_genome(rs523349="GG"))
        assert r["markers_found"] >= 1

    def test_full_profile(self):
        r = profile_hormone_metabolism(_genome(
            rs4646="CT", rs2234693="CC", rs523349="GC", rs6152="AG"))
        assert r["markers_found"] == 4

    def test_result_structure(self):
        r = profile_hormone_metabolism(_genome(rs4646="CC"))
        assert "domains" in r
        assert "markers_found" in r
        assert "summary" in r
        assert "recommendations" in r


# =========================================================================
# EYE HEALTH
# =========================================================================

class TestEyeHealth:
    def test_empty_genome(self):
        r = profile_eye_health({})
        assert r["markers_found"] == 0
        assert "summary" in r

    def test_glaucoma_risk(self):
        r = profile_eye_health(_genome(rs4236601="AA"))
        assert r["markers_found"] >= 1
        assert "conditions" in r
        assert "glaucoma_risk" in r["conditions"]

    def test_myopia_risk(self):
        r = profile_eye_health(_genome(rs524952="AA", rs8027411="TT"))
        assert r["markers_found"] >= 2
        assert "myopia_risk" in r["conditions"]

    def test_myoc_rare_variant(self):
        r = profile_eye_health(_genome(rs74315329="AG"))
        assert r["markers_found"] >= 1
        # MYOC pathogenic variant should flag urgent
        conditions = r["conditions"]
        assert "glaucoma_risk" in conditions

    def test_result_structure(self):
        r = profile_eye_health(_genome(rs524952="GG"))
        assert "conditions" in r
        assert "markers_found" in r
        assert "summary" in r
        assert "recommendations" in r


# =========================================================================
# ALCOHOL PROFILE
# =========================================================================

class TestAlcoholProfile:
    def test_empty_genome(self):
        r = profile_alcohol({})
        assert r["markers_found"] == 0
        assert "summary" in r

    def test_normal_metabolism(self):
        r = profile_alcohol(_genome(rs1229984="CC", rs671="GG"))
        assert r["metabolism_speed"] == "normal"
        assert r["flush_risk"] == "none"
        assert r["cancer_risk"] == "average"

    def test_fast_metabolism(self):
        r = profile_alcohol(_genome(rs1229984="TT"))
        assert r["metabolism_speed"] == "fast"

    def test_flush_risk(self):
        r = profile_alcohol(_genome(rs671="GA"))
        assert r["flush_risk"] in ("mild", "severe")

    def test_severe_flush(self):
        r = profile_alcohol(_genome(rs671="AA"))
        assert r["flush_risk"] == "severe"
        assert r["cancer_risk"] in ("elevated", "high")

    def test_result_structure(self):
        r = profile_alcohol(_genome(rs1229984="CC"))
        assert "metabolism_speed" in r
        assert "flush_risk" in r
        assert "cancer_risk" in r
        assert "markers_found" in r
        assert "summary" in r
        assert "recommendations" in r
