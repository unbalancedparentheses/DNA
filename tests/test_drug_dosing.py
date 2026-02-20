"""Tests for drug-specific dosing recommendations."""

from genetic_health.drug_dosing import generate_drug_dosing


class TestGenerateDrugDosing:
    def _star_alleles(self, **overrides):
        base = {
            "CYP2C9": {"phenotype": "normal", "diplotype": "*1/*1"},
            "CYP2C19": {"phenotype": "normal", "diplotype": "*1/*1"},
            "CYP2D6": {"phenotype": "normal", "diplotype": "*1/*1"},
            "DPYD": {"phenotype": "normal", "diplotype": "*1/*1"},
            "TPMT": {"phenotype": "normal", "diplotype": "*1/*1"},
            "CYP3A5": {"phenotype": "normal", "diplotype": "*1/*1"},
            "NAT2": {"phenotype": "normal", "diplotype": "*1/*1"},
            "CYP1A2": {"phenotype": "normal", "diplotype": "*1/*1F"},
            "SLCO1B1": {"phenotype": "normal", "diplotype": "*1/*1"},
        }
        base.update(overrides)
        return base

    def test_empty_star_alleles(self):
        result = generate_drug_dosing({})
        assert result["recommendations"] == []
        assert "No pharmacogenomic data" in result["summary"]

    def test_none_star_alleles(self):
        result = generate_drug_dosing(None)
        assert result["recommendations"] == []

    def test_normal_metabolizers_still_get_recs(self):
        result = generate_drug_dosing(self._star_alleles())
        # CYP1A2 normal/rapid gets caffeine recommendation
        drugs = [r["drug"] for r in result["recommendations"]]
        assert "Caffeine" in drugs

    def test_cyp2c9_intermediate_warfarin(self):
        sa = self._star_alleles(CYP2C9={"phenotype": "intermediate", "diplotype": "*1/*2"})
        result = generate_drug_dosing(sa)
        warfarin_recs = [r for r in result["recommendations"] if r["drug"] == "Warfarin"]
        assert len(warfarin_recs) >= 1
        assert "reduce" in warfarin_recs[0]["action"].lower() or "25-50%" in warfarin_recs[0]["action"]

    def test_cyp2c19_poor_clopidogrel(self):
        sa = self._star_alleles(CYP2C19={"phenotype": "poor", "diplotype": "*2/*2"})
        result = generate_drug_dosing(sa)
        clop = [r for r in result["recommendations"] if r["drug"] == "Clopidogrel"]
        assert len(clop) >= 1
        assert "prasugrel" in clop[0]["action"].lower() or "ticagrelor" in clop[0]["action"].lower()

    def test_cyp2d6_poor_codeine(self):
        sa = self._star_alleles(CYP2D6={"phenotype": "poor", "diplotype": "*4/*4"})
        result = generate_drug_dosing(sa)
        codeine = [r for r in result["recommendations"] if r["drug"] == "Codeine"]
        assert len(codeine) >= 1
        assert "not" in codeine[0]["action"].lower() or "alternative" in codeine[0]["action"].lower()

    def test_dpyd_poor_is_critical(self):
        sa = self._star_alleles(DPYD={"phenotype": "poor", "diplotype": "*2A/*2A"})
        result = generate_drug_dosing(sa)
        fluoro = [r for r in result["recommendations"] if r["drug"] == "Fluoropyrimidines"]
        assert len(fluoro) >= 1
        assert len(result["warnings"]) >= 1
        assert "FATAL" in fluoro[0]["action"] or "CONTRAINDICATED" in fluoro[0]["action"]

    def test_slco1b1_intermediate_simvastatin(self):
        sa = self._star_alleles(SLCO1B1={"phenotype": "intermediate", "diplotype": "*1/*5"})
        result = generate_drug_dosing(sa)
        statin = [r for r in result["recommendations"] if r["drug"] == "Simvastatin"]
        assert len(statin) >= 1
        assert "rosuvastatin" in statin[0]["dose_guidance"].lower() or "20" in statin[0]["dose_guidance"]

    def test_summary_text(self):
        sa = self._star_alleles(CYP2D6={"phenotype": "poor", "diplotype": "*4/*4"})
        result = generate_drug_dosing(sa)
        assert "dosing recommendation" in result["summary"].lower()

    def test_result_structure(self):
        result = generate_drug_dosing(self._star_alleles())
        assert "recommendations" in result
        assert "warnings" in result
        assert "summary" in result
        for rec in result["recommendations"]:
            assert "drug" in rec
            assert "category" in rec
            assert "genes" in rec
            assert "action" in rec
            assert "dose_guidance" in rec
            assert "source" in rec
