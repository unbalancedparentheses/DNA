"""Tests for preventive care timeline generator."""

from genetic_health.preventive_care import generate_preventive_timeline


class TestGeneratePreventiveTimeline:
    MOCK_PRS_ELEVATED = {
        "hypertension": {
            "name": "Hypertension", "percentile": 90.0,
            "risk_category": "elevated", "snps_found": 15, "snps_total": 22,
        },
    }
    MOCK_PRS_HIGH = {
        "breast_cancer": {
            "name": "Breast Cancer", "percentile": 97.0,
            "risk_category": "high", "snps_found": 18, "snps_total": 22,
        },
    }
    MOCK_APOE_E4 = {
        "apoe_type": "e3/e4", "risk_level": "elevated",
        "alzheimer_or": 2.8, "confidence": "high",
    }
    MOCK_ACMG_BRCA = {
        "genes_screened": 81, "genes_with_variants": 1,
        "acmg_findings": [
            {"gene": "BRCA1", "traits": "Breast cancer"},
        ],
    }

    def test_empty_input(self):
        result = generate_preventive_timeline()
        assert "timeline" in result
        assert "summary" in result
        assert len(result["timeline"]) > 0  # Base screenings

    def test_base_screenings_present(self):
        result = generate_preventive_timeline()
        tests = [t["test"] for t in result["timeline"]]
        assert "Blood pressure" in tests
        assert "Lipid panel" in tests

    def test_prs_elevates_screening(self):
        result = generate_preventive_timeline(prs_results=self.MOCK_PRS_ELEVATED)
        assert result["early_screenings"] > 0
        bp = [t for t in result["timeline"]
              if "blood pressure" in t["test"].lower() and t["priority"] != "standard"]
        assert len(bp) >= 1

    def test_apoe_adds_lipid_screening(self):
        result = generate_preventive_timeline(apoe=self.MOCK_APOE_E4)
        lipid = [t for t in result["timeline"]
                 if "lipid" in t["test"].lower() and "APOE" in (t.get("genetic_basis") or "")]
        assert len(lipid) >= 1

    def test_acmg_brca_adds_screening(self):
        result = generate_preventive_timeline(acmg=self.MOCK_ACMG_BRCA)
        breast = [t for t in result["timeline"] if "breast" in t["test"].lower()]
        assert len(breast) >= 1
        assert any(t["priority"] == "urgent" for t in breast)

    def test_pharmacogenomic_card(self):
        star_alleles = {
            "CYP2C9": {"phenotype": "intermediate"},
            "CYP2C19": {"phenotype": "normal"},
        }
        result = generate_preventive_timeline(star_alleles=star_alleles)
        card = [t for t in result["timeline"] if "pharmacogenomic" in t["test"].lower()]
        assert len(card) == 1
        assert "CYP2C9" in card[0]["reason"]

    def test_no_pharma_card_if_all_normal(self):
        star_alleles = {
            "CYP2C9": {"phenotype": "normal"},
            "CYP2C19": {"phenotype": "normal"},
        }
        result = generate_preventive_timeline(star_alleles=star_alleles)
        card = [t for t in result["timeline"] if "pharmacogenomic" in t["test"].lower()]
        assert len(card) == 0

    def test_timeline_sorted_by_age(self):
        result = generate_preventive_timeline(
            prs_results=self.MOCK_PRS_ELEVATED,
            apoe=self.MOCK_APOE_E4,
        )
        ages = [t["start_age"] for t in result["timeline"]]
        assert ages == sorted(ages)

    def test_result_structure(self):
        result = generate_preventive_timeline()
        assert "timeline" in result
        assert "summary" in result
        assert "early_screenings" in result
        for t in result["timeline"]:
            assert "test" in t
            assert "start_age" in t
            assert "frequency" in t
            assert "reason" in t
            assert "priority" in t
