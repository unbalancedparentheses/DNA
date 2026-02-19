"""Tests for carrier screening module."""

from genetic_health.carrier_screen import organize_carrier_findings


def _make_finding(gene, inheritance="recessive", is_het=True, is_hom=False):
    """Helper to create a ClinVar-like finding."""
    return {
        "gene": gene,
        "rsid": "rs12345",
        "chromosome": "1",
        "position": "100",
        "user_genotype": "AG",
        "ref": "A",
        "alt": "G",
        "traits": "Test condition; other trait",
        "gold_stars": 2,
        "is_homozygous": is_hom,
        "is_heterozygous": is_het,
        "inheritance": inheritance,
    }


class TestOrganizeCarrierFindings:
    def test_empty_input(self):
        result = organize_carrier_findings(None)
        assert result["total_carriers"] == 0
        assert result["carriers"] == []
        assert result["by_system"] == {}

    def test_empty_findings(self):
        result = organize_carrier_findings({
            "pathogenic": [],
            "likely_pathogenic": [],
        })
        assert result["total_carriers"] == 0

    def test_carrier_detected(self):
        findings = {
            "pathogenic": [_make_finding("CFTR", "recessive")],
            "likely_pathogenic": [],
        }
        result = organize_carrier_findings(findings)
        assert result["total_carriers"] == 1
        assert result["carriers"][0]["gene"] == "CFTR"
        assert result["carriers"][0]["condition"] == "Test condition"

    def test_homozygous_not_carrier(self):
        """Homozygous variants should not be classified as carriers."""
        findings = {
            "pathogenic": [_make_finding("CFTR", "recessive", is_het=False, is_hom=True)],
            "likely_pathogenic": [],
        }
        result = organize_carrier_findings(findings)
        assert result["total_carriers"] == 0

    def test_dominant_not_carrier(self):
        """Heterozygous dominant variants are AFFECTED, not carriers."""
        findings = {
            "pathogenic": [_make_finding("FBN1", "dominant")],
            "likely_pathogenic": [],
        }
        result = organize_carrier_findings(findings)
        assert result["total_carriers"] == 0

    def test_grouping_by_system(self):
        findings = {
            "pathogenic": [
                _make_finding("CFTR", "recessive"),
                _make_finding("HBB", "recessive"),
                _make_finding("HEXA", "recessive"),
            ],
            "likely_pathogenic": [],
        }
        result = organize_carrier_findings(findings)
        assert result["total_carriers"] == 3
        assert "Hematologic" in result["by_system"]
        assert "Metabolic" in result["by_system"]

    def test_reproductive_note_recessive(self):
        findings = {
            "pathogenic": [_make_finding("CFTR", "recessive")],
            "likely_pathogenic": [],
        }
        result = organize_carrier_findings(findings)
        assert "25%" in result["carriers"][0]["reproductive_note"]

    def test_couples_relevant(self):
        findings = {
            "pathogenic": [
                _make_finding("CFTR", "recessive"),
                _make_finding("HBB", "recessive"),
            ],
            "likely_pathogenic": [],
        }
        result = organize_carrier_findings(findings)
        assert len(result["couples_relevant"]) == 2

    def test_result_keys(self):
        result = organize_carrier_findings(None)
        assert "carriers" in result
        assert "total_carriers" in result
        assert "by_system" in result
        assert "couples_relevant" in result

    def test_likely_pathogenic_included(self):
        findings = {
            "pathogenic": [],
            "likely_pathogenic": [_make_finding("GBA", "recessive")],
        }
        result = organize_carrier_findings(findings)
        assert result["total_carriers"] == 1

    def test_cancer_genes_classified(self):
        """Cancer predisposition genes should be in Cancer Predisposition system."""
        findings = {
            "pathogenic": [_make_finding("BRCA1", "dominant")],
            "likely_pathogenic": [_make_finding("ATM", "recessive")],
        }
        result = organize_carrier_findings(findings)
        # BRCA1 is dominant so won't be a carrier, ATM recessive het = carrier
        carriers = result["carriers"]
        atm_carriers = [c for c in carriers if c["gene"] == "ATM"]
        assert len(atm_carriers) == 1
        assert atm_carriers[0]["system"] == "Cancer Predisposition"

    def test_cancer_genes_couples_relevant(self):
        """BRCA1/2 and Lynch genes should be flagged as couples-relevant."""
        findings = {
            "pathogenic": [_make_finding("MUTYH", "recessive")],
            "likely_pathogenic": [],
        }
        result = organize_carrier_findings(findings)
        assert result["total_carriers"] == 1
        assert result["carriers"][0]["couples_relevant"] is True
