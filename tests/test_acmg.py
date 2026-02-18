"""Tests for ACMG secondary findings module."""

from genetic_health.acmg import flag_acmg_findings, ACMG_GENES


def _make_finding(gene, category="pathogenic", gold_stars=3):
    """Helper to create a ClinVar-like finding."""
    return {
        "gene": gene,
        "rsid": "rs12345",
        "chromosome": "1",
        "position": "100",
        "user_genotype": "AG",
        "ref": "A",
        "alt": "G",
        "traits": "Test condition",
        "gold_stars": gold_stars,
        "is_homozygous": False,
        "is_heterozygous": True,
        "inheritance": "dominant",
    }


class TestACMGGenes:
    def test_gene_count(self):
        assert len(ACMG_GENES) == 81

    def test_key_genes_present(self):
        for gene in ["BRCA1", "BRCA2", "MLH1", "LDLR", "HFE", "SCN5A", "FBN1"]:
            assert gene in ACMG_GENES


class TestFlagACMGFindings:
    def test_known_acmg_gene_flagged(self):
        findings = {
            "pathogenic": [_make_finding("BRCA1")],
            "likely_pathogenic": [],
            "risk_factor": [],
            "drug_response": [],
            "protective": [],
        }
        result = flag_acmg_findings(findings)
        assert len(result["acmg_findings"]) == 1
        assert result["genes_with_variants"] == 1
        assert "breast" in result["acmg_findings"][0]["acmg_actionability"].lower()

    def test_non_acmg_gene_skipped(self):
        findings = {
            "pathogenic": [_make_finding("MTHFR")],
            "likely_pathogenic": [],
            "risk_factor": [],
            "drug_response": [],
            "protective": [],
        }
        result = flag_acmg_findings(findings)
        assert len(result["acmg_findings"]) == 0
        assert result["genes_with_variants"] == 0

    def test_empty_input(self):
        result = flag_acmg_findings(None)
        assert len(result["acmg_findings"]) == 0
        assert result["genes_screened"] == 81

    def test_empty_findings(self):
        findings = {
            "pathogenic": [],
            "likely_pathogenic": [],
            "risk_factor": [],
            "drug_response": [],
            "protective": [],
        }
        result = flag_acmg_findings(findings)
        assert len(result["acmg_findings"]) == 0
        assert "No pathogenic" in result["summary"]

    def test_likely_pathogenic_included(self):
        findings = {
            "pathogenic": [],
            "likely_pathogenic": [_make_finding("MLH1")],
            "risk_factor": [],
            "drug_response": [],
            "protective": [],
        }
        result = flag_acmg_findings(findings)
        assert len(result["acmg_findings"]) == 1
        assert result["acmg_findings"][0]["acmg_category"] == "likely_pathogenic"

    def test_risk_factor_not_included(self):
        """Risk factors should NOT be flagged as ACMG findings."""
        findings = {
            "pathogenic": [],
            "likely_pathogenic": [],
            "risk_factor": [_make_finding("LDLR")],
            "drug_response": [],
            "protective": [],
        }
        result = flag_acmg_findings(findings)
        assert len(result["acmg_findings"]) == 0

    def test_multiple_genes(self):
        findings = {
            "pathogenic": [_make_finding("BRCA1"), _make_finding("BRCA2"), _make_finding("HFE")],
            "likely_pathogenic": [_make_finding("MLH1")],
            "risk_factor": [],
            "drug_response": [],
            "protective": [],
        }
        result = flag_acmg_findings(findings)
        assert len(result["acmg_findings"]) == 4
        assert result["genes_with_variants"] == 4
        assert "Genetic counseling" in result["summary"]

    def test_result_keys(self):
        result = flag_acmg_findings(None)
        assert "acmg_findings" in result
        assert "genes_screened" in result
        assert "genes_with_variants" in result
        assert "summary" in result

    def test_sorted_by_gold_stars(self):
        findings = {
            "pathogenic": [
                _make_finding("BRCA1", gold_stars=1),
                _make_finding("BRCA2", gold_stars=4),
            ],
            "likely_pathogenic": [],
            "risk_factor": [],
            "drug_response": [],
            "protective": [],
        }
        result = flag_acmg_findings(findings)
        assert result["acmg_findings"][0]["gene"] == "BRCA2"
        assert result["acmg_findings"][1]["gene"] == "BRCA1"
