"""Tests for enhanced HTML report builders."""

from genetic_health.reports.enhanced_html import (
    build_ancestry_section,
    build_prs_section,
    build_epistasis_section,
    svg_ancestry_donut,
    svg_prs_gauge,
)


class TestAncestrySection:
    MOCK_ANCESTRY = {
        "proportions": {"EUR": 0.82, "AFR": 0.05, "EAS": 0.03, "SAS": 0.07, "AMR": 0.03},
        "markers_found": 45,
        "confidence": "high",
        "top_ancestry": "European",
        "details": [],
    }

    def test_builds_html_with_table(self):
        html = build_ancestry_section(self.MOCK_ANCESTRY)
        assert "European" in html
        assert "82.0%" in html
        assert "<table>" in html

    def test_empty_ancestry(self):
        html = build_ancestry_section({})
        assert "No ancestry" in html

    def test_none_ancestry(self):
        html = build_ancestry_section(None)
        assert "No ancestry" in html

    def test_zero_markers(self):
        html = build_ancestry_section({"markers_found": 0, "proportions": {}})
        assert "No ancestry" in html

    def test_confidence_shown(self):
        html = build_ancestry_section(self.MOCK_ANCESTRY)
        assert "High" in html
        assert "45 markers" in html


class TestPRSSection:
    MOCK_PRS = {
        "type2_diabetes": {
            "name": "Type 2 Diabetes",
            "raw_score": 1.5,
            "z_score": 0.8,
            "percentile": 78.8,
            "risk_category": "average",
            "snps_found": 20,
            "snps_total": 25,
            "ancestry_applicable": True,
            "ancestry_warning": "",
            "contributing_snps": [
                {"rsid": "rs7903146", "gene": "TCF7L2", "risk_allele": "T",
                 "copies": 2, "log_or": 0.322, "contribution": 0.644},
            ],
            "reference": "Mahajan 2018",
        },
        "hypertension": {
            "name": "Hypertension",
            "raw_score": 2.1,
            "z_score": 1.9,
            "percentile": 92.0,
            "risk_category": "elevated",
            "snps_found": 18,
            "snps_total": 22,
            "ancestry_applicable": True,
            "ancestry_warning": "",
            "contributing_snps": [],
            "reference": "Evangelou 2018",
        },
    }

    def test_builds_html_with_table(self):
        html = build_prs_section(self.MOCK_PRS)
        assert "Type 2 Diabetes" in html
        assert "Hypertension" in html
        assert "<table" in html

    def test_empty_prs(self):
        html = build_prs_section({})
        assert "No PRS" in html

    def test_none_prs(self):
        html = build_prs_section(None)
        assert "No PRS" in html

    def test_elevated_detail_shown(self):
        html = build_prs_section(self.MOCK_PRS)
        assert "Elevated Risk" in html
        assert "Hypertension" in html
        assert "92" in html

    def test_ancestry_warning_shown(self):
        prs = {
            "test": {
                "name": "Test",
                "raw_score": 0,
                "z_score": 0,
                "percentile": 50.0,
                "risk_category": "average",
                "snps_found": 5,
                "snps_total": 25,
                "ancestry_applicable": False,
                "ancestry_warning": "non-European warning",
                "contributing_snps": [],
                "reference": "Test 2024",
            },
        }
        html = build_prs_section(prs)
        assert "non-European" in html


class TestSVGGenerators:
    def test_ancestry_donut_returns_svg(self):
        data = {
            "proportions": {"EUR": 0.8, "AFR": 0.1, "EAS": 0.05, "SAS": 0.03, "AMR": 0.02},
            "top_ancestry": "European",
            "confidence": "high",
        }
        svg = svg_ancestry_donut(data)
        assert "<svg" in svg
        assert "European" in svg

    def test_ancestry_donut_empty(self):
        assert svg_ancestry_donut(None) == ""
        assert svg_ancestry_donut({}) == ""

    def test_prs_gauge_returns_svg(self):
        svg = svg_prs_gauge("T2D", 78.5, "average")
        assert "<svg" in svg
        assert "T2D" in svg
        assert "78" in svg

    def test_prs_gauge_high(self):
        svg = svg_prs_gauge("CAD", 96.0, "high")
        assert "High" in svg

    def test_prs_gauge_low(self):
        svg = svg_prs_gauge("HTN", 10.0, "low")
        assert "Low" in svg


class TestEpistasisSection:
    MOCK_EPISTASIS = [
        {
            "id": "mthfr_comt_methylation",
            "name": "MTHFR + COMT: Methylation-Catecholamine Interaction",
            "genes_involved": {"MTHFR": ["reduced"], "COMT": ["slow"]},
            "effect": "Dual burden on methylation and catecholamine clearance.",
            "risk_level": "high",
            "mechanism": "MTHFR reduces 5-MTHF, COMT requires SAMe.",
            "actions": ["Start methylfolate at LOW dose", "Prioritize magnesium"],
        },
    ]

    def test_builds_html_with_details(self):
        html = build_epistasis_section(self.MOCK_EPISTASIS)
        assert "MTHFR" in html
        assert "COMT" in html
        assert "HIGH" in html
        assert "<details" in html

    def test_actions_shown(self):
        html = build_epistasis_section(self.MOCK_EPISTASIS)
        assert "methylfolate" in html
        assert "magnesium" in html

    def test_empty_epistasis(self):
        html = build_epistasis_section([])
        assert "No significant" in html

    def test_none_epistasis(self):
        html = build_epistasis_section(None)
        assert "No significant" in html
