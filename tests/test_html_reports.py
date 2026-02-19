"""Tests for enhanced HTML report builders."""

from genetic_health.reports.enhanced_html import (
    build_ancestry_section,
    build_prs_section,
    build_epistasis_section,
    build_disease_risk,
    build_monitoring,
    build_doctor_card,
    build_nutrition_section,
    build_protective,
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


class TestDiseaseRiskSection:
    MOCK_DISEASE = {
        "pathogenic": [
            {"gene": "BRCA1", "traits": "Breast cancer;Ovarian cancer",
             "user_genotype": "AG", "gold_stars": 3, "zygosity": "heterozygous"},
        ],
        "likely_pathogenic": [
            {"gene": "MLH1", "traits": "Lynch syndrome",
             "user_genotype": "CT", "gold_stars": 2, "zygosity": "heterozygous"},
        ],
        "risk_factor": [
            {"gene": "APOE", "traits": "Alzheimer disease",
             "user_genotype": "CT", "gold_stars": 4, "zygosity": "heterozygous"},
        ],
        "drug_response": [
            {"gene": "CYP2C19", "traits": "Clopidogrel response",
             "user_genotype": "GA", "gold_stars": 3, "zygosity": "heterozygous"},
        ],
        "protective": [
            {"gene": "PCSK9", "traits": "Lower LDL cholesterol"},
        ],
        "carriers": [],
    }

    def test_builds_tables(self):
        html = build_disease_risk(self.MOCK_DISEASE)
        assert "BRCA1" in html
        assert "MLH1" in html
        assert "APOE" in html
        assert "CYP2C19" in html
        assert "<table>" in html

    def test_pathogenic_shown(self):
        html = build_disease_risk(self.MOCK_DISEASE)
        assert "Pathogenic Variants" in html
        assert "Breast cancer" in html

    def test_protective_cards(self):
        html = build_disease_risk(self.MOCK_DISEASE)
        assert "Protective Variants" in html
        assert "PCSK9" in html
        assert "good-news-card" in html

    def test_empty_disease(self):
        html = build_disease_risk({})
        assert "No ClinVar" in html

    def test_none_disease(self):
        html = build_disease_risk(None)
        assert "No ClinVar" in html

    def test_partial_data(self):
        data = {"pathogenic": [], "likely_pathogenic": [], "risk_factor": [],
                "drug_response": [], "protective": [], "carriers": []}
        html = build_disease_risk(data)
        assert "ClinVar" in html


class TestMonitoringSection:
    MOCK_RECS = {
        "monitoring_schedule": [
            {"test": "Blood pressure", "frequency": "Weekly (home)", "reason": "Multiple BP genes"},
            {"test": "Lipid panel", "frequency": "Annually", "reason": "APOE e4 carrier"},
            {"test": "Iron panel", "frequency": "Baseline", "reason": "HFE carrier"},
        ],
    }

    def test_builds_table(self):
        html = build_monitoring(self.MOCK_RECS)
        assert "Blood pressure" in html
        assert "Lipid panel" in html
        assert "<table>" in html

    def test_frequency_badges(self):
        html = build_monitoring(self.MOCK_RECS)
        assert "Weekly" in html
        assert "Annually" in html
        assert "mag-badge" in html

    def test_print_callout(self):
        html = build_monitoring(self.MOCK_RECS)
        assert "Print this" in html or "doctor visit" in html

    def test_empty_schedule(self):
        html = build_monitoring({"monitoring_schedule": []})
        assert "No monitoring" in html

    def test_none_data(self):
        html = build_monitoring(None)
        assert "No monitoring" in html


class TestDoctorCardSection:
    MOCK_RECS = {
        "priorities": [
            {"id": "blood_pressure", "title": "Blood Pressure", "priority": "high",
             "why": "AGTR1+AGT", "actions": ["Monitor BP"], "doctor_note": "Check BP meds",
             "monitoring": [], "signal_count": 3},
        ],
        "specialist_referrals": [
            {"specialist": "Cardiologist", "reason": "Multiple CV risk genes", "urgency": "soon"},
        ],
        "monitoring_schedule": [
            {"test": "Blood pressure", "frequency": "Weekly", "reason": "BP genes"},
        ],
    }
    MOCK_STAR = {
        "CYP2C19": {"diplotype": "*1/*2", "phenotype": "intermediate",
                     "snps_found": 4, "snps_total": 5, "clinical_note": "Reduced clopidogrel"},
    }
    MOCK_APOE = {"apoe_type": "e3/e4", "risk_level": "elevated",
                 "alzheimer_or": 2.8, "confidence": "high", "description": "One e4 allele"}
    MOCK_ACMG = {
        "genes_screened": 81, "genes_with_variants": 1,
        "acmg_findings": [
            {"gene": "BRCA2", "traits": "Breast cancer", "user_genotype": "AG",
             "gold_stars": 3, "acmg_actionability": "Enhanced screening"},
        ],
    }

    def test_builds_doctor_card(self):
        html = build_doctor_card(self.MOCK_RECS, self.MOCK_STAR, self.MOCK_APOE, self.MOCK_ACMG)
        assert "Patient Genetic Summary" in html
        assert "Blood Pressure" in html

    def test_pharmacogenomic_table(self):
        html = build_doctor_card(self.MOCK_RECS, self.MOCK_STAR, self.MOCK_APOE, self.MOCK_ACMG)
        assert "CYP2C19" in html
        assert "*1/*2" in html
        assert "Intermediate" in html

    def test_apoe_shown(self):
        html = build_doctor_card(self.MOCK_RECS, self.MOCK_STAR, self.MOCK_APOE, self.MOCK_ACMG)
        assert "e3/e4" in html
        assert "Elevated" in html

    def test_acmg_findings(self):
        html = build_doctor_card(self.MOCK_RECS, self.MOCK_STAR, self.MOCK_APOE, self.MOCK_ACMG)
        assert "BRCA2" in html
        assert "genetic counselor" in html.lower()

    def test_specialist_referrals(self):
        html = build_doctor_card(self.MOCK_RECS, self.MOCK_STAR, self.MOCK_APOE, self.MOCK_ACMG)
        assert "Cardiologist" in html
        assert "soon" in html

    def test_empty_data(self):
        html = build_doctor_card({}, {}, {}, {})
        assert "No significant" in html

    def test_none_data(self):
        html = build_doctor_card(None, None, None, None)
        assert "No significant" in html


class TestNutritionSection:
    MOCK_RECS = {
        "priorities": [
            {"id": "methylation", "title": "Methylation Support", "priority": "high",
             "why": "MTHFR+COMT", "actions": ["Take methylfolate 400mcg"],
             "clinical_actions": ["Monitor homocysteine"], "doctor_note": "",
             "monitoring": [], "signal_count": 2},
            {"id": "blood_pressure", "title": "Blood Pressure", "priority": "high",
             "why": "AGTR1", "actions": ["Reduce sodium"], "doctor_note": "",
             "monitoring": [], "signal_count": 3},
        ],
    }
    MOCK_INSIGHTS = {
        "narratives": [
            {"id": "methylation_choline", "title": "Methylation & Choline",
             "matched_genes": ["MTHFR", "PEMT"],
             "narrative": "Your methylation cycle genes suggest...",
             "references": ["PMID:12345"],
             "practical": "Eat eggs and leafy greens"},
        ],
    }

    def test_builds_nutrition_content(self):
        html = build_nutrition_section(self.MOCK_RECS, self.MOCK_INSIGHTS)
        assert "Methylation" in html
        assert "methylfolate" in html

    def test_filters_nutrition_priorities(self):
        html = build_nutrition_section(self.MOCK_RECS, self.MOCK_INSIGHTS)
        # Blood pressure is not a nutrition group, so should be excluded
        assert "Blood Pressure" not in html

    def test_shows_narratives(self):
        html = build_nutrition_section(self.MOCK_RECS, self.MOCK_INSIGHTS)
        assert "Methylation &amp; Choline" in html or "Methylation & Choline" in html
        assert "eggs" in html

    def test_clinical_actions(self):
        html = build_nutrition_section(self.MOCK_RECS, self.MOCK_INSIGHTS)
        assert "homocysteine" in html

    def test_empty_data(self):
        html = build_nutrition_section({}, {})
        assert "No nutrition" in html

    def test_none_data(self):
        html = build_nutrition_section(None, None)
        assert "No nutrition" in html


class TestProtectiveSection:
    MOCK_RECS = {
        "good_news": [
            {"gene": "FOXO3", "description": "Longevity-associated variant"},
            {"gene": "APOE", "description": "No e4 alleles — average Alzheimer's risk"},
        ],
    }
    MOCK_INSIGHTS = {
        "protective_findings": [
            {"gene": "CHRNA5", "status": "favorable",
             "title": "Reduced nicotine dependence",
             "finding": "Your variant is associated with lower addiction risk.",
             "reference": "Thorgeirsson 2010, PMID: 20418890"},
        ],
    }

    def test_builds_good_news_cards(self):
        html = build_protective(self.MOCK_RECS, self.MOCK_INSIGHTS)
        assert "FOXO3" in html
        assert "Longevity" in html
        assert "good-news-card" in html

    def test_builds_research_backed(self):
        html = build_protective(self.MOCK_RECS, self.MOCK_INSIGHTS)
        assert "CHRNA5" in html
        assert "nicotine" in html
        assert "Thorgeirsson" in html

    def test_empty_data(self):
        html = build_protective({}, {})
        assert "No protective" in html

    def test_none_data(self):
        html = build_protective(None, None)
        assert "No protective" in html

    def test_good_news_only(self):
        html = build_protective(self.MOCK_RECS, {})
        assert "FOXO3" in html
        assert "good-news-card" in html

    def test_protective_findings_only(self):
        html = build_protective({}, self.MOCK_INSIGHTS)
        assert "CHRNA5" in html
