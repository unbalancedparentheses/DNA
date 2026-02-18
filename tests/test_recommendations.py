"""Tests for personalized recommendations synthesis."""

from genetic_health.recommendations import (
    generate_recommendations,
    RISK_GROUPS,
    _compute_priority,
    _deduplicate_monitoring,
)


class TestRiskGroupDefinitions:
    def test_all_groups_have_required_fields(self):
        for gid, group in RISK_GROUPS.items():
            assert "title" in group, f"{gid} missing title"
            assert "genes" in group, f"{gid} missing genes"
            assert "actions" in group, f"{gid} missing actions"
            assert "doctor_note" in group, f"{gid} missing doctor_note"
            assert "monitoring" in group, f"{gid} missing monitoring"
            assert isinstance(group["genes"], set)
            assert isinstance(group["actions"], list)

    def test_monitoring_entries_have_required_fields(self):
        for gid, group in RISK_GROUPS.items():
            for m in group["monitoring"]:
                assert "test" in m, f"{gid} monitoring missing test"
                assert "frequency" in m, f"{gid} monitoring missing frequency"
                assert "reason" in m, f"{gid} monitoring missing reason"


class TestPriorityScoring:
    def test_pathogenic_is_high(self):
        assert _compute_priority(["pathogenic"]) == "high"

    def test_prs_high_is_high(self):
        assert _compute_priority(["prs_high"]) == "high"

    def test_three_gene_signals_is_high(self):
        assert _compute_priority(["gene_signal"] * 3) == "high"

    def test_prs_elevated_is_moderate(self):
        assert _compute_priority(["prs_elevated"]) == "moderate"

    def test_two_gene_signals_is_moderate(self):
        assert _compute_priority(["gene_signal"] * 2) == "moderate"

    def test_high_magnitude_is_moderate(self):
        assert _compute_priority(["gene_signal", "high_magnitude"]) == "moderate"

    def test_epistasis_is_moderate(self):
        assert _compute_priority(["epistasis"]) == "moderate"

    def test_single_gene_signal_is_low(self):
        assert _compute_priority(["gene_signal"]) == "low"


class TestConvergentBPSignals:
    def test_multiple_bp_genes_plus_prs_gives_high(self):
        findings = [
            {"gene": "AGTR1", "status": "elevated", "magnitude": 3,
             "rsid": "rs5186", "genotype": "CC", "category": "Cardiovascular",
             "description": "Elevated BP risk", "note": ""},
            {"gene": "AGT", "status": "elevated", "magnitude": 2,
             "rsid": "rs699", "genotype": "TT", "category": "Cardiovascular",
             "description": "Elevated angiotensinogen", "note": ""},
            {"gene": "ACE", "status": "elevated", "magnitude": 2,
             "rsid": "rs4340", "genotype": "DD", "category": "Cardiovascular",
             "description": "Elevated ACE", "note": ""},
        ]
        prs_results = {
            "hypertension": {
                "name": "Hypertension",
                "percentile": 97,
                "risk_category": "high",
                "snps_found": 10,
                "snps_total": 15,
                "ancestry_applicable": True,
                "reference": "Evangelou 2018",
                "contributing_snps": [],
            }
        }
        recs = generate_recommendations(findings, prs_results=prs_results)
        bp_recs = [p for p in recs["priorities"] if p["id"] == "blood_pressure"]
        assert len(bp_recs) == 1
        assert bp_recs[0]["priority"] == "high"

    def test_single_bp_gene_is_low(self):
        findings = [
            {"gene": "AGTR1", "status": "elevated", "magnitude": 2,
             "rsid": "rs5186", "genotype": "CC", "category": "Cardiovascular",
             "description": "Elevated BP risk", "note": ""},
        ]
        recs = generate_recommendations(findings)
        bp_recs = [p for p in recs["priorities"] if p["id"] == "blood_pressure"]
        assert len(bp_recs) == 1
        assert bp_recs[0]["priority"] == "low"


class TestDrugCard:
    def test_cyp_findings_in_drug_card(self):
        findings = [
            {"gene": "CYP2C9", "status": "poor", "magnitude": 3,
             "rsid": "rs1799853", "genotype": "CT", "category": "Drug Metabolism",
             "description": "Poor metabolizer", "note": ""},
            {"gene": "VKORC1", "status": "sensitive", "magnitude": 3,
             "rsid": "rs9923231", "genotype": "CT", "category": "Drug Metabolism",
             "description": "Warfarin sensitive", "note": ""},
        ]
        recs = generate_recommendations(findings)
        genes_in_card = {entry["gene"] for entry in recs["drug_card"]}
        assert "CYP2C9" in genes_in_card
        assert "VKORC1" in genes_in_card

    def test_clinvar_drug_response_in_card(self):
        disease = {
            "pathogenic": [], "likely_pathogenic": [],
            "risk_factor": [], "protective": [],
            "drug_response": [
                {"gene": "CYP2D6", "rsid": "rs1234", "user_genotype": "AG",
                 "traits": "Codeine metabolism"},
            ],
        }
        recs = generate_recommendations([], disease_findings=disease)
        genes_in_card = {entry["gene"] for entry in recs["drug_card"]}
        assert "CYP2D6" in genes_in_card


class TestMonitoringSchedule:
    def test_deduplication(self):
        priorities = [
            {"monitoring": [
                {"test": "Blood pressure", "frequency": "Weekly (home)", "reason": "A"},
                {"test": "Lipid panel", "frequency": "Annually", "reason": "B"},
            ]},
            {"monitoring": [
                {"test": "Blood pressure", "frequency": "Annually", "reason": "C"},
                {"test": "Ferritin", "frequency": "Annually", "reason": "D"},
            ]},
        ]
        schedule = _deduplicate_monitoring(priorities)
        tests = [m["test"] for m in schedule]
        assert tests.count("Blood pressure") == 1
        # Keep the more frequent one (Weekly)
        bp = next(m for m in schedule if m["test"] == "Blood pressure")
        assert bp["frequency"] == "Weekly (home)"

    def test_sorted_by_urgency(self):
        priorities = [
            {"monitoring": [
                {"test": "Annually test", "frequency": "Annually", "reason": "A"},
                {"test": "Weekly test", "frequency": "Weekly (home)", "reason": "B"},
            ]},
        ]
        schedule = _deduplicate_monitoring(priorities)
        assert schedule[0]["test"] == "Weekly test"


class TestEmptyInputs:
    def test_none_findings(self):
        recs = generate_recommendations(None)
        assert "priorities" in recs
        assert "drug_card" in recs
        assert "monitoring_schedule" in recs
        assert "good_news" in recs
        assert isinstance(recs["priorities"], list)

    def test_empty_findings(self):
        recs = generate_recommendations([])
        assert recs["priorities"] == []
        assert recs["drug_card"] == []

    def test_none_disease_findings(self):
        recs = generate_recommendations([], disease_findings=None)
        assert isinstance(recs["priorities"], list)

    def test_none_prs(self):
        recs = generate_recommendations([], prs_results=None)
        assert isinstance(recs["priorities"], list)


class TestPriorityOrdering:
    def test_high_before_moderate_before_low(self):
        findings = [
            {"gene": "AGTR1", "status": "elevated", "magnitude": 3,
             "rsid": "rs5186", "genotype": "CC", "category": "Cardiovascular",
             "description": "BP risk", "note": ""},
            {"gene": "AGT", "status": "elevated", "magnitude": 2,
             "rsid": "rs699", "genotype": "TT", "category": "Cardiovascular",
             "description": "BP risk", "note": ""},
            {"gene": "ACE", "status": "elevated", "magnitude": 2,
             "rsid": "rs4340", "genotype": "DD", "category": "Cardiovascular",
             "description": "BP risk", "note": ""},
            {"gene": "MC1R", "status": "fair_skin", "magnitude": 2,
             "rsid": "rs1805007", "genotype": "CT", "category": "Skin",
             "description": "Fair skin", "note": ""},
        ]
        recs = generate_recommendations(findings)
        priority_order = {"high": 0, "moderate": 1, "low": 2}
        levels = [priority_order[p["priority"]] for p in recs["priorities"]]
        assert levels == sorted(levels)


class TestGoodNews:
    def test_protective_findings_collected(self):
        findings = [
            {"gene": "FOXO3", "status": "longevity", "magnitude": 1,
             "rsid": "rs2802292", "genotype": "GT", "category": "Longevity",
             "description": "Longevity variant", "note": ""},
        ]
        disease = {
            "pathogenic": [], "likely_pathogenic": [],
            "risk_factor": [], "drug_response": [],
            "protective": [
                {"gene": "APOE", "traits": "Reduced Alzheimer risk"},
            ],
        }
        recs = generate_recommendations(findings, disease_findings=disease)
        genes = [g["gene"] for g in recs["good_news"]]
        assert "FOXO3" in genes
        assert "APOE" in genes
