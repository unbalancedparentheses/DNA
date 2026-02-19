"""Tests for personalized recommendations synthesis."""

from genetic_health.recommendations import (
    generate_recommendations,
    RISK_GROUPS,
    SPECIALIST_REFERRALS,
    _compute_priority,
    _deduplicate_monitoring,
    _build_specialist_referrals,
    _build_clinical_insights,
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
        assert "specialist_referrals" in recs
        assert "clinical_insights" in recs
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


class TestFalseAlarmFixes:
    """Verify that common polymorphisms and low-evidence ClinVar entries
    no longer trigger inappropriate high-priority alarms."""

    def test_tp53_arg72_does_not_trigger_cancer_screening(self):
        """TP53 rs1042522 Arg72Pro (magnitude 1) should NOT trigger cancer_screening."""
        findings = [
            {"gene": "TP53", "status": "arg72", "magnitude": 1,
             "rsid": "rs1042522", "genotype": "GG", "category": "Longevity",
             "description": "Arg72Pro common polymorphism", "note": ""},
        ]
        recs = generate_recommendations(findings)
        cancer = [p for p in recs["priorities"] if p["id"] == "cancer_screening"]
        assert len(cancer) == 0, "TP53 Arg72 should not trigger cancer screening"

    def test_tp53_not_in_cancer_screening_genes(self):
        """TP53 should have been removed from cancer_screening gene list."""
        assert "TP53" not in RISK_GROUPS["cancer_screening"]["genes"]

    def test_magnitude_1_gene_not_counted_as_signal(self):
        """Magnitude 1 findings should not contribute gene_signal."""
        findings = [
            {"gene": "MC1R", "status": "fair_skin", "magnitude": 1,
             "rsid": "rs1805007", "genotype": "CT", "category": "Skin",
             "description": "Fair skin", "note": ""},
        ]
        recs = generate_recommendations(findings)
        skin = [p for p in recs["priorities"] if p["id"] == "skin"]
        assert len(skin) == 0, "Magnitude 1 should not trigger skin group"

    def test_magnitude_2_gene_counts_as_signal(self):
        """Magnitude 2 findings should still contribute gene_signal."""
        findings = [
            {"gene": "MC1R", "status": "fair_skin", "magnitude": 2,
             "rsid": "rs1805007", "genotype": "CT", "category": "Skin",
             "description": "Fair skin", "note": ""},
        ]
        recs = generate_recommendations(findings)
        skin = [p for p in recs["priorities"] if p["id"] == "skin"]
        assert len(skin) == 1

    def test_low_evidence_clinvar_not_triggering_alarm(self):
        """ClinVar entries with gold_stars < 2 should not trigger keyword matching."""
        disease = {
            "pathogenic": [
                {"gene": "ANKRD36", "traits": "malignant tumor of esophagus",
                 "gold_stars": 0, "rsid": "rs9999", "user_genotype": "AG",
                 "chromosome": "2", "position": "100", "ref_allele": "A",
                 "alt_allele": "G", "is_homozygous": False,
                 "is_heterozygous": True, "inheritance": ""},
            ],
            "likely_pathogenic": [],
            "risk_factor": [],
            "drug_response": [],
            "protective": [],
        }
        recs = generate_recommendations([], disease_findings=disease)
        cancer = [p for p in recs["priorities"] if p["id"] == "cancer_screening"]
        assert len(cancer) == 0, (
            "Low-evidence ClinVar entry should not trigger cancer screening"
        )

    def test_high_evidence_clinvar_still_triggers(self):
        """ClinVar entries with gold_stars >= 2 should still trigger keyword matching."""
        disease = {
            "pathogenic": [
                {"gene": "BRCA1", "traits": "breast-ovarian cancer",
                 "gold_stars": 3, "rsid": "rs80357906", "user_genotype": "AG",
                 "chromosome": "17", "position": "100", "ref_allele": "A",
                 "alt_allele": "G", "is_homozygous": False,
                 "is_heterozygous": True, "inheritance": "dominant"},
            ],
            "likely_pathogenic": [],
            "risk_factor": [],
            "drug_response": [],
            "protective": [],
        }
        recs = generate_recommendations([], disease_findings=disease)
        cancer = [p for p in recs["priorities"] if p["id"] == "cancer_screening"]
        assert len(cancer) == 1
        assert cancer[0]["priority"] == "high"

    def test_cancer_keywords_tightened(self):
        """The word 'tumor' alone should not be in cancer_screening keywords."""
        assert "tumor" not in RISK_GROUPS["cancer_screening"]["disease_keywords"]
        assert "malignant neoplasm" in RISK_GROUPS["cancer_screening"]["disease_keywords"]

    def test_serpina1_does_not_trigger_cancer_via_secondary_traits(self):
        """SERPINA1 with cancer in secondary ClinVar traits should NOT
        trigger cancer_screening — only the primary trait is checked."""
        disease = {
            "pathogenic": [
                {"gene": "SERPINA1",
                 "traits": ("PI S|Alpha-1-antitrypsin deficiency|"
                            "not provided|"
                            "Alpha-1-antitrypsin deficiency;COPD|"
                            "Cystic fibrosis|not specified|"
                            "Inborn genetic diseases|COVID-19|"
                            "SERPINA1-related disorder|"
                            "Hepatocellular carcinoma|"
                            "Colorectal cancer"),
                 "gold_stars": 3, "rsid": "rs28929474", "user_genotype": "AT",
                 "chromosome": "14", "position": "94844947",
                 "ref_allele": "A", "alt_allele": "T",
                 "is_homozygous": True, "is_heterozygous": False,
                 "inheritance": "recessive"},
            ],
            "likely_pathogenic": [],
            "risk_factor": [],
            "drug_response": [],
            "protective": [],
        }
        recs = generate_recommendations([], disease_findings=disease)
        cancer = [p for p in recs["priorities"] if p["id"] == "cancer_screening"]
        assert len(cancer) == 0, (
            "SERPINA1 should not trigger cancer_screening via secondary traits"
        )
        # But it SHOULD trigger lung_health
        lung = [p for p in recs["priorities"] if p["id"] == "lung_health"]
        assert len(lung) == 1, "SERPINA1 should trigger lung_health"


class TestSpecialistReferrals:
    def test_return_dict_has_specialist_referrals_key(self):
        recs = generate_recommendations([])
        assert "specialist_referrals" in recs
        assert isinstance(recs["specialist_referrals"], list)

    def test_all_referral_defs_have_required_fields(self):
        for ref_id, ref in SPECIALIST_REFERRALS.items():
            assert "title" in ref, f"{ref_id} missing title"
            assert "triggers" in ref, f"{ref_id} missing triggers"
            assert "reason_template" in ref, f"{ref_id} missing reason_template"

    def test_acmg_triggers_genetic_counselor(self):
        acmg = {
            "acmg_findings": [
                {"gene": "SERPINA1", "traits": "Alpha-1 antitrypsin deficiency",
                 "gold_stars": 3, "rsid": "rs28929474",
                 "acmg_actionability": "Referral"},
            ],
            "genes_screened": 81,
            "genes_with_variants": 1,
            "summary": "1 variant found",
        }
        recs = generate_recommendations([], acmg=acmg)
        specialists = [r["specialist"] for r in recs["specialist_referrals"]]
        assert "Genetic Counselor" in specialists

    def test_acmg_serpina1_triggers_pulmonologist(self):
        acmg = {
            "acmg_findings": [
                {"gene": "SERPINA1", "traits": "Alpha-1 antitrypsin deficiency",
                 "gold_stars": 3, "rsid": "rs28929474",
                 "acmg_actionability": "Referral"},
            ],
            "genes_screened": 81,
            "genes_with_variants": 1,
            "summary": "1 variant found",
        }
        recs = generate_recommendations([], acmg=acmg)
        specialists = [r["specialist"] for r in recs["specialist_referrals"]]
        assert "Pulmonologist" in specialists

    def test_actionable_star_alleles_trigger_pharmacist(self):
        star = {
            "CYP2C9": {
                "gene": "CYP2C9", "diplotype": "*1/*2",
                "phenotype": "intermediate",
                "snps_found": 3, "snps_total": 3,
                "clinical_note": "Lower warfarin dose"},
        }
        recs = generate_recommendations([], star_alleles=star)
        specialists = [r["specialist"] for r in recs["specialist_referrals"]]
        assert "Pharmacist (PGx Review)" in specialists

    def test_hfe_carrier_triggers_hematologist(self):
        findings = [
            {"gene": "HFE", "status": "carrier", "magnitude": 2,
             "rsid": "rs1799945", "genotype": "CG", "category": "Iron",
             "description": "H63D carrier", "note": ""},
        ]
        recs = generate_recommendations(findings)
        specialists = [r["specialist"] for r in recs["specialist_referrals"]]
        assert "Hematologist" in specialists

    def test_skin_priority_triggers_dermatologist(self):
        findings = [
            {"gene": "MC1R", "status": "fair_skin", "magnitude": 2,
             "rsid": "rs1805007", "genotype": "CT", "category": "Skin",
             "description": "Fair skin", "note": ""},
        ]
        recs = generate_recommendations(findings)
        specialists = [r["specialist"] for r in recs["specialist_referrals"]]
        assert "Dermatologist" in specialists

    def test_no_referrals_for_empty_input(self):
        recs = generate_recommendations([])
        assert recs["specialist_referrals"] == []


class TestClinicalInsights:
    def test_return_dict_has_clinical_insights_key(self):
        recs = generate_recommendations([])
        assert "clinical_insights" in recs
        assert isinstance(recs["clinical_insights"], list)

    def test_clinical_context_extracted(self):
        findings = [
            {"gene": "COMT", "status": "slow", "magnitude": 3,
             "rsid": "rs4680", "genotype": "AA", "category": "Neurotransmitters",
             "description": "Slow COMT", "note": ""},
        ]
        recs = generate_recommendations(findings)
        genes = [ci["gene"] for ci in recs["clinical_insights"]]
        assert "COMT" in genes

    def test_clinical_insight_has_mechanism(self):
        findings = [
            {"gene": "COMT", "status": "slow", "magnitude": 3,
             "rsid": "rs4680", "genotype": "AA", "category": "Neurotransmitters",
             "description": "Slow COMT", "note": ""},
        ]
        recs = generate_recommendations(findings)
        comt = next(ci for ci in recs["clinical_insights"] if ci["gene"] == "COMT")
        assert comt["mechanism"], "Should have mechanism from clinical_context"
        assert comt["actions"], "Should have actions from clinical_context"

    def test_low_magnitude_excluded(self):
        findings = [
            {"gene": "COMT", "status": "slow", "magnitude": 0,
             "rsid": "rs4680", "genotype": "AA", "category": "Neurotransmitters",
             "description": "Slow COMT", "note": ""},
        ]
        recs = generate_recommendations(findings)
        genes = [ci["gene"] for ci in recs["clinical_insights"]]
        assert "COMT" not in genes

    def test_clinical_actions_overlaid_on_priorities(self):
        """Clinical context actions should be appended to matching priority groups."""
        findings = [
            {"gene": "AGTR1", "status": "increased", "magnitude": 3,
             "rsid": "rs5186", "genotype": "CC", "category": "Cardiovascular",
             "description": "BP risk", "note": ""},
        ]
        recs = generate_recommendations(findings)
        bp = [p for p in recs["priorities"] if p["id"] == "blood_pressure"]
        assert len(bp) == 1
        # Should have clinical_actions key with AGTR1-specific actions
        assert "clinical_actions" in bp[0]
        assert len(bp[0]["clinical_actions"]) > 0
