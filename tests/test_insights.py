"""Tests for research-backed genomic insights module."""

from genetic_health.insights import (
    generate_insights,
    SINGLE_GENE_INSIGHTS,
    MULTI_GENE_NARRATIVES,
)


class TestSingleGeneInsightDatabase:
    def test_all_entries_have_required_fields(self):
        for (gene, status), entries in SINGLE_GENE_INSIGHTS.items():
            assert isinstance(gene, str) and gene, f"Invalid gene: {gene}"
            assert isinstance(status, str) and status, f"Invalid status for {gene}"
            for entry in entries:
                assert "title" in entry, f"Missing title for ({gene}, {status})"
                assert "finding" in entry, f"Missing finding for ({gene}, {status})"
                assert "reference" in entry, f"Missing reference for ({gene}, {status})"
                assert "practical" in entry, f"Missing practical for ({gene}, {status})"

    def test_database_has_pharmacogenomic_entries(self):
        pgx_genes = {"CYP2C9", "DPYD", "CYP1A2"}
        found = {g for (g, _) in SINGLE_GENE_INSIGHTS.keys()}
        assert pgx_genes & found, "Should have pharmacogenomic gene entries"

    def test_database_has_longevity_entries(self):
        found = {g for (g, _) in SINGLE_GENE_INSIGHTS.keys()}
        assert "CETP" in found or "APOE" in found


class TestMultiGeneNarratives:
    def test_all_narratives_have_required_fields(self):
        for n in MULTI_GENE_NARRATIVES:
            assert "id" in n, f"Missing id"
            assert "title" in n
            assert "required_genes" in n
            assert "optional_genes" in n
            assert "min_matches" in n
            assert "narrative" in n
            assert "references" in n
            assert "practical" in n

    def test_bp_cluster_defined(self):
        ids = {n["id"] for n in MULTI_GENE_NARRATIVES}
        assert "bp_cluster" in ids

    def test_caffeine_profile_defined(self):
        ids = {n["id"] for n in MULTI_GENE_NARRATIVES}
        assert "caffeine_profile" in ids


class TestGenerateInsights:
    def test_empty_input(self):
        result = generate_insights([])
        assert "single_gene" in result
        assert "narratives" in result
        assert "genome_highlights" in result
        assert "protective_findings" in result
        assert isinstance(result["single_gene"], list)
        assert isinstance(result["narratives"], list)

    def test_none_input(self):
        result = generate_insights(None)
        assert isinstance(result["single_gene"], list)

    def test_single_gene_matching(self):
        findings = [
            {"gene": "COMT", "status": "slow", "magnitude": 3,
             "rsid": "rs4680", "genotype": "AA", "category": "Neurotransmitters",
             "description": "Slow COMT", "note": ""},
        ]
        result = generate_insights(findings)
        genes = [e["gene"] for e in result["single_gene"]]
        assert "COMT" in genes

    def test_no_match_for_unrecognized_status(self):
        findings = [
            {"gene": "COMT", "status": "unknown_xyz", "magnitude": 3,
             "rsid": "rs4680", "genotype": "AA", "category": "Neurotransmitters",
             "description": "Unknown", "note": ""},
        ]
        result = generate_insights(findings)
        comt = [e for e in result["single_gene"] if e["gene"] == "COMT"]
        assert len(comt) == 0

    def test_multi_gene_bp_narrative(self):
        findings = [
            {"gene": "AGTR1", "status": "increased", "magnitude": 3,
             "rsid": "rs5186", "genotype": "CC", "category": "Cardiovascular",
             "description": "BP risk", "note": ""},
            {"gene": "ACE", "status": "high", "magnitude": 2,
             "rsid": "rs4340", "genotype": "DD", "category": "Cardiovascular",
             "description": "High ACE", "note": ""},
        ]
        result = generate_insights(findings)
        narrative_ids = [n["id"] for n in result["narratives"]]
        assert "bp_cluster" in narrative_ids

    def test_caffeine_narrative_with_single_required(self):
        findings = [
            {"gene": "CYP1A2", "status": "slow", "magnitude": 2,
             "rsid": "rs762551", "genotype": "AA", "category": "Caffeine",
             "description": "Slow", "note": ""},
        ]
        result = generate_insights(findings)
        narrative_ids = [n["id"] for n in result["narratives"]]
        assert "caffeine_profile" in narrative_ids

    def test_narrative_requires_at_least_one_required_gene(self):
        # AGT is optional for bp_cluster but not required — should NOT trigger
        findings = [
            {"gene": "AGT", "status": "increased", "magnitude": 2,
             "rsid": "rs699", "genotype": "TT", "category": "Cardiovascular",
             "description": "High", "note": ""},
            {"gene": "GNB3", "status": "increased", "magnitude": 2,
             "rsid": "rs5443", "genotype": "TT", "category": "Cardiovascular",
             "description": "High", "note": ""},
        ]
        result = generate_insights(findings)
        narrative_ids = [n["id"] for n in result["narratives"]]
        assert "bp_cluster" not in narrative_ids

    def test_apoe_e2_in_highlights(self):
        apoe = {"apoe_type": "e2/e3", "risk_level": "reduced",
                "confidence": "high", "alzheimer_or": 0.6,
                "description": "Protective"}
        result = generate_insights([], apoe=apoe)
        highlight_titles = [h["title"] for h in result["genome_highlights"]]
        assert any("APOE" in t for t in highlight_titles)

    def test_star_alleles_in_insights(self):
        star = {
            "CYP2C9": {
                "gene": "CYP2C9", "diplotype": "*1/*2",
                "phenotype": "intermediate", "snps_found": 3,
                "snps_total": 3, "clinical_note": "Lower warfarin dose"},
        }
        result = generate_insights([], star_alleles=star)
        genes = [e["gene"] for e in result["single_gene"]]
        assert "CYP2C9" in genes

    def test_protective_findings(self):
        findings = [
            {"gene": "CETP", "status": "favorable", "magnitude": 2,
             "rsid": "rs5882", "genotype": "GG", "category": "Longevity",
             "description": "Favorable lipids", "note": ""},
        ]
        result = generate_insights(findings)
        genes = [p["gene"] for p in result["protective_findings"]]
        assert "CETP" in genes

    def test_genome_highlights_max_five(self):
        # Create many high-mag findings
        findings = [
            {"gene": f"GENE{i}", "status": "high", "magnitude": 4,
             "rsid": f"rs{i}", "genotype": "AA", "category": "Test",
             "description": f"Finding {i}", "note": ""}
            for i in range(10)
        ]
        result = generate_insights(findings)
        assert len(result["genome_highlights"]) <= 5
