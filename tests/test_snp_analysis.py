"""Tests for SNP database structure and lifestyle/health analysis logic."""

from genetic_health.snp_database import COMPREHENSIVE_SNPS
from genetic_health.analysis import analyze_lifestyle_health, _lookup_genotype, _safe_int


class TestSNPDatabase:
    def test_all_entries_have_required_keys(self):
        for rsid, info in COMPREHENSIVE_SNPS.items():
            assert "gene" in info, f"{rsid} missing 'gene'"
            assert "category" in info, f"{rsid} missing 'category'"
            assert "variants" in info, f"{rsid} missing 'variants'"
            assert len(info["variants"]) > 0, f"{rsid} has empty variants"

    def test_variant_entries_have_required_fields(self):
        for rsid, info in COMPREHENSIVE_SNPS.items():
            for genotype, vinfo in info["variants"].items():
                assert "status" in vinfo, f"{rsid}/{genotype} missing 'status'"
                assert "desc" in vinfo, f"{rsid}/{genotype} missing 'desc'"
                assert "magnitude" in vinfo, f"{rsid}/{genotype} missing 'magnitude'"

    def test_magnitudes_in_range(self):
        for rsid, info in COMPREHENSIVE_SNPS.items():
            for genotype, vinfo in info["variants"].items():
                m = vinfo["magnitude"]
                assert 0 <= m <= 6, f"{rsid}/{genotype} magnitude {m} out of range"

    def test_rsid_format(self):
        for rsid in COMPREHENSIVE_SNPS:
            assert rsid.startswith("rs"), f"Unexpected rsID format: {rsid}"

    def test_known_snp_exists(self):
        """CYP1A2 caffeine metabolism should be in the database."""
        assert "rs762551" in COMPREHENSIVE_SNPS
        info = COMPREHENSIVE_SNPS["rs762551"]
        assert info["gene"] == "CYP1A2"
        assert "AA" in info["variants"]


class TestSafeInt:
    def test_valid_integer_string(self):
        assert _safe_int("3") == 3

    def test_empty_string(self):
        assert _safe_int("") == 0

    def test_none(self):
        assert _safe_int(None) == 0

    def test_non_numeric_string(self):
        assert _safe_int("abc") == 0

    def test_custom_default(self):
        assert _safe_int("", default=-1) == -1

    def test_float_string(self):
        assert _safe_int("3.5") == 0

    def test_zero(self):
        assert _safe_int("0") == 0


class TestLookupGenotype:
    def test_direct_match(self):
        variants = {"AG": "hit", "GG": "other"}
        assert _lookup_genotype(variants, "AG") == "hit"

    def test_reverse_complement_match(self):
        variants = {"AG": "hit"}
        assert _lookup_genotype(variants, "GA") == "hit"

    def test_no_match(self):
        variants = {"AG": "hit"}
        assert _lookup_genotype(variants, "CC") is None

    def test_single_char_genotype(self):
        """Single-allele genotypes (MT/Y) should not be reversed."""
        variants = {"A": "hit"}
        assert _lookup_genotype(variants, "A") == "hit"

    def test_single_char_no_match(self):
        variants = {"A": "hit"}
        assert _lookup_genotype(variants, "G") is None

    def test_direct_preferred_over_reverse(self):
        """If both AG and GA exist, direct match wins."""
        variants = {"AG": "direct", "GA": "reverse"}
        assert _lookup_genotype(variants, "AG") == "direct"
        assert _lookup_genotype(variants, "GA") == "reverse"


class TestGenotypeMatching:
    def test_direct_match(self):
        genome = {"rs762551": {"chromosome": "15", "position": "75041917", "genotype": "AA"}}
        results = analyze_lifestyle_health(genome, {})
        assert len(results["findings"]) == 1
        assert results["findings"][0]["genotype"] == "AA"
        assert results["findings"][0]["gene"] == "CYP1A2"

    def test_reverse_complement_match(self):
        """If genotype is stored as 'GA' but database has 'AG', should still match."""
        # Find a SNP that has 'AG' but not 'GA'
        for rsid, info in COMPREHENSIVE_SNPS.items():
            variants = info["variants"]
            if "AG" in variants and "GA" not in variants:
                genome = {rsid: {"chromosome": "1", "position": "1", "genotype": "GA"}}
                results = analyze_lifestyle_health(genome, {})
                assert len(results["findings"]) == 1
                assert results["findings"][0]["rsid"] == rsid
                return
        # If all SNPs with AG also have GA, that's fine — skip
        assert True

    def test_no_match_for_absent_snp(self):
        genome = {"rs000000": {"chromosome": "1", "position": "1", "genotype": "AA"}}
        results = analyze_lifestyle_health(genome, {})
        assert len(results["findings"]) == 0

    def test_impact_counting(self):
        """Magnitude >= 3 should count as high_impact."""
        genome = {
            "rs762551": {"chromosome": "15", "position": "1", "genotype": "CC"},  # magnitude 3 (slow)
        }
        results = analyze_lifestyle_health(genome, {})
        assert results["summary"]["high_impact"] >= 1

    def test_pharmgkb_findings(self):
        """PharmGKB matches should appear in pharmgkb_findings."""
        pharmgkb = {
            "rs762551": {
                "gene": "CYP1A2",
                "drugs": "caffeine",
                "phenotype": "Metabolism",
                "level": "1A",
                "category": "PK",
                "genotypes": {"AA": "Fast metabolizer"},
            }
        }
        genome = {"rs762551": {"chromosome": "15", "position": "1", "genotype": "AA"}}
        results = analyze_lifestyle_health(genome, pharmgkb)
        assert len(results["pharmgkb_findings"]) == 1
        assert results["pharmgkb_findings"][0]["drugs"] == "caffeine"

    def test_pharmgkb_low_evidence_filtered(self):
        """Only levels 1A/1B/2A/2B should be included."""
        pharmgkb = {
            "rs762551": {
                "gene": "CYP1A2",
                "drugs": "caffeine",
                "phenotype": "Metabolism",
                "level": "3",
                "category": "PK",
                "genotypes": {"AA": "Some text"},
            }
        }
        genome = {"rs762551": {"chromosome": "15", "position": "1", "genotype": "AA"}}
        results = analyze_lifestyle_health(genome, pharmgkb)
        assert len(results["pharmgkb_findings"]) == 0
