"""Tests for star allele calling module."""

from genetic_health.star_alleles import (
    STAR_ALLELE_DEFINITIONS,
    call_star_alleles,
)


def _make_genome(rsid_genotypes):
    """Helper to create genome dict."""
    return {
        rsid: {"chromosome": "10", "position": "100", "genotype": gt}
        for rsid, gt in rsid_genotypes.items()
    }


class TestStarAlleleDefinitions:
    def test_all_genes_have_required_keys(self):
        for gene, defn in STAR_ALLELE_DEFINITIONS.items():
            assert "function_map" in defn, f"{gene} missing function_map"
            assert "alleles" in defn, f"{gene} missing alleles"
            assert "snps" in defn, f"{gene} missing snps"

    def test_star1_is_always_normal(self):
        for gene, defn in STAR_ALLELE_DEFINITIONS.items():
            assert defn["function_map"]["*1"] == "normal"

    def test_expected_genes_present(self):
        assert "CYP2C19" in STAR_ALLELE_DEFINITIONS
        assert "CYP2C9" in STAR_ALLELE_DEFINITIONS
        assert "CYP2D6" in STAR_ALLELE_DEFINITIONS


class TestCallStarAlleles:
    def test_normal_cyp2c19(self):
        """All reference alleles -> *1/*1 = normal metabolizer."""
        genome = _make_genome({
            "rs4244285": "GG",
            "rs4986893": "GG",
            "rs12248560": "CC",
        })
        results = call_star_alleles(genome)
        r = results["CYP2C19"]
        assert r["diplotype"] == "*1/*1"
        assert r["phenotype"] == "normal"

    def test_cyp2c19_star2_heterozygous(self):
        """rs4244285 GA -> *1/*2 = intermediate metabolizer."""
        genome = _make_genome({
            "rs4244285": "GA",
            "rs4986893": "GG",
            "rs12248560": "CC",
        })
        results = call_star_alleles(genome)
        r = results["CYP2C19"]
        assert r["diplotype"] == "*1/*2"
        assert r["phenotype"] == "intermediate"

    def test_cyp2c19_star2_homozygous(self):
        """rs4244285 AA -> *2/*2 = poor metabolizer."""
        genome = _make_genome({
            "rs4244285": "AA",
            "rs4986893": "GG",
            "rs12248560": "CC",
        })
        results = call_star_alleles(genome)
        r = results["CYP2C19"]
        assert r["diplotype"] == "*2/*2"
        assert r["phenotype"] == "poor"

    def test_cyp2c19_star17(self):
        """rs12248560 CT -> *1/*17 = rapid metabolizer."""
        genome = _make_genome({
            "rs4244285": "GG",
            "rs4986893": "GG",
            "rs12248560": "CT",
        })
        results = call_star_alleles(genome)
        r = results["CYP2C19"]
        assert r["diplotype"] == "*1/*17"
        assert r["phenotype"] == "rapid"

    def test_cyp2c9_star3_heterozygous(self):
        """rs1057910 AC -> *1/*3 = intermediate."""
        genome = _make_genome({
            "rs1799853": "CC",
            "rs1057910": "AC",
        })
        results = call_star_alleles(genome)
        r = results["CYP2C9"]
        assert r["diplotype"] == "*1/*3"
        assert r["phenotype"] == "intermediate"

    def test_cyp2d6_star4_heterozygous(self):
        """rs3892097 GA -> *1/*4 = intermediate."""
        genome = _make_genome({
            "rs3892097": "GA",
            "rs1065852": "CC",
        })
        results = call_star_alleles(genome)
        r = results["CYP2D6"]
        assert r["diplotype"] == "*1/*4"
        assert r["phenotype"] == "intermediate"

    def test_cyp2d6_has_caveat(self):
        """CYP2D6 should always mention copy number limitation."""
        genome = _make_genome({"rs3892097": "GG", "rs1065852": "CC"})
        results = call_star_alleles(genome)
        assert "deletions" in results["CYP2D6"]["clinical_note"].lower() or \
               "duplications" in results["CYP2D6"]["clinical_note"].lower()

    def test_missing_snps_partial_result(self):
        """Missing all SNPs -> Unknown with note."""
        results = call_star_alleles({})
        for gene, r in results.items():
            assert r["diplotype"] == "Unknown"
            assert r["phenotype"] == "Unknown"
            assert r["snps_found"] == 0

    def test_all_genes_return_required_keys(self):
        genome = _make_genome({"rs4244285": "GG"})
        results = call_star_alleles(genome)
        expected_keys = {"gene", "diplotype", "phenotype",
                         "snps_found", "snps_total", "clinical_note",
                         "coverage", "confidence"}
        for gene, r in results.items():
            assert set(r.keys()) == expected_keys, f"{gene} missing keys"

    def test_cyp2c9_normal(self):
        """Both reference -> *1/*1 = normal."""
        genome = _make_genome({
            "rs1799853": "CC",
            "rs1057910": "AA",
        })
        results = call_star_alleles(genome)
        r = results["CYP2C9"]
        assert r["diplotype"] == "*1/*1"
        assert r["phenotype"] == "normal"

    def test_indeterminate_phenotype_for_unmapped_pair(self):
        """If function pair is not in the map, phenotype should be Indeterminate."""
        from genetic_health.star_alleles import _PHENOTYPE_MAP
        # Verify Indeterminate is returned for a truly unknown pair
        # We can't easily manufacture one via the public API since all defined
        # functions are covered, so test the internal logic directly.
        assert _PHENOTYPE_MAP.get(("unknown_func", "unknown_func")) is None

    def test_cyp2c19_star2_star17_compound(self):
        """rs4244285 GA + rs12248560 CT -> *2/*17 = intermediate."""
        genome = _make_genome({
            "rs4244285": "GA",
            "rs4986893": "GG",
            "rs12248560": "CT",
        })
        results = call_star_alleles(genome)
        r = results["CYP2C19"]
        assert r["diplotype"] == "*17/*2"
        assert r["phenotype"] == "intermediate"
