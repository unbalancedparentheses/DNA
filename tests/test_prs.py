"""Tests for Polygenic Risk Score module."""

from genetic_health.prs import (
    PRS_MODELS,
    calculate_prs,
    _count_risk_allele,
    _z_to_percentile,
    _categorize_percentile,
)


class TestPRSModels:
    def test_all_models_have_required_keys(self):
        for cid, model in PRS_MODELS.items():
            assert "name" in model, f"{cid} missing 'name'"
            assert "reference" in model, f"{cid} missing 'reference'"
            assert "snps" in model, f"{cid} missing 'snps'"
            assert len(model["snps"]) >= 15, f"{cid} has too few SNPs"

    def test_all_snps_have_required_fields(self):
        for cid, model in PRS_MODELS.items():
            for snp in model["snps"]:
                assert "rsid" in snp, f"{cid} snp missing 'rsid'"
                assert "risk_allele" in snp, f"{cid}/{snp.get('rsid')} missing 'risk_allele'"
                assert "log_or" in snp, f"{cid}/{snp.get('rsid')} missing 'log_or'"
                assert "gene" in snp, f"{cid}/{snp.get('rsid')} missing 'gene'"
                assert "eur_freq" in snp, f"{cid}/{snp.get('rsid')} missing 'eur_freq'"

    def test_log_or_positive(self):
        for cid, model in PRS_MODELS.items():
            for snp in model["snps"]:
                assert snp["log_or"] > 0, f"{cid}/{snp['rsid']} log_or should be positive"

    def test_eur_freq_in_range(self):
        for cid, model in PRS_MODELS.items():
            for snp in model["snps"]:
                assert 0.0 < snp["eur_freq"] < 1.0, \
                    f"{cid}/{snp['rsid']} eur_freq out of range"

    def test_risk_allele_single_char(self):
        for cid, model in PRS_MODELS.items():
            for snp in model["snps"]:
                assert len(snp["risk_allele"]) == 1, \
                    f"{cid}/{snp['rsid']} risk_allele should be single char"

    def test_five_conditions_present(self):
        expected = [
            "type2_diabetes",
            "coronary_artery_disease",
            "hypertension",
            "breast_cancer",
            "age_related_macular_degeneration",
        ]
        for cid in expected:
            assert cid in PRS_MODELS, f"Missing condition: {cid}"

    def test_rsid_format(self):
        for cid, model in PRS_MODELS.items():
            for snp in model["snps"]:
                assert snp["rsid"].startswith("rs"), \
                    f"{cid}/{snp['rsid']} invalid rsID format"


class TestCountRiskAllele:
    def test_homozygous_risk(self):
        assert _count_risk_allele("TT", "T") == 2

    def test_heterozygous(self):
        assert _count_risk_allele("AT", "T") == 1

    def test_no_risk(self):
        assert _count_risk_allele("AA", "T") == 0

    def test_single_allele(self):
        assert _count_risk_allele("T", "T") == 1


class TestZToPercentile:
    def test_zero_gives_50(self):
        assert abs(_z_to_percentile(0.0) - 50.0) < 0.1

    def test_positive_above_50(self):
        assert _z_to_percentile(1.0) > 50.0

    def test_negative_below_50(self):
        assert _z_to_percentile(-1.0) < 50.0

    def test_extreme_positive(self):
        assert _z_to_percentile(3.0) > 99.0

    def test_extreme_negative(self):
        assert _z_to_percentile(-3.0) < 1.0


class TestCategorizePercentile:
    def test_low(self):
        assert _categorize_percentile(10.0) == "low"

    def test_average(self):
        assert _categorize_percentile(50.0) == "average"

    def test_elevated(self):
        assert _categorize_percentile(85.0) == "elevated"

    def test_high(self):
        assert _categorize_percentile(96.0) == "high"

    def test_boundaries(self):
        assert _categorize_percentile(19.9) == "low"
        assert _categorize_percentile(20.0) == "average"
        assert _categorize_percentile(79.9) == "average"
        assert _categorize_percentile(80.0) == "elevated"
        assert _categorize_percentile(94.9) == "elevated"
        assert _categorize_percentile(95.0) == "high"


class TestCalculatePRS:
    def _make_genome(self, rsid_genotypes):
        return {
            rsid: {"chromosome": "1", "position": "100", "genotype": gt}
            for rsid, gt in rsid_genotypes.items()
        }

    def test_empty_genome(self):
        results = calculate_prs({})
        assert len(results) == len(PRS_MODELS)
        for cid, r in results.items():
            assert r["snps_found"] == 0
            assert r["percentile"] == 50.0

    def test_result_keys(self):
        results = calculate_prs({})
        for cid, r in results.items():
            assert "name" in r
            assert "raw_score" in r
            assert "z_score" in r
            assert "percentile" in r
            assert "risk_category" in r
            assert "snps_found" in r
            assert "snps_total" in r
            assert "ancestry_applicable" in r
            assert "contributing_snps" in r
            assert "reference" in r

    def test_homozygous_risk_alleles_increase_score(self):
        """Loading up on risk alleles should increase raw score."""
        # Use T2D model, give homozygous risk for TCF7L2
        genome = self._make_genome({"rs7903146": "TT"})
        results = calculate_prs(genome)
        r = results["type2_diabetes"]
        assert r["snps_found"] == 1
        assert r["raw_score"] > 0

    def test_no_risk_alleles_low_score(self):
        """No risk alleles should give below-average score."""
        genome = self._make_genome({"rs7903146": "CC"})
        results = calculate_prs(genome)
        r = results["type2_diabetes"]
        assert r["snps_found"] == 1
        # raw_score should be 0 (no risk alleles counted)
        assert r["raw_score"] == 0.0

    def test_ancestry_warning_non_european(self):
        """Non-EUR ancestry should flag models as not fully applicable."""
        genome = self._make_genome({"rs7903146": "TT"})
        ancestry = {"EUR": 0.3, "AFR": 0.5, "EAS": 0.1, "SAS": 0.05, "AMR": 0.05}
        results = calculate_prs(genome, ancestry_proportions=ancestry)
        for r in results.values():
            assert r["ancestry_applicable"] is False
            assert "EUR ancestry" in r["ancestry_warning"]

    def test_ancestry_ok_european(self):
        """EUR-dominant ancestry should keep models as applicable."""
        genome = self._make_genome({"rs7903146": "TT"})
        ancestry = {"EUR": 0.85, "AFR": 0.05, "EAS": 0.03, "SAS": 0.04, "AMR": 0.03}
        results = calculate_prs(genome, ancestry_proportions=ancestry)
        for r in results.values():
            assert r["ancestry_applicable"] is True

    def test_contributing_snps_populated(self):
        """Contributing SNPs should include details when risk alleles found."""
        genome = self._make_genome({
            "rs7903146": "TT",
            "rs5219": "TT",
        })
        results = calculate_prs(genome)
        r = results["type2_diabetes"]
        assert len(r["contributing_snps"]) == 2
        snp = r["contributing_snps"][0]
        assert "rsid" in snp
        assert "gene" in snp
        assert "copies" in snp
        assert "contribution" in snp

    def test_contributing_sorted_by_contribution(self):
        """Contributing SNPs should be sorted by contribution descending."""
        genome = self._make_genome({
            "rs7903146": "TT",  # log_or 0.322 * 2 = 0.644
            "rs5219": "TT",    # log_or 0.148 * 2 = 0.296
        })
        results = calculate_prs(genome)
        contribs = results["type2_diabetes"]["contributing_snps"]
        for i in range(len(contribs) - 1):
            assert contribs[i]["contribution"] >= contribs[i + 1]["contribution"]

    def test_multiple_conditions_scored(self):
        """A SNP shared between models should contribute to both."""
        # rs3184504 appears in both CAD and hypertension
        genome = self._make_genome({"rs3184504": "TT"})
        results = calculate_prs(genome)
        assert results["coronary_artery_disease"]["snps_found"] >= 1
        assert results["hypertension"]["snps_found"] >= 1

    def test_percentile_range(self):
        """Percentile should always be between 0 and 100."""
        # Give a mix of genotypes
        genome = self._make_genome({
            "rs7903146": "TT",
            "rs1111875": "CC",
            "rs10811661": "TT",
        })
        results = calculate_prs(genome)
        for r in results.values():
            assert 0.0 <= r["percentile"] <= 100.0
