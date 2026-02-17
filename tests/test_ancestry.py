"""Tests for ancestry estimation module."""

from genetic_health.ancestry import (
    ANCESTRY_MARKERS,
    POPULATIONS,
    POPULATION_LABELS,
    POPULATION_NOTES,
    estimate_ancestry,
    get_population_warnings,
    _count_allele,
    _softmax,
)


class TestAncestryMarkers:
    def test_all_markers_have_required_keys(self):
        for rsid, info in ANCESTRY_MARKERS.items():
            assert "gene" in info, f"{rsid} missing 'gene'"
            assert "description" in info, f"{rsid} missing 'description'"
            assert "allele" in info, f"{rsid} missing 'allele'"
            assert "frequencies" in info, f"{rsid} missing 'frequencies'"

    def test_all_markers_have_all_populations(self):
        for rsid, info in ANCESTRY_MARKERS.items():
            for pop in POPULATIONS:
                assert pop in info["frequencies"], f"{rsid} missing freq for {pop}"

    def test_frequencies_in_valid_range(self):
        for rsid, info in ANCESTRY_MARKERS.items():
            for pop, freq in info["frequencies"].items():
                assert 0.0 <= freq <= 1.0, f"{rsid}/{pop} freq {freq} out of range"

    def test_allele_is_single_char(self):
        for rsid, info in ANCESTRY_MARKERS.items():
            assert len(info["allele"]) == 1, f"{rsid} allele should be single char"

    def test_known_marker_exists(self):
        assert "rs1426654" in ANCESTRY_MARKERS
        info = ANCESTRY_MARKERS["rs1426654"]
        assert info["gene"] == "SLC24A5"
        assert info["frequencies"]["EUR"] > 0.9
        assert info["frequencies"]["AFR"] < 0.1

    def test_marker_count(self):
        assert len(ANCESTRY_MARKERS) >= 50, "Should have at least 50 AIMs"

    def test_population_labels_complete(self):
        for pop in POPULATIONS:
            assert pop in POPULATION_LABELS


class TestCountAllele:
    def test_homozygous_match(self):
        assert _count_allele("AA", "A") == 2

    def test_heterozygous(self):
        assert _count_allele("AG", "A") == 1

    def test_no_match(self):
        assert _count_allele("GG", "A") == 0

    def test_single_allele(self):
        assert _count_allele("A", "A") == 1

    def test_single_no_match(self):
        assert _count_allele("G", "A") == 0


class TestSoftmax:
    def test_equal_likelihoods_give_uniform(self):
        lls = {"A": 0.0, "B": 0.0, "C": 0.0}
        result = _softmax(lls)
        for v in result.values():
            assert abs(v - 1.0 / 3) < 1e-6

    def test_sums_to_one(self):
        lls = {"A": -10, "B": -5, "C": -20}
        result = _softmax(lls)
        assert abs(sum(result.values()) - 1.0) < 1e-6

    def test_highest_wins(self):
        lls = {"A": -1, "B": -100, "C": -100}
        result = _softmax(lls)
        assert result["A"] > result["B"]
        assert result["A"] > result["C"]


class TestEstimateAncestry:
    def _make_genome(self, rsid_genotypes):
        """Helper to create a genome dict."""
        return {
            rsid: {"chromosome": "1", "position": "100", "genotype": gt}
            for rsid, gt in rsid_genotypes.items()
        }

    def test_empty_genome_returns_uniform(self):
        result = estimate_ancestry({})
        assert result["markers_found"] == 0
        assert result["confidence"] == "none"
        for pop in POPULATIONS:
            assert abs(result["proportions"][pop] - 0.2) < 1e-6

    def test_european_markers(self):
        """SLC24A5 AA + HERC2 GG should strongly indicate European."""
        genome = self._make_genome({
            "rs1426654": "AA",   # SLC24A5 - EUR ~0.98
            "rs12913832": "GG",  # HERC2 - EUR ~0.72
            "rs16891982": "GG",  # SLC45A2 - EUR ~0.94
            "rs1042602": "AA",   # TYR - EUR ~0.38, very rare elsewhere
            "rs1805007": "CT",   # MC1R - EUR ~0.10, rare elsewhere
        })
        result = estimate_ancestry(genome)
        assert result["proportions"]["EUR"] > 0.5
        assert result["top_ancestry"] == "European"

    def test_east_asian_markers(self):
        """EDAR AA + ALDH2 markers should indicate East Asian."""
        genome = self._make_genome({
            "rs3827760": "AA",    # EDAR - EAS ~0.87
            "rs1229984": "TT",    # ADH1B - EAS ~0.70
            "rs1800414": "CC",    # OCA2 - EAS ~0.63
            "rs17822931": "TT",   # ABCC11 - EAS ~0.92
        })
        result = estimate_ancestry(genome)
        assert result["proportions"]["EAS"] > 0.5
        assert result["top_ancestry"] == "East Asian"

    def test_african_markers(self):
        """DARC + MFSD12 markers should indicate African."""
        genome = self._make_genome({
            "rs2814778": "TT",    # DARC - AFR has very low C freq
            "rs6497268": "TT",    # MFSD12 - AFR ~0.53
            "rs2065160": "CC",    # LOXL1 - AFR ~0.90
            "rs174537": "TT",     # FADS1 - AFR ~0.96
            "rs2250072": "TT",    # HBB - AFR ~0.45
        })
        result = estimate_ancestry(genome)
        assert result["proportions"]["AFR"] > 0.3

    def test_proportions_sum_to_one(self):
        genome = self._make_genome({"rs1426654": "AG", "rs762551": "AC"})
        result = estimate_ancestry(genome)
        assert abs(sum(result["proportions"].values()) - 1.0) < 1e-6

    def test_confidence_high_with_many_markers(self):
        """40+ markers should give high confidence."""
        genome = {}
        for rsid, info in list(ANCESTRY_MARKERS.items())[:45]:
            allele = info["allele"]
            genome[rsid] = {
                "chromosome": "1", "position": "100",
                "genotype": allele + allele,
            }
        result = estimate_ancestry(genome)
        assert result["confidence"] == "high"
        assert result["markers_found"] >= 40

    def test_confidence_moderate(self):
        """20-39 markers should give moderate confidence."""
        genome = {}
        for rsid, info in list(ANCESTRY_MARKERS.items())[:25]:
            allele = info["allele"]
            genome[rsid] = {
                "chromosome": "1", "position": "100",
                "genotype": allele + allele,
            }
        result = estimate_ancestry(genome)
        assert result["confidence"] == "moderate"

    def test_confidence_low(self):
        """<20 markers should give low confidence."""
        genome = self._make_genome({"rs1426654": "AA"})
        result = estimate_ancestry(genome)
        assert result["confidence"] == "low"

    def test_details_populated(self):
        genome = self._make_genome({"rs1426654": "AG", "rs762551": "AA"})
        result = estimate_ancestry(genome)
        assert len(result["details"]) == 2
        detail = result["details"][0]
        assert "rsid" in detail
        assert "gene" in detail
        assert "allele_count" in detail

    def test_result_keys(self):
        result = estimate_ancestry({})
        assert "proportions" in result
        assert "markers_found" in result
        assert "confidence" in result
        assert "top_ancestry" in result
        assert "details" in result


class TestPopulationWarnings:
    def test_known_warning(self):
        warnings = get_population_warnings("ALDH2", "reduced")
        assert len(warnings) == 1
        assert "East Asian" in warnings[0]

    def test_no_warning(self):
        warnings = get_population_warnings("UNKNOWN_GENE", "unknown_status")
        assert warnings == []

    def test_hfe_carrier_warning(self):
        warnings = get_population_warnings("HFE", "carrier")
        assert len(warnings) == 1
        assert "European" in warnings[0]
