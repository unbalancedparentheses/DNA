"""Tests for mitochondrial haplogroup estimation module."""

from genetic_health.mt_haplogroup import (
    MT_HAPLOGROUP_TREE,
    estimate_mt_haplogroup,
)


def _make_genome(rsid_genotypes):
    """Helper to create genome dict."""
    return {
        rsid: {"chromosome": "MT", "position": "100", "genotype": gt}
        for rsid, gt in rsid_genotypes.items()
    }


class TestMTHaplogroupTree:
    def test_all_entries_have_required_keys(self):
        for entry in MT_HAPLOGROUP_TREE:
            assert "rsid" in entry
            assert "allele" in entry
            assert "haplogroup" in entry
            assert "description" in entry

    def test_alleles_are_single_char(self):
        for entry in MT_HAPLOGROUP_TREE:
            assert len(entry["allele"]) == 1

    def test_minimum_markers(self):
        assert len(MT_HAPLOGROUP_TREE) >= 25


class TestEstimateMTHaplogroup:
    def test_european_haplogroup_h(self):
        genome = _make_genome({"rs2032658": "G"})
        result = estimate_mt_haplogroup(genome)
        assert result["haplogroup"] == "H"
        assert "European" in result["lineage"]

    def test_african_l2(self):
        genome = _make_genome({"rs28358571": "T"})
        result = estimate_mt_haplogroup(genome)
        assert result["haplogroup"] == "L2"
        assert "African" in result["lineage"]

    def test_east_asian_haplogroup_d(self):
        genome = _make_genome({"rs28357968": "T"})
        result = estimate_mt_haplogroup(genome)
        assert result["haplogroup"] == "D"

    def test_most_specific_match_wins(self):
        """H1 is more specific than H — should pick H1."""
        genome = _make_genome({
            "rs2032658": "G",   # H
            "rs2853825": "T",   # H1 (more specific)
        })
        result = estimate_mt_haplogroup(genome)
        assert result["haplogroup"] == "H1"

    def test_no_mt_snps_graceful_fallback(self):
        result = estimate_mt_haplogroup({})
        assert result["haplogroup"] == "Unknown"
        assert result["confidence"] == "none"
        assert result["markers_found"] == 0

    def test_low_confidence_with_few_markers(self):
        genome = _make_genome({"rs2032658": "G"})
        result = estimate_mt_haplogroup(genome)
        assert result["confidence"] == "low"

    def test_moderate_confidence(self):
        # Provide 5+ markers that exist (but not all match)
        genome = {}
        for i, entry in enumerate(MT_HAPLOGROUP_TREE[:8]):
            genome[entry["rsid"]] = {
                "chromosome": "MT",
                "position": str(i),
                "genotype": entry["allele"],
            }
        result = estimate_mt_haplogroup(genome)
        assert result["confidence"] in ("moderate", "high")
        assert result["markers_found"] >= 5

    def test_result_keys(self):
        result = estimate_mt_haplogroup({})
        expected_keys = {
            "haplogroup", "description", "confidence",
            "markers_found", "markers_tested", "lineage", "details",
        }
        assert set(result.keys()) == expected_keys

    def test_details_populated(self):
        genome = _make_genome({"rs2032658": "G", "rs2853825": "T"})
        result = estimate_mt_haplogroup(genome)
        assert len(result["details"]) >= 1
        detail = result["details"][0]
        assert "rsid" in detail
        assert "haplogroup" in detail
        assert "allele" in detail

    def test_non_matching_allele_not_counted(self):
        """If the genotype doesn't contain the defining allele, no match."""
        genome = _make_genome({"rs2032658": "A"})  # Defining allele is G
        result = estimate_mt_haplogroup(genome)
        assert result["haplogroup"] == "Unknown"
