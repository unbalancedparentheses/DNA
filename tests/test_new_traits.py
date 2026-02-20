"""Tests for the 10 new trait predictions added to traits.py."""

from genetic_health.traits import predict_traits


def _genome(**snps):
    return {rsid: {"genotype": gt, "chromosome": "1", "position": "1"}
            for rsid, gt in snps.items()}


class TestLactoseTolerance:
    def test_tolerant_aa(self):
        r = predict_traits(_genome(rs4988235="AA"))
        assert "tolerant" in r["lactose_tolerance"]["prediction"].lower()
        assert r["lactose_tolerance"]["confidence"] == "high"

    def test_intolerant_gg(self):
        r = predict_traits(_genome(rs4988235="GG"))
        assert "intolerant" in r["lactose_tolerance"]["prediction"].lower()

    def test_heterozygous(self):
        r = predict_traits(_genome(rs4988235="AG"))
        assert "tolerant" in r["lactose_tolerance"]["prediction"].lower()

    def test_missing(self):
        r = predict_traits({})
        assert r["lactose_tolerance"]["prediction"] == "Unknown"


class TestBitterTaste:
    def test_supertaster(self):
        r = predict_traits(_genome(rs713598="GG", rs1726866="CC", rs10246939="CC"))
        assert "strong" in r["bitter_taste"]["prediction"].lower() or "super" in r["bitter_taste"]["prediction"].lower()

    def test_nontaster(self):
        r = predict_traits(_genome(rs713598="CC", rs1726866="TT", rs10246939="TT"))
        assert "non-taster" in r["bitter_taste"]["prediction"].lower()

    def test_partial_snps(self):
        r = predict_traits(_genome(rs713598="GC"))
        assert r["bitter_taste"]["prediction"] != "Unknown"

    def test_missing(self):
        r = predict_traits({})
        assert r["bitter_taste"]["prediction"] == "Unknown"


class TestCilantroTaste:
    def test_soapy_cc(self):
        r = predict_traits(_genome(rs72921001="CC"))
        assert "soapy" in r["cilantro_taste"]["prediction"].lower()

    def test_normal_tt(self):
        r = predict_traits(_genome(rs72921001="TT"))
        assert "normal" in r["cilantro_taste"]["prediction"].lower()

    def test_carrier(self):
        r = predict_traits(_genome(rs72921001="CT"))
        assert "mild" in r["cilantro_taste"]["prediction"].lower() or "possible" in r["cilantro_taste"]["prediction"].lower()

    def test_missing(self):
        r = predict_traits({})
        assert r["cilantro_taste"]["prediction"] == "Unknown"


class TestAsparaguSmell:
    def test_can_smell(self):
        r = predict_traits(_genome(rs4481887="GG"))
        assert "can smell" in r["asparagus_smell"]["prediction"].lower()

    def test_cannot_smell(self):
        r = predict_traits(_genome(rs4481887="AA"))
        assert "cannot" in r["asparagus_smell"]["prediction"].lower()

    def test_missing(self):
        r = predict_traits({})
        assert r["asparagus_smell"]["prediction"] == "Unknown"


class TestMuscleFiberType:
    def test_power_cc(self):
        r = predict_traits(_genome(rs1815739="CC"))
        assert "power" in r["muscle_fiber_type"]["prediction"].lower()

    def test_endurance_tt(self):
        r = predict_traits(_genome(rs1815739="TT"))
        assert "endurance" in r["muscle_fiber_type"]["prediction"].lower()

    def test_mixed_ct(self):
        r = predict_traits(_genome(rs1815739="CT"))
        assert "mixed" in r["muscle_fiber_type"]["prediction"].lower()

    def test_missing(self):
        r = predict_traits({})
        assert r["muscle_fiber_type"]["prediction"] == "Unknown"


class TestSecretorStatus:
    def test_nonsecretor_aa(self):
        r = predict_traits(_genome(rs601338="AA"))
        assert "non-secretor" in r["secretor_status"]["prediction"].lower()

    def test_secretor_gg(self):
        r = predict_traits(_genome(rs601338="GG"))
        assert "secretor" in r["secretor_status"]["prediction"].lower()
        assert "non" not in r["secretor_status"]["prediction"].lower()

    def test_carrier_ga(self):
        r = predict_traits(_genome(rs601338="GA"))
        assert "secretor" in r["secretor_status"]["prediction"].lower()

    def test_missing(self):
        r = predict_traits({})
        assert r["secretor_status"]["prediction"] == "Unknown"


class TestPhoticSneeze:
    def test_sneezer(self):
        r = predict_traits(_genome(rs10427255="CC"))
        assert "sneezer" in r["photic_sneeze"]["prediction"].lower()

    def test_not_sneezer(self):
        r = predict_traits(_genome(rs10427255="TT"))
        assert "unlikely" in r["photic_sneeze"]["prediction"].lower()

    def test_missing(self):
        r = predict_traits({})
        assert r["photic_sneeze"]["prediction"] == "Unknown"


class TestHairCurl:
    def test_straight_aa(self):
        r = predict_traits(_genome(rs11803731="AA"))
        assert "straight" in r["hair_curl"]["prediction"].lower()

    def test_curly_tt(self):
        r = predict_traits(_genome(rs11803731="TT"))
        assert "curly" in r["hair_curl"]["prediction"].lower() or "wavy" in r["hair_curl"]["prediction"].lower()

    def test_missing(self):
        r = predict_traits({})
        assert r["hair_curl"]["prediction"] == "Unknown"


class TestBaldnessRisk:
    def test_higher_cc(self):
        r = predict_traits(_genome(rs2180439="CC"))
        assert "higher" in r["baldness_risk"]["prediction"].lower()

    def test_lower_tt(self):
        r = predict_traits(_genome(rs2180439="TT"))
        assert "lower" in r["baldness_risk"]["prediction"].lower()

    def test_missing(self):
        r = predict_traits({})
        assert r["baldness_risk"]["prediction"] == "Unknown"


class TestUnibrowTendency:
    def test_higher_cc(self):
        r = predict_traits(_genome(rs12651896="CC"))
        assert "higher" in r["unibrow_tendency"]["prediction"].lower()

    def test_lower_tt(self):
        r = predict_traits(_genome(rs12651896="TT"))
        assert "lower" in r["unibrow_tendency"]["prediction"].lower()

    def test_missing(self):
        r = predict_traits({})
        assert r["unibrow_tendency"]["prediction"] == "Unknown"


class TestAllTraitsPresent:
    def test_all_14_traits_returned(self):
        r = predict_traits({})
        expected = {
            "eye_color", "hair_color", "earwax_type", "freckling",
            "lactose_tolerance", "bitter_taste", "cilantro_taste",
            "asparagus_smell", "muscle_fiber_type", "secretor_status",
            "photic_sneeze", "hair_curl", "baldness_risk", "unibrow_tendency",
        }
        assert set(r.keys()) == expected

    def test_all_have_required_keys(self):
        r = predict_traits({})
        for trait_id, trait in r.items():
            assert "prediction" in trait, f"{trait_id} missing prediction"
            assert "confidence" in trait, f"{trait_id} missing confidence"
            assert "snps_used" in trait, f"{trait_id} missing snps_used"
            assert "description" in trait, f"{trait_id} missing description"
