"""Visible trait predictions from well-established SNP associations.

Predicts eye color, hair color, earwax type, and freckling/sun sensitivity
from SNPs that are typically available in 23andMe data.
"""


def _get_genotype(genome_by_rsid, rsid):
    """Get genotype string for an rsID, or empty string if missing."""
    entry = genome_by_rsid.get(rsid, {})
    return entry.get("genotype", "")


def _predict_eye_color(genome_by_rsid):
    """Predict eye color from HERC2 (rs12913832) and OCA2 (rs1800407)."""
    herc2 = _get_genotype(genome_by_rsid, "rs12913832")
    oca2 = _get_genotype(genome_by_rsid, "rs1800407")

    snps_used = []
    if herc2:
        snps_used.append(f"rs12913832 (HERC2): {herc2}")
    if oca2:
        snps_used.append(f"rs1800407 (OCA2): {oca2}")

    if not herc2:
        return {
            "prediction": "Unknown",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "Primary eye color SNP (rs12913832) not available.",
        }

    # rs12913832: GG = likely brown, AG = green/hazel, AA = likely blue
    herc2_aa = herc2.count("A")

    if herc2_aa == 2:
        # AA — likely blue
        base = "blue"
        # OCA2 rs1800407 can modify toward green/hazel
        if oca2 and "A" in oca2:
            return {
                "prediction": "Blue or green",
                "confidence": "moderate",
                "snps_used": snps_used,
                "description": "HERC2 AA predicts blue, but OCA2 variant shifts toward green/hazel.",
            }
        return {
            "prediction": "Likely blue",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "HERC2 AA is the strongest predictor of blue eyes (~70-80% accuracy).",
        }
    elif herc2_aa == 1:
        # AG — green/hazel most common
        if oca2 and oca2.count("A") >= 1:
            return {
                "prediction": "Green or hazel",
                "confidence": "moderate",
                "snps_used": snps_used,
                "description": "Heterozygous HERC2 with OCA2 variant favors green/hazel.",
            }
        return {
            "prediction": "Green, hazel, or light brown",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "Heterozygous HERC2 — intermediate pigmentation. Most commonly green or hazel.",
        }
    else:
        # GG — likely brown
        return {
            "prediction": "Likely brown",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "HERC2 GG is strongly associated with brown eyes (~90% in most populations).",
        }


def _predict_hair_color(genome_by_rsid):
    """Predict hair color from MC1R variants (rs1805007, rs1805008)."""
    mc1r_7 = _get_genotype(genome_by_rsid, "rs1805007")
    mc1r_8 = _get_genotype(genome_by_rsid, "rs1805008")

    snps_used = []
    if mc1r_7:
        snps_used.append(f"rs1805007 (MC1R R151C): {mc1r_7}")
    if mc1r_8:
        snps_used.append(f"rs1805008 (MC1R R160W): {mc1r_8}")

    if not mc1r_7 and not mc1r_8:
        return {
            "prediction": "Unknown",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "MC1R SNPs not available for hair color prediction.",
        }

    # Count MC1R variant alleles (T alleles for both SNPs)
    red_alleles = 0
    if mc1r_7:
        red_alleles += mc1r_7.count("T")
    if mc1r_8:
        red_alleles += mc1r_8.count("T")

    if red_alleles >= 2:
        return {
            "prediction": "Likely red or auburn",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "Two or more MC1R loss-of-function alleles strongly predict red hair.",
        }
    elif red_alleles == 1:
        return {
            "prediction": "Possible red tint or auburn highlights",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "One MC1R variant — may have red highlights, strawberry blonde, or be a carrier without visible effect.",
        }
    else:
        return {
            "prediction": "Non-red (blonde, brown, or black — not determinable from MC1R alone)",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "No MC1R red hair variants detected. Hair color depends on many other genes.",
        }


def _predict_earwax(genome_by_rsid):
    """Predict earwax type from ABCC11 (rs17822931)."""
    abcc11 = _get_genotype(genome_by_rsid, "rs17822931")

    snps_used = []
    if abcc11:
        snps_used.append(f"rs17822931 (ABCC11): {abcc11}")

    if not abcc11:
        return {
            "prediction": "Unknown",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "ABCC11 SNP not available.",
        }

    t_count = abcc11.count("T")

    if t_count == 2:
        return {
            "prediction": "Dry earwax",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "TT genotype — dry, flaky earwax. Common in East Asian populations (~80-95%). Also associated with reduced body odor.",
        }
    elif t_count == 1:
        return {
            "prediction": "Wet earwax",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "CT genotype — wet earwax (C is dominant). One copy of the dry earwax allele.",
        }
    else:
        return {
            "prediction": "Wet earwax",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "CC genotype — wet, sticky earwax. Most common in European and African populations.",
        }


def _predict_freckling(genome_by_rsid):
    """Predict freckling/sun sensitivity from MC1R variants."""
    mc1r_7 = _get_genotype(genome_by_rsid, "rs1805007")
    mc1r_8 = _get_genotype(genome_by_rsid, "rs1805008")

    snps_used = []
    if mc1r_7:
        snps_used.append(f"rs1805007 (MC1R R151C): {mc1r_7}")
    if mc1r_8:
        snps_used.append(f"rs1805008 (MC1R R160W): {mc1r_8}")

    if not mc1r_7 and not mc1r_8:
        return {
            "prediction": "Unknown",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "MC1R SNPs not available for sun sensitivity prediction.",
        }

    variant_alleles = 0
    if mc1r_7:
        variant_alleles += mc1r_7.count("T")
    if mc1r_8:
        variant_alleles += mc1r_8.count("T")

    if variant_alleles >= 2:
        return {
            "prediction": "High freckling tendency, increased sun sensitivity",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "Multiple MC1R variants — higher UV sensitivity, increased freckling, and elevated skin cancer risk. Daily SPF recommended.",
        }
    elif variant_alleles == 1:
        return {
            "prediction": "Moderate freckling tendency",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "One MC1R variant — mildly increased sun sensitivity and freckling tendency.",
        }
    else:
        return {
            "prediction": "Typical sun sensitivity",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "No MC1R risk variants — typical melanin production. Sun protection still recommended.",
        }


def _predict_lactose(genome_by_rsid):
    """Predict lactose tolerance from MCM6/LCT (rs4988235)."""
    lct = _get_genotype(genome_by_rsid, "rs4988235")

    snps_used = []
    if lct:
        snps_used.append(f"rs4988235 (MCM6/LCT): {lct}")

    if not lct:
        return {
            "prediction": "Unknown",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "Lactase persistence SNP (rs4988235) not available.",
        }

    # rs4988235: AA = lactase persistent (tolerant), AG = likely tolerant, GG = likely intolerant
    a_count = lct.count("A")

    if a_count == 2:
        return {
            "prediction": "Lactose tolerant",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "AA genotype — you produce lactase into adulthood. Dairy digestion is normal. Common in Northern European descent.",
        }
    elif a_count == 1:
        return {
            "prediction": "Likely lactose tolerant",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "AG genotype — one copy of the lactase persistence allele. Most people with this genotype tolerate dairy well.",
        }
    else:
        return {
            "prediction": "Likely lactose intolerant",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "GG genotype — lactase production typically declines after childhood. You may experience bloating, gas, or diarrhea from dairy. Common in East Asian, African, and Southern European populations.",
        }


def _predict_bitter_taste(genome_by_rsid):
    """Predict bitter taste (PTC/PROP) sensitivity from TAS2R38 variants."""
    snp1 = _get_genotype(genome_by_rsid, "rs713598")
    snp2 = _get_genotype(genome_by_rsid, "rs1726866")
    snp3 = _get_genotype(genome_by_rsid, "rs10246939")

    snps_used = []
    if snp1:
        snps_used.append(f"rs713598 (TAS2R38 A49P): {snp1}")
    if snp2:
        snps_used.append(f"rs1726866 (TAS2R38 V262A): {snp2}")
    if snp3:
        snps_used.append(f"rs10246939 (TAS2R38 I296V): {snp3}")

    if not snp1 and not snp2 and not snp3:
        return {
            "prediction": "Unknown",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "TAS2R38 SNPs not available for bitter taste prediction.",
        }

    # Count taster alleles: G for rs713598, C for rs1726866, C for rs10246939
    # These form the PAV (taster) haplotype vs AVI (non-taster)
    taster_alleles = 0
    total_checked = 0
    if snp1:
        taster_alleles += snp1.count("G")
        total_checked += 2
    if snp2:
        taster_alleles += snp2.count("C")
        total_checked += 2
    if snp3:
        taster_alleles += snp3.count("C")
        total_checked += 2

    if total_checked == 0:
        ratio = 0
    else:
        ratio = taster_alleles / total_checked

    if ratio >= 0.8:
        return {
            "prediction": "Strong bitter taster (supertaster)",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "PAV/PAV haplotype — you strongly taste bitter compounds like PTC and PROP. Broccoli, Brussels sprouts, and dark coffee may taste more bitter to you.",
        }
    elif ratio >= 0.4:
        return {
            "prediction": "Medium bitter taster",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "PAV/AVI haplotype — you have intermediate bitter taste sensitivity. Most bitter foods are noticeable but not overwhelming.",
        }
    else:
        return {
            "prediction": "Non-taster (low bitter sensitivity)",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "AVI/AVI haplotype — you have low sensitivity to bitter compounds. You may enjoy bitter vegetables and black coffee more than average.",
        }


def _predict_cilantro(genome_by_rsid):
    """Predict cilantro/coriander taste aversion from OR6A2 (rs72921001)."""
    or6a2 = _get_genotype(genome_by_rsid, "rs72921001")

    snps_used = []
    if or6a2:
        snps_used.append(f"rs72921001 (OR6A2): {or6a2}")

    if not or6a2:
        return {
            "prediction": "Unknown",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "OR6A2 SNP not available for cilantro taste prediction.",
        }

    # rs72921001: CC = likely soapy taste, CT = mild aversion possible, TT = normal taste
    c_count = or6a2.count("C")

    if c_count == 2:
        return {
            "prediction": "Likely perceives cilantro as soapy",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "CC genotype at OR6A2 — associated with detecting aldehyde compounds in cilantro that taste soapy or unpleasant. ~15% of people of European descent report this.",
        }
    elif c_count == 1:
        return {
            "prediction": "Possible mild cilantro aversion",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "CT genotype — one copy of the aversion allele. Some people with this genotype notice a slight soapy taste; others do not.",
        }
    else:
        return {
            "prediction": "Normal cilantro taste (no soapy perception)",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "TT genotype — no OR6A2 aversion allele. Cilantro likely tastes herbal and pleasant to you.",
        }


def _predict_asparagus_smell(genome_by_rsid):
    """Predict ability to smell asparagus metabolites from rs4481887."""
    snp = _get_genotype(genome_by_rsid, "rs4481887")

    snps_used = []
    if snp:
        snps_used.append(f"rs4481887 (near OR2M7): {snp}")

    if not snp:
        return {
            "prediction": "Unknown",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "rs4481887 not available for asparagus smell prediction.",
        }

    # rs4481887: GG = likely can smell, GA = intermediate, AA = likely cannot smell
    g_count = snp.count("G")

    if g_count == 2:
        return {
            "prediction": "Likely can smell asparagus metabolites",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "GG genotype — associated with ability to detect the characteristic smell of asparagus in urine. Almost everyone produces the metabolites, but not everyone can smell them.",
        }
    elif g_count == 1:
        return {
            "prediction": "May partially detect asparagus smell",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "GA genotype — intermediate ability to detect asparagus metabolites in urine.",
        }
    else:
        return {
            "prediction": "Likely cannot smell asparagus metabolites",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "AA genotype — associated with asparagus anosmia (inability to smell asparagus metabolites). ~40% of people have this.",
        }


def _predict_muscle_fiber(genome_by_rsid):
    """Predict muscle fiber type from ACTN3 (rs1815739)."""
    actn3 = _get_genotype(genome_by_rsid, "rs1815739")

    snps_used = []
    if actn3:
        snps_used.append(f"rs1815739 (ACTN3 R577X): {actn3}")

    if not actn3:
        return {
            "prediction": "Unknown",
            "confidence": "low",
            "snps_used": snps_used,
            "description": "ACTN3 SNP (rs1815739) not available.",
        }

    # rs1815739: CC = RR (power), CT = RX (mixed), TT = XX (endurance)
    t_count = actn3.count("T")

    if t_count == 0:
        return {
            "prediction": "Power/sprint oriented (RR)",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "CC (RR) genotype — full alpha-actinin-3 expression in fast-twitch muscle fibers. Found in ~97% of Olympic sprinters. Favors explosive power sports.",
        }
    elif t_count == 1:
        return {
            "prediction": "Mixed power/endurance (RX)",
            "confidence": "moderate",
            "snps_used": snps_used,
            "description": "CT (RX) genotype — intermediate fast-twitch fiber composition. Versatile athletic profile suitable for both power and endurance activities.",
        }
    else:
        return {
            "prediction": "Endurance oriented (XX)",
            "confidence": "high",
            "snps_used": snps_used,
            "description": "TT (XX) genotype — no alpha-actinin-3 in fast-twitch fibers. Associated with endurance performance and better muscle recovery. ~18% of Europeans, ~25% of Asians.",
        }


def predict_traits(genome_by_rsid):
    """Predict visible traits from well-established SNP associations.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict mapping trait name to prediction dict with keys:
        prediction, confidence, snps_used, description.
    """
    return {
        "eye_color": _predict_eye_color(genome_by_rsid),
        "hair_color": _predict_hair_color(genome_by_rsid),
        "earwax_type": _predict_earwax(genome_by_rsid),
        "freckling": _predict_freckling(genome_by_rsid),
        "lactose_tolerance": _predict_lactose(genome_by_rsid),
        "bitter_taste": _predict_bitter_taste(genome_by_rsid),
        "cilantro_taste": _predict_cilantro(genome_by_rsid),
        "asparagus_smell": _predict_asparagus_smell(genome_by_rsid),
        "muscle_fiber_type": _predict_muscle_fiber(genome_by_rsid),
    }
