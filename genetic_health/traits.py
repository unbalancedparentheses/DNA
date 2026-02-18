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
    }
