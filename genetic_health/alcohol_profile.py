"""Unified alcohol metabolism profile from ADH1B, ALDH2, and CYP2E1 variants.

Alcohol metabolism follows a two-step pathway:
1. ADH1B converts ethanol to acetaldehyde (toxic)
2. ALDH2 converts acetaldehyde to acetate (harmless)

Variants in these enzymes affect metabolism speed, flush reaction,
and alcohol-related cancer risk.
"""


def _get_genotype(genome_by_rsid, rsid):
    """Get genotype string for an rsID, or empty string if missing."""
    entry = genome_by_rsid.get(rsid, {})
    return entry.get("genotype", "")


def _assess_adh1b(genome_by_rsid):
    """Assess alcohol dehydrogenase speed from ADH1B rs1229984 (Arg48His).

    The His48 variant (T allele) produces a hyperactive ADH1B enzyme
    that converts ethanol to acetaldehyde ~40x faster.
    """
    gt = _get_genotype(genome_by_rsid, "rs1229984")
    if not gt:
        return None

    t_count = gt.count("T")
    if t_count == 2:
        return {
            "rsid": "rs1229984",
            "gene": "ADH1B",
            "genotype": gt,
            "finding": "Very fast ethanol metabolism",
            "description": "TT (His/His) genotype -- hyperactive ADH1B produces acetaldehyde ~40x faster. Rapid acetaldehyde buildup causes unpleasant effects and is protective against alcoholism. Common in East Asian populations (~70%).",
            "speed_contribution": "very_fast",
            "cancer_modifier": 1,
        }
    elif t_count == 1:
        return {
            "rsid": "rs1229984",
            "gene": "ADH1B",
            "genotype": gt,
            "finding": "Fast ethanol metabolism",
            "description": "CT (Arg/His) genotype -- one hyperactive ADH1B allele. Faster-than-average ethanol-to-acetaldehyde conversion. Some protective effect against alcohol dependence.",
            "speed_contribution": "fast",
            "cancer_modifier": 0,
        }
    else:
        return {
            "rsid": "rs1229984",
            "gene": "ADH1B",
            "genotype": gt,
            "finding": "Normal ethanol metabolism",
            "description": "CC (Arg/Arg) genotype -- standard ADH1B activity. Normal rate of ethanol-to-acetaldehyde conversion. Most common in European populations (~90%).",
            "speed_contribution": "normal",
            "cancer_modifier": 0,
        }


def _assess_aldh2(genome_by_rsid):
    """Assess aldehyde dehydrogenase activity from ALDH2 rs671 (Glu504Lys).

    The Lys504 variant (A allele) causes a near-inactive ALDH2 enzyme,
    leading to acetaldehyde accumulation and the alcohol flush reaction.
    """
    gt = _get_genotype(genome_by_rsid, "rs671")
    if not gt:
        return None

    a_count = gt.count("A")
    if a_count == 2:
        return {
            "rsid": "rs671",
            "gene": "ALDH2",
            "genotype": gt,
            "finding": "Very reduced acetaldehyde clearance (severe flush)",
            "description": "AA (Lys/Lys) genotype -- near-complete ALDH2 deficiency. Unable to clear acetaldehyde effectively. Severe facial flushing, nausea, and tachycardia with even small amounts of alcohol. Dramatically elevated esophageal cancer risk if alcohol is consumed.",
            "flush_contribution": "severe",
            "cancer_modifier": 3,
        }
    elif a_count == 1:
        return {
            "rsid": "rs671",
            "gene": "ALDH2",
            "genotype": gt,
            "finding": "Reduced acetaldehyde clearance (mild flush)",
            "description": "GA (Glu/Lys) genotype -- ~60-70% reduced ALDH2 activity. Mild to moderate facial flushing with alcohol. Significantly elevated esophageal and head/neck cancer risk if drinking despite flush. Common in ~30-40% of East Asians.",
            "flush_contribution": "mild",
            "cancer_modifier": 2,
        }
    else:
        return {
            "rsid": "rs671",
            "gene": "ALDH2",
            "genotype": gt,
            "finding": "Normal acetaldehyde clearance (no flush)",
            "description": "GG (Glu/Glu) genotype -- fully functional ALDH2 enzyme. Efficient acetaldehyde clearance. No genetic predisposition to alcohol flush reaction.",
            "flush_contribution": "none",
            "cancer_modifier": 0,
        }


def _assess_cyp2e1(genome_by_rsid):
    """Assess CYP2E1 activity from rs2031920.

    CYP2E1 is an alternative ethanol oxidation pathway (MEOS) that
    becomes more active with chronic alcohol use.
    """
    gt = _get_genotype(genome_by_rsid, "rs2031920")
    if not gt:
        return None

    t_count = gt.count("T")
    if t_count == 2:
        return {
            "rsid": "rs2031920",
            "gene": "CYP2E1",
            "genotype": gt,
            "finding": "Reduced CYP2E1 metabolism",
            "description": "TT genotype -- reduced CYP2E1 transcription. Lower alternative ethanol oxidation pathway activity. May slightly reduce reactive oxygen species production from alcohol metabolism.",
            "speed_contribution": "reduced",
            "cancer_modifier": -1,
        }
    elif t_count == 1:
        return {
            "rsid": "rs2031920",
            "gene": "CYP2E1",
            "genotype": gt,
            "finding": "Reduced CYP2E1 metabolism",
            "description": "CT genotype -- mildly reduced CYP2E1 activity. Slight reduction in alternative ethanol oxidation.",
            "speed_contribution": "reduced",
            "cancer_modifier": 0,
        }
    else:
        return {
            "rsid": "rs2031920",
            "gene": "CYP2E1",
            "genotype": gt,
            "finding": "Normal CYP2E1 metabolism",
            "description": "CC genotype -- standard CYP2E1 expression and activity. Normal microsomal ethanol oxidizing system.",
            "speed_contribution": "normal",
            "cancer_modifier": 0,
        }


def profile_alcohol(genome_by_rsid):
    """Profile alcohol metabolism from ADH1B, ALDH2, and CYP2E1 variants.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict with keys:
        metabolism_speed : str ("slow", "normal", or "fast")
        flush_risk : str ("none", "mild", or "severe")
        cancer_risk : str ("average", "elevated", or "high")
        markers_found : int
        markers_tested : int
        gene_results : list of dicts
        summary : str
        recommendations : list of str
    """
    assessments = {
        "adh1b": _assess_adh1b(genome_by_rsid),
        "aldh2": _assess_aldh2(genome_by_rsid),
        "cyp2e1": _assess_cyp2e1(genome_by_rsid),
    }

    gene_results = [a for a in assessments.values() if a is not None]
    markers_found = len(gene_results)
    markers_tested = len(assessments)

    # Determine metabolism speed
    if assessments["adh1b"]:
        speed_val = assessments["adh1b"]["speed_contribution"]
        if speed_val == "very_fast":
            metabolism_speed = "fast"
        elif speed_val == "fast":
            metabolism_speed = "fast"
        else:
            metabolism_speed = "normal"
    else:
        metabolism_speed = "normal"

    # CYP2E1 can modify speed slightly
    if assessments["cyp2e1"] and assessments["cyp2e1"]["speed_contribution"] == "reduced":
        if metabolism_speed == "normal":
            metabolism_speed = "slow"

    # Determine flush risk
    if assessments["aldh2"]:
        flush_risk = assessments["aldh2"]["flush_contribution"]
    else:
        flush_risk = "unknown"

    # Determine cancer risk from alcohol
    cancer_score = sum(a.get("cancer_modifier", 0) for a in gene_results)
    if cancer_score >= 3:
        cancer_risk = "high"
    elif cancer_score >= 1:
        cancer_risk = "elevated"
    else:
        cancer_risk = "average"

    # Summary
    if markers_found == 0:
        summary = "No alcohol metabolism markers available in genotype data."
    else:
        parts = []
        if metabolism_speed == "fast":
            parts.append("fast ethanol-to-acetaldehyde conversion (ADH1B)")
        elif metabolism_speed == "slow":
            parts.append("slower ethanol metabolism")

        if flush_risk == "severe":
            parts.append("severe alcohol flush reaction (ALDH2 deficiency)")
        elif flush_risk == "mild":
            parts.append("mild alcohol flush risk (ALDH2 partial deficiency)")

        if cancer_risk == "high":
            parts.append("high alcohol-related cancer risk")
        elif cancer_risk == "elevated":
            parts.append("elevated alcohol-related cancer risk")

        if parts:
            summary = "Your alcohol metabolism profile shows: " + "; ".join(parts) + "."
        else:
            summary = "Your alcohol metabolism profile is typical. No notable genetic variants affecting alcohol processing detected."

    # Recommendations
    recommendations = []
    if flush_risk == "severe":
        recommendations.extend([
            "ALDH2 deficiency detected -- strongly consider avoiding alcohol entirely.",
            "Even small amounts of alcohol cause toxic acetaldehyde accumulation, increasing esophageal cancer risk up to 10x.",
            "Do not take medications that inhibit ALDH (e.g., disulfiram) without medical supervision.",
        ])
    elif flush_risk == "mild":
        recommendations.extend([
            "Partial ALDH2 deficiency detected -- limit alcohol consumption significantly.",
            "Drinking despite flush reaction substantially increases esophageal and head/neck cancer risk.",
            "If you choose to drink, limit to very small quantities and never on an empty stomach.",
        ])

    if metabolism_speed == "fast" and flush_risk in ("none", "unknown"):
        recommendations.append(
            "Fast ADH1B metabolism provides some natural protection against alcohol dependence, but does not eliminate health risks from heavy drinking."
        )

    if cancer_risk in ("elevated", "high"):
        recommendations.append(
            "Discuss alcohol-related cancer screening (esophageal, head/neck) with your physician, especially if you have a history of alcohol consumption."
        )

    if not recommendations and markers_found > 0:
        recommendations.append(
            "Standard alcohol guidelines apply. No genetic variants requiring special precautions detected."
        )

    return {
        "metabolism_speed": metabolism_speed,
        "flush_risk": flush_risk,
        "cancer_risk": cancer_risk,
        "markers_found": markers_found,
        "markers_tested": markers_tested,
        "gene_results": gene_results,
        "summary": summary,
        "recommendations": recommendations,
    }
