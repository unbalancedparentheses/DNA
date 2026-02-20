"""Profile eye health genetics beyond age-related macular degeneration (covered in PRS).

Evaluates variants affecting:
- Glaucoma risk: MYOC (rs74315329), CAV1 (rs4236601), LOXL1 (rs1048661)
- Myopia risk: GJD2 (rs524952), RASGRF1 (rs8027411)
"""


def _get_genotype(genome_by_rsid, rsid):
    """Get genotype string for an rsID, or empty string if missing."""
    entry = genome_by_rsid.get(rsid, {})
    return entry.get("genotype", "")


def _assess_myoc(genome_by_rsid):
    """Assess glaucoma risk from MYOC rs74315329.

    MYOC encodes myocilin. The Gln368Ter variant (rare) causes
    primary open-angle glaucoma (POAG) with high penetrance.
    """
    gt = _get_genotype(genome_by_rsid, "rs74315329")
    if not gt:
        return None

    # This is a rare pathogenic variant; A = stop-gain (Ter), G = normal (Gln)
    a_count = gt.count("A")
    if a_count >= 1:
        return {
            "rsid": "rs74315329",
            "gene": "MYOC",
            "genotype": gt,
            "finding": "Myocilin glaucoma variant detected",
            "description": f"{'Homozygous' if a_count == 2 else 'Heterozygous'} for MYOC Gln368Ter -- a high-penetrance glaucoma variant. Strongly associated with juvenile/early-onset primary open-angle glaucoma. Ophthalmology referral recommended.",
            "condition": "glaucoma_risk",
            "severity": 3,
        }
    else:
        return {
            "rsid": "rs74315329",
            "gene": "MYOC",
            "genotype": gt,
            "finding": "No MYOC glaucoma variant",
            "description": "GG genotype -- no myocilin Gln368Ter variant detected. This does not rule out other glaucoma risk factors.",
            "condition": "glaucoma_risk",
            "severity": 0,
        }


def _assess_cav1(genome_by_rsid):
    """Assess glaucoma risk from CAV1 rs4236601.

    CAV1 encodes caveolin-1, involved in aqueous humor outflow regulation.
    """
    gt = _get_genotype(genome_by_rsid, "rs4236601")
    if not gt:
        return None

    a_count = gt.count("A")
    if a_count == 2:
        return {
            "rsid": "rs4236601",
            "gene": "CAV1",
            "genotype": gt,
            "finding": "Higher glaucoma risk",
            "description": "AA genotype -- strongest association with primary open-angle glaucoma at the CAV1/CAV2 locus. Increased intraocular pressure risk.",
            "condition": "glaucoma_risk",
            "severity": 2,
        }
    elif a_count == 1:
        return {
            "rsid": "rs4236601",
            "gene": "CAV1",
            "genotype": gt,
            "finding": "Moderate glaucoma risk",
            "description": "AC genotype -- one risk allele for glaucoma at the CAV1 locus. Mildly increased susceptibility.",
            "condition": "glaucoma_risk",
            "severity": 1,
        }
    else:
        return {
            "rsid": "rs4236601",
            "gene": "CAV1",
            "genotype": gt,
            "finding": "Average glaucoma risk (CAV1)",
            "description": "CC genotype -- no CAV1 risk alleles. Average population risk for open-angle glaucoma at this locus.",
            "condition": "glaucoma_risk",
            "severity": 0,
        }


def _assess_loxl1(genome_by_rsid):
    """Assess exfoliation glaucoma risk from LOXL1 rs1048661 (R141L).

    LOXL1 variants are the strongest known genetic risk factor for
    exfoliation syndrome/glaucoma.
    """
    gt = _get_genotype(genome_by_rsid, "rs1048661")
    if not gt:
        return None

    g_count = gt.count("G")
    if g_count == 2:
        return {
            "rsid": "rs1048661",
            "gene": "LOXL1",
            "genotype": gt,
            "finding": "Higher exfoliation glaucoma risk",
            "description": "GG genotype -- homozygous for LOXL1 risk allele. Significantly increased risk for exfoliation syndrome and exfoliation glaucoma, especially after age 60.",
            "condition": "glaucoma_risk",
            "severity": 2,
        }
    elif g_count == 1:
        return {
            "rsid": "rs1048661",
            "gene": "LOXL1",
            "genotype": gt,
            "finding": "Moderate exfoliation glaucoma risk",
            "description": "GT genotype -- one LOXL1 risk allele. Moderately increased risk for exfoliation syndrome.",
            "condition": "glaucoma_risk",
            "severity": 1,
        }
    else:
        return {
            "rsid": "rs1048661",
            "gene": "LOXL1",
            "genotype": gt,
            "finding": "Lower exfoliation glaucoma risk",
            "description": "TT genotype -- no LOXL1 risk alleles at this position. Lower risk for exfoliation syndrome/glaucoma.",
            "condition": "glaucoma_risk",
            "severity": 0,
        }


def _assess_gjd2(genome_by_rsid):
    """Assess myopia risk from GJD2 rs524952.

    GJD2 encodes connexin-36, involved in retinal signaling and
    eye growth regulation.
    """
    gt = _get_genotype(genome_by_rsid, "rs524952")
    if not gt:
        return None

    a_count = gt.count("A")
    if a_count == 2:
        return {
            "rsid": "rs524952",
            "gene": "GJD2",
            "genotype": gt,
            "finding": "Higher myopia risk (GJD2)",
            "description": "AA genotype -- associated with increased myopia susceptibility. GJD2 affects retinal gap junction signaling involved in eye growth.",
            "condition": "myopia_risk",
            "severity": 2,
        }
    elif a_count == 1:
        return {
            "rsid": "rs524952",
            "gene": "GJD2",
            "genotype": gt,
            "finding": "Moderate myopia risk (GJD2)",
            "description": "AG genotype -- one risk allele for myopia at GJD2. Mildly increased susceptibility to near-sightedness.",
            "condition": "myopia_risk",
            "severity": 1,
        }
    else:
        return {
            "rsid": "rs524952",
            "gene": "GJD2",
            "genotype": gt,
            "finding": "Average myopia risk (GJD2)",
            "description": "GG genotype -- no GJD2 risk alleles. Average genetic risk for myopia at this locus.",
            "condition": "myopia_risk",
            "severity": 0,
        }


def _assess_rasgrf1(genome_by_rsid):
    """Assess myopia risk from RASGRF1 rs8027411.

    RASGRF1 is involved in retinal signaling and axial eye growth.
    """
    gt = _get_genotype(genome_by_rsid, "rs8027411")
    if not gt:
        return None

    t_count = gt.count("T")
    if t_count == 2:
        return {
            "rsid": "rs8027411",
            "gene": "RASGRF1",
            "genotype": gt,
            "finding": "Higher myopia risk (RASGRF1)",
            "description": "TT genotype -- associated with increased myopia susceptibility through altered retinal RAS/MAPK signaling and eye growth.",
            "condition": "myopia_risk",
            "severity": 2,
        }
    elif t_count == 1:
        return {
            "rsid": "rs8027411",
            "gene": "RASGRF1",
            "genotype": gt,
            "finding": "Moderate myopia risk (RASGRF1)",
            "description": "GT genotype -- one risk allele for myopia at RASGRF1. Mildly increased near-sightedness susceptibility.",
            "condition": "myopia_risk",
            "severity": 1,
        }
    else:
        return {
            "rsid": "rs8027411",
            "gene": "RASGRF1",
            "genotype": gt,
            "finding": "Average myopia risk (RASGRF1)",
            "description": "GG genotype -- no RASGRF1 risk alleles. Average genetic myopia risk at this locus.",
            "condition": "myopia_risk",
            "severity": 0,
        }


def profile_eye_health(genome_by_rsid):
    """Profile eye health genetics for glaucoma and myopia risk.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict with keys:
        conditions : dict with glaucoma_risk and myopia_risk sub-dicts
        markers_found : int
        markers_tested : int
        gene_results : list of dicts
        summary : str
        recommendations : list of str
    """
    assessments = [
        _assess_myoc(genome_by_rsid),
        _assess_cav1(genome_by_rsid),
        _assess_loxl1(genome_by_rsid),
        _assess_gjd2(genome_by_rsid),
        _assess_rasgrf1(genome_by_rsid),
    ]

    gene_results = [a for a in assessments if a is not None]
    markers_found = len(gene_results)
    markers_tested = len(assessments)

    # Group by condition
    condition_scores = {}
    condition_details = {}
    for result in gene_results:
        cond = result["condition"]
        if cond not in condition_scores:
            condition_scores[cond] = []
            condition_details[cond] = []
        condition_scores[cond].append(result["severity"])
        condition_details[cond].append({
            "gene": result["gene"],
            "genotype": result["genotype"],
            "finding": result["finding"],
        })

    def _classify_condition(scores):
        if not scores:
            return "unknown"
        max_score = max(scores)
        avg_score = sum(scores) / len(scores)
        # A single high-penetrance variant (severity 3) dominates
        if max_score >= 3:
            return "high"
        if avg_score >= 1.5:
            return "elevated"
        elif avg_score >= 0.8:
            return "moderate"
        return "average"

    conditions = {}
    for cond_name in ["glaucoma_risk", "myopia_risk"]:
        scores = condition_scores.get(cond_name, [])
        conditions[cond_name] = {
            "level": _classify_condition(scores),
            "markers": condition_details.get(cond_name, []),
        }

    # Summary
    if markers_found == 0:
        summary = "No eye health markers available in genotype data."
    else:
        findings = []
        glaucoma_level = conditions["glaucoma_risk"]["level"]
        myopia_level = conditions["myopia_risk"]["level"]

        if glaucoma_level == "high":
            findings.append("high-penetrance glaucoma variant detected (MYOC) -- urgent ophthalmology referral recommended")
        elif glaucoma_level == "elevated":
            findings.append("elevated glaucoma risk")
        elif glaucoma_level == "moderate":
            findings.append("moderate glaucoma risk")

        if myopia_level in ("elevated", "high"):
            findings.append("elevated myopia susceptibility")
        elif myopia_level == "moderate":
            findings.append("moderate myopia susceptibility")

        if findings:
            summary = "Your eye health profile shows: " + "; ".join(findings) + "."
        else:
            summary = "Your eye health genetic profile appears average across all tested markers."

    # Recommendations
    recommendations = []
    glaucoma_level = conditions["glaucoma_risk"]["level"]
    myopia_level = conditions["myopia_risk"]["level"]

    if glaucoma_level == "high":
        recommendations.extend([
            "Urgent: MYOC glaucoma variant detected. Schedule an ophthalmology evaluation for intraocular pressure and optic nerve assessment.",
            "Annual comprehensive eye exams with IOP measurement and visual field testing.",
            "First-degree relatives should also be screened for this variant.",
        ])
    elif glaucoma_level in ("elevated", "moderate"):
        recommendations.extend([
            "Regular comprehensive eye exams including intraocular pressure measurement, especially after age 40.",
            "Report any vision changes (halos, peripheral vision loss) to your ophthalmologist promptly.",
        ])

    if myopia_level in ("elevated", "high"):
        recommendations.extend([
            "For children: encourage outdoor time (at least 2 hours/day) -- the strongest modifiable factor for myopia prevention.",
            "Follow the 20-20-20 rule: every 20 minutes of near work, look at something 20 feet away for 20 seconds.",
            "Regular eye exams to monitor refractive changes and screen for myopia complications (retinal detachment, macular degeneration).",
        ])
    elif myopia_level == "moderate":
        recommendations.append(
            "Moderate myopia genetic risk -- regular vision screening and outdoor time recommended."
        )

    if not recommendations and markers_found > 0:
        recommendations.append(
            "No specific eye health interventions indicated based on tested genetic markers. Routine eye exams still recommended."
        )

    return {
        "conditions": conditions,
        "markers_found": markers_found,
        "markers_tested": markers_tested,
        "gene_results": gene_results,
        "summary": summary,
        "recommendations": recommendations,
    }
