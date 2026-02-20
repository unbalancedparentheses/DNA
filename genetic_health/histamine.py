"""Profile histamine intolerance risk from DAO and HNMT enzyme genetics.

Histamine intolerance arises when degradation capacity is overwhelmed.
Two enzymes handle histamine clearance:
- DAO (diamine oxidase, encoded by AOC1) -- extracellular, gut-focused
- HNMT (histamine N-methyltransferase) -- intracellular

Variants reducing either enzyme can cause symptoms after histamine-rich foods.
"""


def _get_genotype(genome_by_rsid, rsid):
    """Get genotype string for an rsID, or empty string if missing."""
    entry = genome_by_rsid.get(rsid, {})
    return entry.get("genotype", "")


def _assess_dao_rs10156191(genome_by_rsid):
    """Assess DAO activity from AOC1 rs10156191."""
    gt = _get_genotype(genome_by_rsid, "rs10156191")
    if not gt:
        return None

    t_count = gt.count("T")
    if t_count == 2:
        return {
            "rsid": "rs10156191",
            "gene": "AOC1 (DAO)",
            "genotype": gt,
            "finding": "Low DAO activity",
            "description": "TT genotype -- substantially reduced diamine oxidase activity. Impaired histamine degradation in the gut.",
            "risk_points": 2,
        }
    elif t_count == 1:
        return {
            "rsid": "rs10156191",
            "gene": "AOC1 (DAO)",
            "genotype": gt,
            "finding": "Reduced DAO activity",
            "description": "CT genotype -- moderately reduced DAO enzyme activity. Partial impairment of gut histamine clearance.",
            "risk_points": 1,
        }
    else:
        return {
            "rsid": "rs10156191",
            "gene": "AOC1 (DAO)",
            "genotype": gt,
            "finding": "Normal DAO activity",
            "description": "CC genotype -- normal diamine oxidase expression and activity.",
            "risk_points": 0,
        }


def _assess_dao_rs1049793(genome_by_rsid):
    """Assess DAO activity from AOC1 rs1049793."""
    gt = _get_genotype(genome_by_rsid, "rs1049793")
    if not gt:
        return None

    c_count = gt.count("C")
    if c_count == 2:
        return {
            "rsid": "rs1049793",
            "gene": "AOC1 (DAO)",
            "genotype": gt,
            "finding": "Low DAO activity",
            "description": "CC genotype -- reduced DAO enzyme stability. Combined with other DAO variants, significantly impairs histamine clearance.",
            "risk_points": 2,
        }
    elif c_count == 1:
        return {
            "rsid": "rs1049793",
            "gene": "AOC1 (DAO)",
            "genotype": gt,
            "finding": "Reduced DAO activity",
            "description": "GC genotype -- mildly reduced DAO stability. Some impairment of histamine degradation capacity.",
            "risk_points": 1,
        }
    else:
        return {
            "rsid": "rs1049793",
            "gene": "AOC1 (DAO)",
            "genotype": gt,
            "finding": "Normal DAO activity",
            "description": "GG genotype -- normal DAO enzyme stability and function.",
            "risk_points": 0,
        }


def _assess_hnmt(genome_by_rsid):
    """Assess HNMT activity from rs11558538."""
    gt = _get_genotype(genome_by_rsid, "rs11558538")
    if not gt:
        return None

    t_count = gt.count("T")
    if t_count == 2:
        return {
            "rsid": "rs11558538",
            "gene": "HNMT",
            "genotype": gt,
            "finding": "Low HNMT activity",
            "description": "TT genotype -- substantially reduced histamine N-methyltransferase activity. Impaired intracellular histamine clearance.",
            "risk_points": 2,
        }
    elif t_count == 1:
        return {
            "rsid": "rs11558538",
            "gene": "HNMT",
            "genotype": gt,
            "finding": "Reduced HNMT activity",
            "description": "CT genotype -- moderately reduced HNMT activity. Some impairment of intracellular histamine degradation.",
            "risk_points": 1,
        }
    else:
        return {
            "rsid": "rs11558538",
            "gene": "HNMT",
            "genotype": gt,
            "finding": "Normal HNMT activity",
            "description": "CC genotype -- normal HNMT enzyme function. Efficient intracellular histamine clearance.",
            "risk_points": 0,
        }


def profile_histamine(genome_by_rsid):
    """Profile histamine intolerance risk from DAO and HNMT variants.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict with keys:
        risk_level : str ("low", "moderate", or "elevated")
        markers_found : int
        markers_tested : int
        gene_results : list of dicts
        summary : str
        recommendations : list of str
        foods_to_watch : list of str
    """
    assessments = [
        _assess_dao_rs10156191(genome_by_rsid),
        _assess_dao_rs1049793(genome_by_rsid),
        _assess_hnmt(genome_by_rsid),
    ]

    gene_results = [a for a in assessments if a is not None]
    markers_found = len(gene_results)
    markers_tested = len(assessments)

    # Total risk score (0-6)
    total_risk = sum(a["risk_points"] for a in gene_results)

    if markers_found == 0:
        risk_level = "unknown"
    elif total_risk >= 4:
        risk_level = "elevated"
    elif total_risk >= 2:
        risk_level = "moderate"
    else:
        risk_level = "low"

    # Summary
    if markers_found == 0:
        summary = "No histamine metabolism markers available in genotype data."
    elif risk_level == "elevated":
        summary = (
            "Your genetic profile suggests elevated histamine intolerance risk. "
            "Both DAO and/or HNMT pathways show reduced activity, which may lead to "
            "symptoms like headaches, flushing, nasal congestion, or GI issues after "
            "histamine-rich foods."
        )
    elif risk_level == "moderate":
        summary = (
            "Your genetic profile suggests moderate histamine intolerance risk. "
            "Some reduction in histamine-degrading enzyme activity detected. "
            "You may notice symptoms with large amounts of histamine-rich foods."
        )
    else:
        summary = (
            "Your genetic profile suggests low histamine intolerance risk. "
            "Histamine-degrading enzyme activity appears normal based on tested markers."
        )

    # Recommendations
    recommendations = []
    if risk_level == "elevated":
        recommendations.extend([
            "Consider a low-histamine diet trial (2-4 weeks) to assess symptom improvement.",
            "DAO supplementation before meals may help with food-triggered symptoms.",
            "Track symptoms with a food diary to identify personal triggers.",
            "Discuss with your physician -- symptoms overlap with allergies and mast cell disorders.",
        ])
    elif risk_level == "moderate":
        recommendations.extend([
            "Be aware of cumulative histamine load -- spacing high-histamine meals may help.",
            "Freshly prepared foods have lower histamine than leftovers or fermented foods.",
            "Consider DAO supplementation if you notice post-meal flushing or headaches.",
        ])
    else:
        recommendations.append(
            "No specific dietary modifications indicated based on histamine genetics alone."
        )

    # Foods to watch (relevant for moderate/elevated risk)
    foods_to_watch = [
        "Aged cheeses (parmesan, cheddar, gouda)",
        "Fermented foods (sauerkraut, kimchi, kombucha)",
        "Cured meats (salami, prosciutto, bacon)",
        "Alcoholic beverages (especially red wine and beer)",
        "Smoked fish and canned tuna",
        "Tomatoes, spinach, and eggplant",
        "Vinegar and soy sauce",
        "Leftover/reheated meals (histamine increases with time)",
    ]

    return {
        "risk_level": risk_level,
        "markers_found": markers_found,
        "markers_tested": markers_tested,
        "gene_results": gene_results,
        "summary": summary,
        "recommendations": recommendations,
        "foods_to_watch": foods_to_watch,
    }
