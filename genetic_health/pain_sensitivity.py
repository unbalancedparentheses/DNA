"""Profile pain sensitivity from well-characterized SNP associations.

Evaluates opioid response (OPRM1), pain threshold (COMT warrior/worrier),
sodium channel sensitivity (SCN9A), and capsaicin sensitivity (TRPV1)
to produce an overall pain sensitivity score.
"""


def _get_genotype(genome_by_rsid, rsid):
    """Get genotype string for an rsID, or empty string if missing."""
    entry = genome_by_rsid.get(rsid, {})
    return entry.get("genotype", "")


def _assess_oprm1(genome_by_rsid):
    """Assess opioid response from OPRM1 rs1799971 (A118G).

    A>G substitution reduces mu-opioid receptor expression and
    alters opioid binding affinity.
    """
    gt = _get_genotype(genome_by_rsid, "rs1799971")
    if not gt:
        return None

    g_count = gt.count("G")
    if g_count == 0:
        return {
            "rsid": "rs1799971",
            "gene": "OPRM1",
            "genotype": gt,
            "label": "Normal opioid response",
            "description": "AA genotype -- normal mu-opioid receptor expression. Standard response to opioid medications.",
            "sensitivity_contribution": 50,
        }
    elif g_count == 1:
        return {
            "rsid": "rs1799971",
            "gene": "OPRM1",
            "genotype": gt,
            "label": "Reduced opioid response",
            "description": "AG genotype -- reduced mu-opioid receptor binding. May require higher opioid doses for equivalent analgesia. Also associated with reduced placebo response.",
            "sensitivity_contribution": 35,
        }
    else:
        return {
            "rsid": "rs1799971",
            "gene": "OPRM1",
            "genotype": gt,
            "label": "Significantly reduced opioid response",
            "description": "GG genotype -- substantially reduced receptor expression. Markedly lower opioid efficacy; alternative pain management strategies may be needed.",
            "sensitivity_contribution": 20,
        }


def _assess_comt(genome_by_rsid):
    """Assess pain threshold from COMT rs4680 (Val158Met).

    Val/Val (GG) = rapid catecholamine clearance = lower pain sensitivity (warrior).
    Met/Met (AA) = slow clearance = higher pain sensitivity (worrier).
    """
    gt = _get_genotype(genome_by_rsid, "rs4680")
    if not gt:
        return None

    a_count = gt.count("A")
    if a_count == 0:
        return {
            "rsid": "rs4680",
            "gene": "COMT",
            "genotype": gt,
            "label": "Low pain sensitivity (Val/Val warrior)",
            "description": "GG genotype -- high COMT activity rapidly clears catecholamines. Lower pain sensitivity, better stress resilience, but lower dopamine availability at baseline.",
            "sensitivity_contribution": 25,
        }
    elif a_count == 1:
        return {
            "rsid": "rs4680",
            "gene": "COMT",
            "genotype": gt,
            "label": "Intermediate pain sensitivity (Val/Met)",
            "description": "AG genotype -- intermediate COMT enzyme activity. Balanced pain sensitivity and stress tolerance.",
            "sensitivity_contribution": 50,
        }
    else:
        return {
            "rsid": "rs4680",
            "gene": "COMT",
            "genotype": gt,
            "label": "High pain sensitivity (Met/Met worrier)",
            "description": "AA genotype -- low COMT activity, slower catecholamine clearance. Higher pain sensitivity but also higher cognitive performance under low-stress conditions.",
            "sensitivity_contribution": 80,
        }


def _assess_scn9a(genome_by_rsid):
    """Assess pain sensitivity from SCN9A rs6746030.

    SCN9A encodes the Nav1.7 sodium channel critical for pain signaling.
    """
    gt = _get_genotype(genome_by_rsid, "rs6746030")
    if not gt:
        return None

    a_count = gt.count("A")
    if a_count == 2:
        return {
            "rsid": "rs6746030",
            "gene": "SCN9A",
            "genotype": gt,
            "label": "Higher pain sensitivity",
            "description": "AA genotype -- associated with increased pain perception through enhanced Nav1.7 sodium channel activity.",
            "sensitivity_contribution": 80,
        }
    elif a_count == 1:
        return {
            "rsid": "rs6746030",
            "gene": "SCN9A",
            "genotype": gt,
            "label": "Intermediate pain sensitivity",
            "description": "AG genotype -- one copy of the high-sensitivity allele. Mildly increased pain perception.",
            "sensitivity_contribution": 60,
        }
    else:
        return {
            "rsid": "rs6746030",
            "gene": "SCN9A",
            "genotype": gt,
            "label": "Normal pain sensitivity",
            "description": "GG genotype -- typical Nav1.7 channel function. Standard pain signaling.",
            "sensitivity_contribution": 45,
        }


def _assess_trpv1(genome_by_rsid):
    """Assess capsaicin sensitivity from TRPV1 rs8065080.

    TRPV1 is the capsaicin receptor; variants affect heat/spice sensitivity.
    """
    gt = _get_genotype(genome_by_rsid, "rs8065080")
    if not gt:
        return None

    t_count = gt.count("T")
    if t_count == 2:
        return {
            "rsid": "rs8065080",
            "gene": "TRPV1",
            "genotype": gt,
            "label": "Higher capsaicin sensitivity",
            "description": "TT genotype -- altered TRPV1 receptor. Increased sensitivity to capsaicin, heat, and spicy foods.",
            "sensitivity_contribution": 80,
        }
    elif t_count == 1:
        return {
            "rsid": "rs8065080",
            "gene": "TRPV1",
            "genotype": gt,
            "label": "Intermediate capsaicin sensitivity",
            "description": "CT genotype -- one variant allele. Moderately increased capsaicin and heat sensitivity.",
            "sensitivity_contribution": 55,
        }
    else:
        return {
            "rsid": "rs8065080",
            "gene": "TRPV1",
            "genotype": gt,
            "label": "Normal capsaicin sensitivity",
            "description": "CC genotype -- typical TRPV1 receptor function. Standard response to capsaicin and heat stimuli.",
            "sensitivity_contribution": 40,
        }


def profile_pain_sensitivity(genome_by_rsid):
    """Profile pain sensitivity from multiple SNP domains.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict with keys:
        pain_sensitivity_score : int (0-100, higher = more sensitive)
        domains : dict with opioid_response, pain_threshold, capsaicin_sensitivity
        markers_found : int
        markers_tested : int
        summary : str
        recommendations : list of str
    """
    assessments = {
        "oprm1": _assess_oprm1(genome_by_rsid),
        "comt": _assess_comt(genome_by_rsid),
        "scn9a": _assess_scn9a(genome_by_rsid),
        "trpv1": _assess_trpv1(genome_by_rsid),
    }

    markers_found = sum(1 for v in assessments.values() if v is not None)
    markers_tested = len(assessments)

    # Build domains
    domains = {}

    if assessments["oprm1"]:
        domains["opioid_response"] = {
            "gene": "OPRM1",
            "genotype": assessments["oprm1"]["genotype"],
            "finding": assessments["oprm1"]["label"],
            "description": assessments["oprm1"]["description"],
        }

    if assessments["comt"]:
        domains["pain_threshold"] = {
            "gene": "COMT",
            "genotype": assessments["comt"]["genotype"],
            "finding": assessments["comt"]["label"],
            "description": assessments["comt"]["description"],
        }

    if assessments["scn9a"]:
        domains["pain_signaling"] = {
            "gene": "SCN9A",
            "genotype": assessments["scn9a"]["genotype"],
            "finding": assessments["scn9a"]["label"],
            "description": assessments["scn9a"]["description"],
        }

    if assessments["trpv1"]:
        domains["capsaicin_sensitivity"] = {
            "gene": "TRPV1",
            "genotype": assessments["trpv1"]["genotype"],
            "finding": assessments["trpv1"]["label"],
            "description": assessments["trpv1"]["description"],
        }

    # Compute overall score (0-100)
    contributions = [a["sensitivity_contribution"] for a in assessments.values() if a]
    if contributions:
        pain_sensitivity_score = round(sum(contributions) / len(contributions))
    else:
        pain_sensitivity_score = 50  # neutral default

    # Generate summary
    if markers_found == 0:
        summary = "No pain sensitivity markers available in genotype data."
    elif pain_sensitivity_score >= 70:
        summary = "Your genetic profile suggests higher-than-average pain sensitivity. You may benefit from proactive pain management strategies and should discuss opioid dosing with your physician."
    elif pain_sensitivity_score >= 50:
        summary = "Your genetic profile suggests average pain sensitivity across the markers tested."
    elif pain_sensitivity_score >= 30:
        summary = "Your genetic profile suggests lower-than-average pain sensitivity. You may have higher tolerance to pain stimuli but could need adjusted opioid dosing."
    else:
        summary = "Your genetic profile suggests notably low pain sensitivity. Opioid medications may be less effective; discuss alternative analgesics with your physician."

    # Recommendations
    recommendations = []
    if markers_found > 0:
        if assessments["oprm1"] and assessments["oprm1"]["sensitivity_contribution"] <= 35:
            recommendations.append(
                "Inform your anesthesiologist about OPRM1 AG/GG status -- you may require higher opioid doses or alternative analgesics."
            )
        if assessments["comt"] and assessments["comt"]["sensitivity_contribution"] >= 70:
            recommendations.append(
                "COMT Met/Met carriers may benefit from stress-reduction techniques (meditation, yoga) to manage heightened pain perception."
            )
        if assessments["trpv1"] and assessments["trpv1"]["sensitivity_contribution"] >= 70:
            recommendations.append(
                "Higher capsaicin sensitivity detected -- consider gradual spice exposure and topical capsaicin patches at lower concentrations."
            )
        recommendations.append(
            "Share your pain sensitivity profile with your healthcare provider for personalized pain management planning."
        )

    return {
        "pain_sensitivity_score": pain_sensitivity_score,
        "domains": domains,
        "markers_found": markers_found,
        "markers_tested": markers_tested,
        "summary": summary,
        "recommendations": recommendations,
    }
