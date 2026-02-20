"""Profile hormone metabolism genetics from aromatase, estrogen receptor, and androgen pathway SNPs.

Evaluates variants in:
- CYP19A1 rs4646 (aromatase): estrogen synthesis rate
- ESR1 rs2234693 (estrogen receptor alpha): estrogen sensitivity
- SRD5A2 rs523349 (5-alpha reductase type 2): DHT conversion
- AR rs6152 (androgen receptor): androgen sensitivity proxy
"""


def _get_genotype(genome_by_rsid, rsid):
    """Get genotype string for an rsID, or empty string if missing."""
    entry = genome_by_rsid.get(rsid, {})
    return entry.get("genotype", "")


def _assess_cyp19a1(genome_by_rsid):
    """Assess aromatase activity from CYP19A1 rs4646.

    CYP19A1 encodes aromatase, which converts androgens to estrogens.
    """
    gt = _get_genotype(genome_by_rsid, "rs4646")
    if not gt:
        return None

    t_count = gt.count("T")
    if t_count == 2:
        return {
            "rsid": "rs4646",
            "gene": "CYP19A1",
            "genotype": gt,
            "finding": "Higher aromatase activity",
            "description": "TT genotype -- associated with higher aromatase expression and increased estrogen production from androgen precursors.",
            "domain": "estrogen_metabolism",
            "estrogen_effect": "increased",
        }
    elif t_count == 1:
        return {
            "rsid": "rs4646",
            "gene": "CYP19A1",
            "genotype": gt,
            "finding": "Intermediate aromatase activity",
            "description": "CT genotype -- intermediate aromatase activity. Normal estrogen synthesis from androgen precursors.",
            "domain": "estrogen_metabolism",
            "estrogen_effect": "normal",
        }
    else:
        return {
            "rsid": "rs4646",
            "gene": "CYP19A1",
            "genotype": gt,
            "finding": "Lower aromatase activity",
            "description": "CC genotype -- associated with reduced aromatase activity. Lower rate of androgen-to-estrogen conversion.",
            "domain": "estrogen_metabolism",
            "estrogen_effect": "decreased",
        }


def _assess_esr1(genome_by_rsid):
    """Assess estrogen receptor sensitivity from ESR1 rs2234693 (PvuII).

    The PvuII polymorphism affects ESR1 transcription and estrogen signaling.
    """
    gt = _get_genotype(genome_by_rsid, "rs2234693")
    if not gt:
        return None

    c_count = gt.count("C")
    if c_count == 2:
        return {
            "rsid": "rs2234693",
            "gene": "ESR1",
            "genotype": gt,
            "finding": "Higher estrogen receptor sensitivity",
            "description": "CC genotype -- associated with increased ESR1 expression and enhanced estrogen signaling. May influence bone density, cardiovascular protection, and breast tissue sensitivity.",
            "domain": "estrogen_metabolism",
            "estrogen_effect": "enhanced_signaling",
        }
    elif c_count == 1:
        return {
            "rsid": "rs2234693",
            "gene": "ESR1",
            "genotype": gt,
            "finding": "Intermediate estrogen receptor sensitivity",
            "description": "CT genotype -- intermediate estrogen receptor expression. Typical estrogen signaling.",
            "domain": "estrogen_metabolism",
            "estrogen_effect": "normal",
        }
    else:
        return {
            "rsid": "rs2234693",
            "gene": "ESR1",
            "genotype": gt,
            "finding": "Lower estrogen receptor sensitivity",
            "description": "TT genotype -- lower ESR1 expression. May reduce estrogen-mediated effects on bone density and cardiovascular health.",
            "domain": "estrogen_metabolism",
            "estrogen_effect": "reduced_signaling",
        }


def _assess_srd5a2(genome_by_rsid):
    """Assess 5-alpha reductase activity from SRD5A2 rs523349 (V89L).

    SRD5A2 converts testosterone to dihydrotestosterone (DHT).
    """
    gt = _get_genotype(genome_by_rsid, "rs523349")
    if not gt:
        return None

    g_count = gt.count("G")
    if g_count == 2:
        return {
            "rsid": "rs523349",
            "gene": "SRD5A2",
            "genotype": gt,
            "finding": "Higher DHT conversion",
            "description": "GG (Val/Val) genotype -- normal/higher 5-alpha reductase activity. Higher DHT production from testosterone. Associated with male pattern baldness and prostate growth.",
            "domain": "androgen_metabolism",
            "androgen_effect": "increased_dht",
        }
    elif g_count == 1:
        return {
            "rsid": "rs523349",
            "gene": "SRD5A2",
            "genotype": gt,
            "finding": "Intermediate DHT conversion",
            "description": "GC (Val/Leu) genotype -- intermediate 5-alpha reductase activity. Moderate DHT production.",
            "domain": "androgen_metabolism",
            "androgen_effect": "normal",
        }
    else:
        return {
            "rsid": "rs523349",
            "gene": "SRD5A2",
            "genotype": gt,
            "finding": "Lower DHT conversion",
            "description": "CC (Leu/Leu) genotype -- reduced 5-alpha reductase activity (~40% lower). Lower DHT production may reduce androgenetic alopecia and benign prostate hyperplasia risk.",
            "domain": "androgen_metabolism",
            "androgen_effect": "decreased_dht",
        }


def _assess_ar(genome_by_rsid):
    """Assess androgen receptor sensitivity from AR rs6152 (CAG repeat proxy).

    This SNP serves as a tag for the CAG repeat length polymorphism
    in the androgen receptor.
    """
    gt = _get_genotype(genome_by_rsid, "rs6152")
    if not gt:
        return None

    g_count = gt.count("G")
    if g_count == 2:
        return {
            "rsid": "rs6152",
            "gene": "AR",
            "genotype": gt,
            "finding": "Standard androgen receptor sensitivity",
            "description": "GG genotype -- associated with typical CAG repeat length and normal androgen receptor sensitivity.",
            "domain": "androgen_metabolism",
            "androgen_effect": "normal_sensitivity",
        }
    elif g_count == 1:
        return {
            "rsid": "rs6152",
            "gene": "AR",
            "genotype": gt,
            "finding": "Intermediate androgen receptor sensitivity",
            "description": "AG genotype -- intermediate androgen receptor signaling based on CAG repeat proxy.",
            "domain": "androgen_metabolism",
            "androgen_effect": "intermediate_sensitivity",
        }
    else:
        return {
            "rsid": "rs6152",
            "gene": "AR",
            "genotype": gt,
            "finding": "Potentially altered androgen receptor sensitivity",
            "description": "AA genotype -- may tag longer CAG repeats associated with reduced androgen receptor transactivation. Lower androgen sensitivity despite normal hormone levels.",
            "domain": "androgen_metabolism",
            "androgen_effect": "reduced_sensitivity",
        }


def profile_hormone_metabolism(genome_by_rsid):
    """Profile hormone metabolism genetics from aromatase and androgen pathway SNPs.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict with keys:
        domains : dict with estrogen_metabolism and androgen_metabolism sub-dicts
        markers_found : int
        markers_tested : int
        gene_results : list of dicts
        summary : str
        recommendations : list of str
    """
    assessments = {
        "cyp19a1": _assess_cyp19a1(genome_by_rsid),
        "esr1": _assess_esr1(genome_by_rsid),
        "srd5a2": _assess_srd5a2(genome_by_rsid),
        "ar": _assess_ar(genome_by_rsid),
    }

    gene_results = [a for a in assessments.values() if a is not None]
    markers_found = len(gene_results)
    markers_tested = len(assessments)

    # Build domain summaries
    estrogen_markers = [a for a in gene_results if a["domain"] == "estrogen_metabolism"]
    androgen_markers = [a for a in gene_results if a["domain"] == "androgen_metabolism"]

    def _summarize_estrogen(markers):
        if not markers:
            return {"level": "unknown", "markers": []}
        findings = []
        for m in markers:
            findings.append({
                "gene": m["gene"],
                "genotype": m["genotype"],
                "finding": m["finding"],
            })
        effects = [m.get("estrogen_effect", "normal") for m in markers]
        if "increased" in effects and "enhanced_signaling" in effects:
            level = "high estrogen activity"
        elif "increased" in effects or "enhanced_signaling" in effects:
            level = "moderately increased estrogen activity"
        elif "decreased" in effects or "reduced_signaling" in effects:
            level = "reduced estrogen activity"
        else:
            level = "normal estrogen activity"
        return {"level": level, "markers": findings}

    def _summarize_androgen(markers):
        if not markers:
            return {"level": "unknown", "markers": []}
        findings = []
        for m in markers:
            findings.append({
                "gene": m["gene"],
                "genotype": m["genotype"],
                "finding": m["finding"],
            })
        effects = [m.get("androgen_effect", "normal") for m in markers]
        if "increased_dht" in effects:
            level = "elevated DHT pathway"
        elif "decreased_dht" in effects:
            level = "reduced DHT pathway"
        else:
            level = "normal androgen metabolism"
        return {"level": level, "markers": findings}

    domains = {
        "estrogen_metabolism": _summarize_estrogen(estrogen_markers),
        "androgen_metabolism": _summarize_androgen(androgen_markers),
    }

    # Determine overall profile
    e_level = domains["estrogen_metabolism"]["level"]
    a_level = domains["androgen_metabolism"]["level"]

    if markers_found == 0:
        overall = "unknown"
    else:
        notable = []
        if "increased" in e_level or "high" in e_level:
            notable.append("higher estrogen activity")
        elif "reduced" in e_level:
            notable.append("lower estrogen activity")
        if "elevated" in a_level:
            notable.append("elevated DHT conversion")
        elif "reduced" in a_level:
            notable.append("reduced DHT conversion")
        overall = ", ".join(notable) if notable else "balanced hormone metabolism"

    domains["overall"] = overall

    # Summary
    if markers_found == 0:
        summary = "No hormone metabolism markers available in genotype data."
    else:
        summary = f"Your hormone metabolism profile indicates {overall}."
        if overall == "balanced hormone metabolism":
            summary += " No notable deviations detected in the tested markers."

    # Recommendations
    recommendations = []
    if "increased" in e_level or "high" in e_level:
        recommendations.extend([
            "Higher aromatase activity may increase estrogen levels. Cruciferous vegetables (broccoli, cauliflower) support healthy estrogen metabolism via DIM/I3C.",
            "Maintain a healthy body fat percentage -- adipose tissue is a significant source of aromatase activity.",
            "Discuss estrogen-related cancer screening frequency with your physician.",
        ])
    if "reduced" in e_level:
        recommendations.extend([
            "Lower estrogen activity may affect bone density. Ensure adequate calcium, vitamin D, and weight-bearing exercise.",
            "Consider bone density screening, especially post-menopause or if other osteoporosis risk factors are present.",
        ])
    if "elevated" in a_level:
        recommendations.extend([
            "Higher DHT conversion may increase androgenetic alopecia and prostate growth risk.",
            "Saw palmetto and green tea (EGCG) are natural 5-alpha reductase inhibitors -- discuss with your physician.",
        ])
    if not recommendations and markers_found > 0:
        recommendations.append(
            "No specific hormone interventions indicated based on tested genetic markers."
        )

    return {
        "domains": domains,
        "markers_found": markers_found,
        "markers_tested": markers_tested,
        "gene_results": gene_results,
        "summary": summary,
        "recommendations": recommendations,
    }
