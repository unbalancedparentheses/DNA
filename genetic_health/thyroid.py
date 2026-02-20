"""Profile thyroid genetics from SNPs affecting autoimmune risk, hormone conversion, and cancer risk.

Evaluates variants in:
- FOXE1 (rs965513): thyroid cancer susceptibility
- TPO (rs2071403): autoimmune thyroid disease (Hashimoto's, Graves')
- TSHR (rs179247): Graves' disease susceptibility
- DIO1 (rs11206244): type 1 deiodinase (T4 to T3 conversion)
- DIO2 (rs225014): type 2 deiodinase (T4 to T3 conversion, Thr92Ala)
"""


def _get_genotype(genome_by_rsid, rsid):
    """Get genotype string for an rsID, or empty string if missing."""
    entry = genome_by_rsid.get(rsid, {})
    return entry.get("genotype", "")


def _assess_foxe1(genome_by_rsid):
    """Assess thyroid cancer risk from FOXE1 rs965513."""
    gt = _get_genotype(genome_by_rsid, "rs965513")
    if not gt:
        return None

    a_count = gt.count("A")
    if a_count == 2:
        return {
            "rsid": "rs965513",
            "gene": "FOXE1",
            "genotype": gt,
            "finding": "Higher thyroid cancer risk",
            "description": "AA genotype -- strongest association with thyroid cancer susceptibility at the FOXE1/PTCSC2 locus. Regular thyroid monitoring recommended.",
            "domain": "cancer_risk",
            "severity": 2,
        }
    elif a_count == 1:
        return {
            "rsid": "rs965513",
            "gene": "FOXE1",
            "genotype": gt,
            "finding": "Moderate thyroid cancer risk",
            "description": "AG genotype -- one risk allele at the FOXE1 locus. Mildly increased thyroid cancer susceptibility.",
            "domain": "cancer_risk",
            "severity": 1,
        }
    else:
        return {
            "rsid": "rs965513",
            "gene": "FOXE1",
            "genotype": gt,
            "finding": "Average thyroid cancer risk",
            "description": "GG genotype -- no risk alleles at the FOXE1 locus. Average population risk for thyroid cancer.",
            "domain": "cancer_risk",
            "severity": 0,
        }


def _assess_tpo(genome_by_rsid):
    """Assess autoimmune thyroid risk from TPO rs2071403."""
    gt = _get_genotype(genome_by_rsid, "rs2071403")
    if not gt:
        return None

    c_count = gt.count("C")
    if c_count == 2:
        return {
            "rsid": "rs2071403",
            "gene": "TPO",
            "genotype": gt,
            "finding": "Higher autoimmune thyroid risk",
            "description": "CC genotype -- associated with increased anti-TPO antibody production. Higher susceptibility to Hashimoto's thyroiditis.",
            "domain": "autoimmune_risk",
            "severity": 2,
        }
    elif c_count == 1:
        return {
            "rsid": "rs2071403",
            "gene": "TPO",
            "genotype": gt,
            "finding": "Moderate autoimmune thyroid risk",
            "description": "CT genotype -- one risk allele for autoimmune thyroid disease. Moderate susceptibility.",
            "domain": "autoimmune_risk",
            "severity": 1,
        }
    else:
        return {
            "rsid": "rs2071403",
            "gene": "TPO",
            "genotype": gt,
            "finding": "Average autoimmune thyroid risk",
            "description": "TT genotype -- no risk alleles at this TPO locus. Average risk for autoimmune thyroid disease.",
            "domain": "autoimmune_risk",
            "severity": 0,
        }


def _assess_tshr(genome_by_rsid):
    """Assess Graves' disease susceptibility from TSHR rs179247."""
    gt = _get_genotype(genome_by_rsid, "rs179247")
    if not gt:
        return None

    a_count = gt.count("A")
    if a_count == 2:
        return {
            "rsid": "rs179247",
            "gene": "TSHR",
            "genotype": gt,
            "finding": "Higher Graves' disease susceptibility",
            "description": "AA genotype -- TSH receptor variant associated with increased Graves' disease risk through altered receptor expression.",
            "domain": "autoimmune_risk",
            "severity": 2,
        }
    elif a_count == 1:
        return {
            "rsid": "rs179247",
            "gene": "TSHR",
            "genotype": gt,
            "finding": "Moderate Graves' disease susceptibility",
            "description": "AG genotype -- one risk allele for Graves' disease. Mildly increased susceptibility.",
            "domain": "autoimmune_risk",
            "severity": 1,
        }
    else:
        return {
            "rsid": "rs179247",
            "gene": "TSHR",
            "genotype": gt,
            "finding": "Average Graves' disease susceptibility",
            "description": "GG genotype -- no TSHR risk alleles. Average population risk for Graves' disease.",
            "domain": "autoimmune_risk",
            "severity": 0,
        }


def _assess_dio1(genome_by_rsid):
    """Assess T4-to-T3 conversion from DIO1 rs11206244."""
    gt = _get_genotype(genome_by_rsid, "rs11206244")
    if not gt:
        return None

    t_count = gt.count("T")
    if t_count == 2:
        return {
            "rsid": "rs11206244",
            "gene": "DIO1",
            "genotype": gt,
            "finding": "Impaired T4-to-T3 conversion (DIO1)",
            "description": "TT genotype -- reduced type 1 deiodinase activity. May have higher T4 and lower T3 levels, potentially requiring combination T4/T3 therapy.",
            "domain": "conversion_efficiency",
            "severity": 2,
        }
    elif t_count == 1:
        return {
            "rsid": "rs11206244",
            "gene": "DIO1",
            "genotype": gt,
            "finding": "Reduced T4-to-T3 conversion (DIO1)",
            "description": "CT genotype -- mildly reduced DIO1 activity. Slight decrease in peripheral T4-to-T3 conversion.",
            "domain": "conversion_efficiency",
            "severity": 1,
        }
    else:
        return {
            "rsid": "rs11206244",
            "gene": "DIO1",
            "genotype": gt,
            "finding": "Normal T4-to-T3 conversion (DIO1)",
            "description": "CC genotype -- normal type 1 deiodinase activity. Efficient peripheral T4-to-T3 conversion.",
            "domain": "conversion_efficiency",
            "severity": 0,
        }


def _assess_dio2(genome_by_rsid):
    """Assess T4-to-T3 conversion from DIO2 rs225014 (Thr92Ala)."""
    gt = _get_genotype(genome_by_rsid, "rs225014")
    if not gt:
        return None

    t_count = gt.count("T")
    if t_count == 2:
        return {
            "rsid": "rs225014",
            "gene": "DIO2",
            "genotype": gt,
            "finding": "Slower T4-to-T3 conversion (DIO2 Ala/Ala)",
            "description": "TT (Ala/Ala) genotype -- reduced type 2 deiodinase activity. May impair local T3 production in brain and other tissues. Some patients report better well-being on combination T4/T3.",
            "domain": "conversion_efficiency",
            "severity": 2,
        }
    elif t_count == 1:
        return {
            "rsid": "rs225014",
            "gene": "DIO2",
            "genotype": gt,
            "finding": "Reduced T4-to-T3 conversion (DIO2 Thr/Ala)",
            "description": "CT (Thr/Ala) genotype -- mildly reduced DIO2 activity. Minor impact on local T3 availability.",
            "domain": "conversion_efficiency",
            "severity": 1,
        }
    else:
        return {
            "rsid": "rs225014",
            "gene": "DIO2",
            "genotype": gt,
            "finding": "Normal T4-to-T3 conversion (DIO2 Thr/Thr)",
            "description": "CC (Thr/Thr) genotype -- normal type 2 deiodinase activity. Efficient local T4-to-T3 conversion.",
            "domain": "conversion_efficiency",
            "severity": 0,
        }


def profile_thyroid(genome_by_rsid):
    """Profile thyroid genetics from autoimmune, conversion, and cancer risk SNPs.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict with keys:
        risk_profile : dict with autoimmune_risk, conversion_efficiency, cancer_risk
        markers_found : int
        markers_tested : int
        gene_results : list of dicts
        summary : str
        recommendations : list of str
    """
    assessments = [
        _assess_foxe1(genome_by_rsid),
        _assess_tpo(genome_by_rsid),
        _assess_tshr(genome_by_rsid),
        _assess_dio1(genome_by_rsid),
        _assess_dio2(genome_by_rsid),
    ]

    gene_results = [a for a in assessments if a is not None]
    markers_found = len(gene_results)
    markers_tested = len(assessments)

    # Build risk profile by domain
    domain_scores = {}
    domain_details = {}
    for result in gene_results:
        domain = result["domain"]
        if domain not in domain_scores:
            domain_scores[domain] = []
            domain_details[domain] = []
        domain_scores[domain].append(result["severity"])
        domain_details[domain].append({
            "gene": result["gene"],
            "genotype": result["genotype"],
            "finding": result["finding"],
        })

    def _classify_domain(scores):
        if not scores:
            return "unknown"
        avg = sum(scores) / len(scores)
        if avg >= 1.5:
            return "elevated"
        elif avg >= 0.8:
            return "moderate"
        return "average"

    risk_profile = {}
    for domain_name in ["autoimmune_risk", "conversion_efficiency", "cancer_risk"]:
        scores = domain_scores.get(domain_name, [])
        risk_profile[domain_name] = {
            "level": _classify_domain(scores),
            "markers": domain_details.get(domain_name, []),
        }

    # Summary
    summaries = []
    if markers_found == 0:
        summary = "No thyroid-related markers available in genotype data."
    else:
        ai_level = risk_profile["autoimmune_risk"]["level"]
        conv_level = risk_profile["conversion_efficiency"]["level"]
        cancer_level = risk_profile["cancer_risk"]["level"]

        if ai_level == "elevated":
            summaries.append("elevated autoimmune thyroid risk (Hashimoto's/Graves')")
        elif ai_level == "moderate":
            summaries.append("moderate autoimmune thyroid risk")

        if conv_level == "elevated":
            summaries.append("impaired T4-to-T3 conversion in both DIO1 and DIO2")
        elif conv_level == "moderate":
            summaries.append("mildly reduced T4-to-T3 conversion")

        if cancer_level == "elevated":
            summaries.append("elevated thyroid cancer susceptibility")
        elif cancer_level == "moderate":
            summaries.append("moderate thyroid cancer susceptibility")

        if summaries:
            summary = "Your thyroid genetic profile shows: " + "; ".join(summaries) + "."
        else:
            summary = "Your thyroid genetic profile appears average across all tested markers."

    # Recommendations
    recommendations = []
    if risk_profile["autoimmune_risk"]["level"] in ("elevated", "moderate"):
        recommendations.extend([
            "Monitor thyroid function annually (TSH, free T4, anti-TPO antibodies).",
            "Ensure adequate selenium intake (200 mcg/day) -- supports both deiodinase function and may reduce anti-TPO antibodies.",
            "Maintain adequate iodine intake but avoid excess supplementation, which can trigger autoimmune thyroid disease.",
        ])

    if risk_profile["conversion_efficiency"]["level"] in ("elevated", "moderate"):
        recommendations.extend([
            "If hypothyroid, discuss combination T4/T3 therapy with your endocrinologist -- DIO variants may impair T4-only response.",
            "Selenium and zinc support deiodinase enzyme function.",
            "Monitor both free T4 and free T3 levels, not just TSH.",
        ])

    if risk_profile["cancer_risk"]["level"] in ("elevated", "moderate"):
        recommendations.append(
            "Discuss thyroid ultrasound screening with your physician, especially if there is a family history of thyroid cancer."
        )

    if not recommendations and markers_found > 0:
        recommendations.append(
            "No specific thyroid interventions indicated based on genetic markers tested."
        )

    return {
        "risk_profile": risk_profile,
        "markers_found": markers_found,
        "markers_tested": markers_tested,
        "gene_results": gene_results,
        "summary": summary,
        "recommendations": recommendations,
    }
