"""Nutrigenomics profiling — micronutrient needs based on genotype.

Synthesizes vitamin/mineral metabolism variants into a personalized
nutrient profile showing likely deficiencies and conversion efficiencies.
"""


# Each nutrient has associated genes, their relevant SNPs/statuses,
# and the impact on that nutrient's metabolism.
NUTRIENT_PROFILES = {
    "folate_methylfolate": {
        "name": "Folate / Methylfolate (B9)",
        "genes": {
            "MTHFR": {
                "risk_statuses": ["reduced", "severely_reduced"],
                "impact": "Reduced conversion of folic acid to active 5-MTHF",
                "severity_map": {"reduced": 0.5, "severely_reduced": 1.0},
            },
            "MTRR": {
                "risk_statuses": ["reduced", "severely_reduced"],
                "impact": "Impaired B12 recycling reduces methionine synthase activity",
                "severity_map": {"reduced": 0.3, "severely_reduced": 0.7},
            },
            "MTR": {
                "risk_statuses": ["reduced"],
                "impact": "Reduced methionine synthase activity",
                "severity_map": {"reduced": 0.4},
            },
        },
        "food_sources": "Leafy greens, legumes, liver, fortified grains",
        "supplement_form": "L-methylfolate (5-MTHF) — NOT folic acid",
        "dose_range": "400-800 mcg methylfolate/day",
        "testing": "Serum folate, homocysteine, RBC folate",
    },
    "vitamin_b12": {
        "name": "Vitamin B12 (Cobalamin)",
        "genes": {
            "MTRR": {
                "risk_statuses": ["reduced", "severely_reduced"],
                "impact": "Impaired B12 recycling — may need higher intake",
                "severity_map": {"reduced": 0.3, "severely_reduced": 0.6},
            },
            "FUT2": {
                "risk_statuses": ["non_secretor"],
                "impact": "Non-secretors have ~35% lower B12 absorption from food",
                "severity_map": {"non_secretor": 0.5},
            },
        },
        "food_sources": "Meat, fish, eggs, dairy, fortified nutritional yeast",
        "supplement_form": "Methylcobalamin or adenosylcobalamin (sublingual)",
        "dose_range": "1000-5000 mcg methylcobalamin/day if deficient",
        "testing": "Serum B12, methylmalonic acid (MMA), homocysteine",
    },
    "vitamin_d": {
        "name": "Vitamin D",
        "genes": {
            "VDR": {
                "risk_statuses": ["reduced", "low_sensitivity"],
                "impact": "Reduced vitamin D receptor sensitivity",
                "severity_map": {"reduced": 0.4, "low_sensitivity": 0.5},
            },
            "GC": {
                "risk_statuses": ["reduced", "low_binding"],
                "impact": "Lower vitamin D binding protein — reduced transport",
                "severity_map": {"reduced": 0.5, "low_binding": 0.6},
            },
            "CYP2R1": {
                "risk_statuses": ["reduced"],
                "impact": "Reduced 25-hydroxylation of vitamin D in liver",
                "severity_map": {"reduced": 0.5},
            },
        },
        "food_sources": "Fatty fish, egg yolks, fortified foods, sunlight exposure",
        "supplement_form": "Vitamin D3 (cholecalciferol) with K2 (MK-7)",
        "dose_range": "2000-5000 IU D3/day (adjust by serum levels)",
        "testing": "Serum 25-OH vitamin D (target 40-60 ng/mL)",
    },
    "iron": {
        "name": "Iron",
        "genes": {
            "HFE": {
                "risk_statuses": ["carrier", "homozygous"],
                "impact": "Increased iron absorption — risk of iron overload",
                "severity_map": {"carrier": -0.3, "homozygous": -0.8},  # negative = excess risk
            },
        },
        "food_sources": "Red meat, organ meats, shellfish, legumes, dark leafy greens",
        "supplement_form": "Iron bisglycinate (if deficient) — AVOID if HFE homozygous",
        "dose_range": "18-30 mg elemental iron/day (only if deficient)",
        "testing": "Ferritin, serum iron, TIBC, transferrin saturation",
    },
    "choline": {
        "name": "Choline",
        "genes": {
            "PEMT": {
                "risk_statuses": ["reduced"],
                "impact": "Reduced endogenous choline synthesis — higher dietary need",
                "severity_map": {"reduced": 0.6},
            },
            "MTHFR": {
                "risk_statuses": ["reduced", "severely_reduced"],
                "impact": "Impaired methylation increases choline demand as alternate methyl donor",
                "severity_map": {"reduced": 0.3, "severely_reduced": 0.5},
            },
        },
        "food_sources": "Eggs (yolks), liver, beef, fish, cruciferous vegetables",
        "supplement_form": "CDP-choline (citicoline) or phosphatidylcholine",
        "dose_range": "250-500 mg citicoline/day",
        "testing": "No standard clinical test; assess dietary intake",
    },
    "omega3": {
        "name": "Omega-3 (EPA/DHA)",
        "genes": {
            "FADS1": {
                "risk_statuses": ["reduced", "low_conversion"],
                "impact": "Reduced ALA→EPA→DHA conversion (FADS1/2 desaturase)",
                "severity_map": {"reduced": 0.5, "low_conversion": 0.7},
            },
        },
        "food_sources": "Fatty fish (salmon, sardines, mackerel), algae",
        "supplement_form": "Fish oil (EPA+DHA) or algal DHA",
        "dose_range": "1000-2000 mg combined EPA+DHA/day",
        "testing": "Omega-3 index blood test (target >8%)",
    },
    "magnesium": {
        "name": "Magnesium",
        "genes": {
            "COMT": {
                "risk_statuses": ["slow"],
                "impact": "Slow COMT needs magnesium as cofactor; increased demand",
                "severity_map": {"slow": 0.4},
            },
        },
        "food_sources": "Dark chocolate, nuts, seeds, leafy greens, avocado",
        "supplement_form": "Magnesium glycinate or threonate (best bioavailability)",
        "dose_range": "300-400 mg elemental magnesium/day",
        "testing": "RBC magnesium (not serum — serum is unreliable)",
    },
    "vitamin_a": {
        "name": "Vitamin A (Retinol)",
        "genes": {
            "BCMO1": {
                "risk_statuses": ["reduced", "poor_converter"],
                "impact": "Reduced beta-carotene to retinol conversion (up to 70% lower)",
                "severity_map": {"reduced": 0.5, "poor_converter": 0.8},
            },
        },
        "food_sources": "Liver, egg yolks, dairy (preformed); carrots, sweet potato (provitamin A)",
        "supplement_form": "Retinyl palmitate or preformed retinol (if poor converter)",
        "dose_range": "2500-5000 IU retinol/day (preformed, not beta-carotene)",
        "testing": "Serum retinol, retinol-binding protein",
    },
}


def profile_nutrigenomics(genome_by_rsid, lifestyle_findings=None):
    """Build a nutrigenomics profile from genome data and lifestyle findings.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict.
    lifestyle_findings : list or None
        Findings from analyze_lifestyle_health().

    Returns
    -------
    dict with keys: nutrient_needs, top_deficiency_risks, testing_recommendations,
                    supplement_priorities, summary.
    """
    if lifestyle_findings is None:
        lifestyle_findings = []

    # Build gene -> status lookup
    gene_status = {}
    for f in lifestyle_findings:
        gene = f.get("gene", "")
        status = f.get("status", "")
        if gene and status:
            gene_status[gene] = status

    nutrient_needs = []
    top_risks = []
    testing_recs = []
    supplement_priorities = []

    for nutrient_id, profile in NUTRIENT_PROFILES.items():
        total_severity = 0.0
        gene_impacts = []
        has_signal = False

        for gene, gene_info in profile["genes"].items():
            status = gene_status.get(gene, "")
            if status in gene_info["risk_statuses"]:
                severity = gene_info["severity_map"].get(status, 0.3)
                total_severity += severity
                has_signal = True
                gene_impacts.append({
                    "gene": gene,
                    "status": status,
                    "impact": gene_info["impact"],
                    "severity": severity,
                })

        if not has_signal:
            need_level = "normal"
            recommendation = "Standard dietary intake likely sufficient"
        elif total_severity >= 0.8:
            need_level = "high"
            recommendation = f"Supplementation recommended: {profile['supplement_form']}"
            top_risks.append(profile["name"])
            supplement_priorities.append({
                "nutrient": profile["name"],
                "form": profile["supplement_form"],
                "dose": profile["dose_range"],
                "priority": "high",
            })
        elif total_severity >= 0.4:
            need_level = "moderate"
            recommendation = f"Dietary focus + consider supplementation: {profile['supplement_form']}"
            supplement_priorities.append({
                "nutrient": profile["name"],
                "form": profile["supplement_form"],
                "dose": profile["dose_range"],
                "priority": "moderate",
            })
        elif total_severity < 0:
            need_level = "caution_excess"
            recommendation = "Caution: genetic tendency toward excess. Avoid supplementation unless deficient."
        else:
            need_level = "low"
            recommendation = "Minor genetic signal — dietary optimization sufficient"

        testing_recs.append({
            "nutrient": profile["name"],
            "test": profile["testing"],
        })

        nutrient_needs.append({
            "id": nutrient_id,
            "name": profile["name"],
            "need_level": need_level,
            "severity": round(total_severity, 2),
            "gene_impacts": gene_impacts,
            "food_sources": profile["food_sources"],
            "recommendation": recommendation,
        })

    # Sort by severity (highest need first)
    nutrient_needs.sort(key=lambda n: -abs(n["severity"]))
    supplement_priorities.sort(key=lambda s: {"high": 0, "moderate": 1, "low": 2}.get(s["priority"], 3))

    high_count = sum(1 for n in nutrient_needs if n["need_level"] == "high")
    mod_count = sum(1 for n in nutrient_needs if n["need_level"] == "moderate")

    if high_count >= 2:
        summary = f"Your genetics suggest increased needs for {high_count} key nutrients."
    elif high_count == 1:
        summary = f"One nutrient shows high genetic need; {mod_count} moderate."
    else:
        summary = "No major nutrient deficiency risks detected from genetics."

    return {
        "nutrient_needs": nutrient_needs,
        "top_deficiency_risks": top_risks,
        "testing_recommendations": testing_recs,
        "supplement_priorities": supplement_priorities,
        "summary": summary,
    }
