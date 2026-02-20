"""Preventive care timeline generator.

Creates age-based screening recommendations personalized to genetic risk
profile. Combines PRS risk categories, APOE status, ACMG findings,
carrier status, and pharmacogenomic results into a single timeline.
"""


# Base screening schedule (general population guidelines)
# Each entry: {test, start_age, frequency, condition, source}

BASE_SCREENINGS = [
    {"test": "Blood pressure", "start_age": 18, "frequency": "Every 1-2 years", "condition": "hypertension", "source": "USPSTF"},
    {"test": "Lipid panel", "start_age": 35, "frequency": "Every 5 years", "condition": "cardiovascular", "source": "USPSTF"},
    {"test": "Fasting glucose / HbA1c", "start_age": 35, "frequency": "Every 3 years", "condition": "type2_diabetes", "source": "USPSTF/ADA"},
    {"test": "Colorectal cancer screening", "start_age": 45, "frequency": "Every 10 years (colonoscopy)", "condition": "colorectal_cancer", "source": "USPSTF"},
    {"test": "Mammography", "start_age": 50, "frequency": "Every 2 years", "condition": "breast_cancer", "source": "USPSTF"},
    {"test": "PSA test (discuss with doctor)", "start_age": 55, "frequency": "Shared decision", "condition": "prostate_cancer", "source": "USPSTF"},
    {"test": "Eye exam (dilated)", "start_age": 60, "frequency": "Every 1-2 years", "condition": "age_related_macular_degeneration", "source": "AAO"},
    {"test": "DEXA bone density scan", "start_age": 65, "frequency": "Every 2 years", "condition": "osteoporosis", "source": "USPSTF"},
]

# Genetic risk modifiers: how genetic findings change the screening schedule
# Each: {trigger, test_modification, reason}
GENETIC_MODIFIERS = {
    # PRS-based earlier screening
    "prs_elevated": {
        "hypertension": {"start_age_delta": -10, "frequency": "Annually", "reason": "Elevated genetic risk for hypertension"},
        "type2_diabetes": {"start_age_delta": -10, "frequency": "Annually", "reason": "Elevated genetic risk for type 2 diabetes"},
        "coronary_artery_disease": {"start_age_delta": -10, "frequency": "Every 2 years", "reason": "Elevated genetic risk for coronary artery disease"},
        "breast_cancer": {"start_age_delta": -10, "frequency": "Annually", "reason": "Elevated genetic risk for breast cancer"},
        "colorectal_cancer": {"start_age_delta": -5, "frequency": "Every 5 years", "reason": "Elevated genetic risk for colorectal cancer"},
        "prostate_cancer": {"start_age_delta": -10, "frequency": "Every 2 years", "reason": "Elevated genetic risk for prostate cancer"},
        "age_related_macular_degeneration": {"start_age_delta": -10, "frequency": "Annually", "reason": "Elevated genetic risk for AMD"},
        "osteoporosis": {"start_age_delta": -10, "frequency": "Every 2 years", "reason": "Elevated genetic risk for osteoporosis"},
        "parkinsons_disease": {"start_age_delta": 0, "frequency": "Discuss with neurologist", "reason": "Elevated genetic risk for Parkinson's disease"},
        "celiac_disease": {"start_age_delta": 0, "frequency": "If symptomatic", "reason": "Elevated genetic risk for celiac disease — test tTG-IgA if GI symptoms"},
        "rheumatoid_arthritis": {"start_age_delta": 0, "frequency": "If symptomatic", "reason": "Elevated genetic risk for RA — test anti-CCP if joint symptoms"},
    },
    "prs_high": {
        "hypertension": {"start_age_delta": -15, "frequency": "Every 6 months", "reason": "High genetic risk for hypertension"},
        "type2_diabetes": {"start_age_delta": -15, "frequency": "Every 6 months", "reason": "High genetic risk for type 2 diabetes"},
        "coronary_artery_disease": {"start_age_delta": -15, "frequency": "Annually", "reason": "High genetic risk for CAD — consider calcium score at 40"},
        "breast_cancer": {"start_age_delta": -15, "frequency": "Annually + MRI", "reason": "High genetic risk for breast cancer"},
        "colorectal_cancer": {"start_age_delta": -10, "frequency": "Every 3 years", "reason": "High genetic risk for colorectal cancer"},
    },
}

# APOE-specific modifications
APOE_MODIFIERS = {
    "elevated": [
        {"test": "Lipid panel", "start_age": 25, "frequency": "Every 2 years", "reason": "APOE e4 carrier — earlier lipid monitoring recommended"},
        {"test": "Cognitive screening", "start_age": 55, "frequency": "Every 2 years", "reason": "APOE e4 — proactive cognitive health monitoring"},
    ],
    "high": [
        {"test": "Lipid panel", "start_age": 20, "frequency": "Annually", "reason": "APOE e4/e4 — aggressive lipid management recommended"},
        {"test": "Cognitive screening", "start_age": 50, "frequency": "Annually", "reason": "APOE e4/e4 — early cognitive screening advised"},
        {"test": "Amyloid PET (discuss)", "start_age": 55, "frequency": "Discuss with neurologist", "reason": "APOE e4/e4 — consider amyloid imaging"},
    ],
}


def generate_preventive_timeline(prs_results=None, apoe=None, acmg=None,
                                  star_alleles=None, carrier_screen=None):
    """Generate personalized preventive care timeline.

    Parameters
    ----------
    prs_results : dict, optional
        PRS results from calculate_prs().
    apoe : dict, optional
        APOE haplotype from call_apoe_haplotype().
    acmg : dict, optional
        ACMG findings from flag_acmg_findings().
    star_alleles : dict, optional
        Star allele results.
    carrier_screen : dict, optional
        Carrier screening results.

    Returns
    -------
    dict with keys:
        timeline: list of {test, start_age, frequency, reason, priority, genetic_basis}
        summary: text summary
        early_screenings: count of screenings moved earlier
    """
    timeline = []
    early_count = 0

    # Start with base screenings
    base_map = {}
    for s in BASE_SCREENINGS:
        base_map[s["condition"]] = {
            "test": s["test"],
            "start_age": s["start_age"],
            "frequency": s["frequency"],
            "reason": f"General population guideline ({s['source']})",
            "priority": "standard",
            "genetic_basis": None,
        }

    # Apply PRS modifiers
    if prs_results:
        for cid, result in prs_results.items():
            cat = result.get("risk_category", "average")
            modifier_key = None
            if cat == "high":
                modifier_key = "prs_high"
            elif cat == "elevated":
                modifier_key = "prs_elevated"

            if modifier_key and modifier_key in GENETIC_MODIFIERS:
                mods = GENETIC_MODIFIERS[modifier_key]
                if cid in mods:
                    mod = mods[cid]
                    if cid in base_map:
                        base_map[cid]["start_age"] = max(18, base_map[cid]["start_age"] + mod["start_age_delta"])
                        base_map[cid]["frequency"] = mod["frequency"]
                        base_map[cid]["reason"] = mod["reason"]
                        base_map[cid]["priority"] = "high" if cat == "high" else "elevated"
                        base_map[cid]["genetic_basis"] = f"PRS {result['percentile']:.0f}th percentile ({cat})"
                        early_count += 1
                    else:
                        # Add new screening for conditions not in base
                        timeline.append({
                            "test": f"{result['name']} screening",
                            "start_age": 30,
                            "frequency": mod["frequency"],
                            "reason": mod["reason"],
                            "priority": "elevated",
                            "genetic_basis": f"PRS {result['percentile']:.0f}th percentile ({cat})",
                        })
                        early_count += 1

    # Apply APOE modifiers
    if apoe and apoe.get("risk_level") in APOE_MODIFIERS:
        for mod in APOE_MODIFIERS[apoe["risk_level"]]:
            timeline.append({
                "test": mod["test"],
                "start_age": mod["start_age"],
                "frequency": mod["frequency"],
                "reason": mod["reason"],
                "priority": "high" if apoe["risk_level"] == "high" else "elevated",
                "genetic_basis": f"APOE {apoe['apoe_type']}",
            })
            early_count += 1

    # ACMG findings — add specific screenings
    if acmg:
        for finding in acmg.get("acmg_findings", []):
            gene = finding.get("gene", "")
            if gene in ("BRCA1", "BRCA2"):
                timeline.append({
                    "test": "Breast MRI + mammography",
                    "start_age": 25,
                    "frequency": "Every 6 months (alternating)",
                    "reason": f"Pathogenic {gene} variant — enhanced breast cancer screening per NCCN",
                    "priority": "urgent",
                    "genetic_basis": f"{gene} pathogenic variant",
                })
                timeline.append({
                    "test": "Ovarian cancer screening (CA-125 + ultrasound)",
                    "start_age": 30,
                    "frequency": "Every 6-12 months",
                    "reason": f"Pathogenic {gene} variant — discuss risk-reducing surgery",
                    "priority": "urgent",
                    "genetic_basis": f"{gene} pathogenic variant",
                })
                early_count += 2
            elif gene in ("MLH1", "MSH2", "MSH6", "PMS2"):
                timeline.append({
                    "test": "Colonoscopy",
                    "start_age": 20,
                    "frequency": "Every 1-2 years",
                    "reason": f"Lynch syndrome ({gene}) — early colonoscopy per NCCN",
                    "priority": "urgent",
                    "genetic_basis": f"{gene} pathogenic variant",
                })
                early_count += 1

    # Pharmacogenomic card reminder
    if star_alleles:
        non_normal = [g for g, r in star_alleles.items()
                      if r.get("phenotype") not in ("normal", "Unknown")]
        if non_normal:
            timeline.append({
                "test": "Pharmacogenomic card review (bring to every prescriber)",
                "start_age": 18,
                "frequency": "Every visit",
                "reason": f"Non-standard metabolizer for: {', '.join(non_normal)}",
                "priority": "ongoing",
                "genetic_basis": "Pharmacogenomic profile",
            })

    # Add base screenings to timeline
    for cid, entry in base_map.items():
        timeline.append(entry)

    # Sort by start_age, then priority
    priority_order = {"urgent": 0, "high": 1, "elevated": 2, "ongoing": 3, "standard": 4}
    timeline.sort(key=lambda x: (x["start_age"], priority_order.get(x["priority"], 5)))

    if early_count:
        summary = (f"{early_count} screening(s) recommended earlier or more frequently "
                   f"than general population guidelines based on your genetic profile.")
    else:
        summary = "Your genetic profile does not indicate need for earlier screenings beyond standard guidelines."

    return {
        "timeline": timeline,
        "summary": summary,
        "early_screenings": early_count,
    }
