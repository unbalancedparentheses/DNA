"""Drug-specific dosing recommendations based on pharmacogenomic profile.

Maps genotype combinations to specific drug dosing guidance using CPIC
and DPWG guidelines. Each recommendation includes the drug, the relevant
genes, the genotype-specific action, and the evidence source.
"""


# Drug dosing rules: keyed by drug name
# Each rule has: genes (list of relevant genes), rules (list of condition -> recommendation),
# and a source citation.

DRUG_DOSING_RULES = {
    "warfarin": {
        "genes": ["CYP2C9", "VKORC1"],
        "category": "anticoagulant",
        "rules": [
            {
                "condition": lambda sa, f: (
                    sa.get("CYP2C9", {}).get("phenotype") in ("intermediate", "poor")
                ),
                "action": "Reduce warfarin starting dose by 25-50%. CYP2C9 intermediate/poor metabolizers clear warfarin more slowly, increasing bleeding risk.",
                "dose_guidance": "Consider starting at 2-3 mg/day instead of standard 5 mg/day.",
            },
            {
                "condition": lambda sa, f: any(
                    fi.get("gene") == "VKORC1" and fi.get("status") in ("high_sensitivity", "sensitive")
                    for fi in f
                ),
                "action": "VKORC1 sensitivity detected — warfarin target dose likely 1-3 mg/day. Requires careful INR monitoring.",
                "dose_guidance": "Start low (1-2 mg/day), titrate to INR 2.0-3.0.",
            },
        ],
        "source": "CPIC Guideline for Warfarin (Johnson et al. 2017, PMID: 28198005)",
    },
    "clopidogrel": {
        "genes": ["CYP2C19"],
        "category": "antiplatelet",
        "rules": [
            {
                "condition": lambda sa, f: (
                    sa.get("CYP2C19", {}).get("phenotype") in ("poor", "intermediate")
                ),
                "action": "Reduced CYP2C19 function — clopidogrel may not be adequately activated. Consider prasugrel or ticagrelor as alternatives.",
                "dose_guidance": "Avoid clopidogrel if possible. Use prasugrel 10mg or ticagrelor 90mg twice daily.",
            },
            {
                "condition": lambda sa, f: (
                    sa.get("CYP2C19", {}).get("phenotype") == "ultrarapid"
                ),
                "action": "Ultrarapid CYP2C19 — standard clopidogrel dosing is appropriate. Enhanced activation may increase bleeding risk at higher doses.",
                "dose_guidance": "Standard 75 mg/day is appropriate.",
            },
        ],
        "source": "CPIC Guideline for Clopidogrel (Scott et al. 2013, PMID: 23698643)",
    },
    "codeine": {
        "genes": ["CYP2D6"],
        "category": "analgesic",
        "rules": [
            {
                "condition": lambda sa, f: (
                    sa.get("CYP2D6", {}).get("phenotype") == "poor"
                ),
                "action": "CYP2D6 poor metabolizer — codeine will NOT be converted to morphine. No analgesic effect expected. Use alternative pain medication.",
                "dose_guidance": "Avoid codeine entirely. Use non-opioid analgesics or morphine/oxycodone at standard doses.",
            },
            {
                "condition": lambda sa, f: (
                    sa.get("CYP2D6", {}).get("phenotype") == "ultrarapid"
                ),
                "action": "CYP2D6 ultrarapid metabolizer — codeine is converted to morphine too rapidly, risking respiratory depression. AVOID codeine.",
                "dose_guidance": "CONTRAINDICATED. Use non-opioid alternatives or morphine at reduced doses with monitoring.",
            },
            {
                "condition": lambda sa, f: (
                    sa.get("CYP2D6", {}).get("phenotype") == "intermediate"
                ),
                "action": "CYP2D6 intermediate metabolizer — reduced codeine activation. May have suboptimal pain relief.",
                "dose_guidance": "Consider alternative analgesics if inadequate pain relief at standard doses.",
            },
        ],
        "source": "CPIC Guideline for Codeine (Crews et al. 2014, PMID: 24458010)",
    },
    "simvastatin": {
        "genes": ["SLCO1B1"],
        "category": "statin",
        "rules": [
            {
                "condition": lambda sa, f: (
                    sa.get("SLCO1B1", {}).get("phenotype") in ("intermediate", "poor")
                ),
                "action": "SLCO1B1 decreased function — elevated simvastatin myopathy risk (5-17x). Avoid simvastatin >20mg or switch statin.",
                "dose_guidance": "Use simvastatin ≤20 mg/day, or switch to rosuvastatin/pravastatin which are less affected.",
            },
        ],
        "source": "CPIC Guideline for Simvastatin (Ramsey et al. 2014, PMID: 24918167)",
    },
    "fluoropyrimidines": {
        "genes": ["DPYD"],
        "category": "chemotherapy",
        "rules": [
            {
                "condition": lambda sa, f: (
                    sa.get("DPYD", {}).get("phenotype") == "poor"
                ),
                "action": "DPYD deficient — fluoropyrimidine chemotherapy (5-FU, capecitabine) can be FATAL. Complete dose avoidance or alternative regimen required.",
                "dose_guidance": "CONTRAINDICATED at full dose. If no alternative, start at ≤25% dose with intensive monitoring.",
            },
            {
                "condition": lambda sa, f: (
                    sa.get("DPYD", {}).get("phenotype") == "intermediate"
                ),
                "action": "DPYD intermediate metabolizer — start fluoropyrimidines at 50% reduced dose with dose escalation based on tolerance.",
                "dose_guidance": "Start at 50% dose. Escalate in subsequent cycles if tolerated.",
            },
        ],
        "source": "CPIC Guideline for Fluoropyrimidines (Amstutz et al. 2018, PMID: 29152729)",
    },
    "thiopurines": {
        "genes": ["TPMT"],
        "category": "immunosuppressant",
        "rules": [
            {
                "condition": lambda sa, f: (
                    sa.get("TPMT", {}).get("phenotype") == "poor"
                ),
                "action": "TPMT deficient — standard thiopurine doses (azathioprine, 6-MP) cause life-threatening myelosuppression. Drastic dose reduction required.",
                "dose_guidance": "Reduce dose to 10% of standard. Start azathioprine at 0.5 mg/kg/day (vs standard 2-3 mg/kg/day).",
            },
            {
                "condition": lambda sa, f: (
                    sa.get("TPMT", {}).get("phenotype") == "intermediate"
                ),
                "action": "TPMT intermediate metabolizer — reduce thiopurine starting dose by 30-50%.",
                "dose_guidance": "Start at 50% dose (azathioprine ~1-1.5 mg/kg/day). Monitor CBC weekly for first month.",
            },
        ],
        "source": "CPIC Guideline for Thiopurines (Relling et al. 2019, PMID: 30447069)",
    },
    "tacrolimus": {
        "genes": ["CYP3A5"],
        "category": "immunosuppressant",
        "rules": [
            {
                "condition": lambda sa, f: (
                    sa.get("CYP3A5", {}).get("phenotype") == "normal"
                ),
                "action": "CYP3A5 expresser — requires higher tacrolimus doses to reach therapeutic levels. Standard dosing may be subtherapeutic.",
                "dose_guidance": "Start at 0.3 mg/kg/day (1.5-2x standard). Monitor trough levels closely.",
            },
            {
                "condition": lambda sa, f: (
                    sa.get("CYP3A5", {}).get("phenotype") in ("intermediate", "poor")
                ),
                "action": "CYP3A5 non-expresser — standard tacrolimus dosing appropriate.",
                "dose_guidance": "Start at 0.15 mg/kg/day (standard). Monitor trough levels per protocol.",
            },
        ],
        "source": "CPIC Guideline for Tacrolimus (Birdwell et al. 2015, PMID: 25801146)",
    },
    "isoniazid": {
        "genes": ["NAT2"],
        "category": "antibiotic",
        "rules": [
            {
                "condition": lambda sa, f: (
                    sa.get("NAT2", {}).get("phenotype") in ("poor", "intermediate")
                ),
                "action": "NAT2 slow acetylator — increased risk of isoniazid-induced hepatotoxicity and peripheral neuropathy. Supplement with pyridoxine.",
                "dose_guidance": "Standard dose but ADD pyridoxine (vitamin B6) 25-50 mg/day to prevent neuropathy. Monitor liver function.",
            },
        ],
        "source": "DPWG Guideline for Isoniazid (Swen et al. 2011, PMID: 21412232)",
    },
    "caffeine": {
        "genes": ["CYP1A2"],
        "category": "dietary",
        "rules": [
            {
                "condition": lambda sa, f: (
                    sa.get("CYP1A2", {}).get("phenotype") in ("poor", "intermediate")
                ),
                "action": "CYP1A2 slow metabolizer — caffeine stays active longer. Limit to 1-2 cups before noon to protect sleep and cardiovascular health.",
                "dose_guidance": "≤200 mg caffeine/day (~2 cups coffee), none after 12 PM.",
            },
            {
                "condition": lambda sa, f: (
                    sa.get("CYP1A2", {}).get("phenotype") in ("rapid", "ultrarapid", "normal")
                ),
                "action": "CYP1A2 normal/rapid metabolizer — standard caffeine clearance. Moderate intake (3-4 cups/day) generally well tolerated.",
                "dose_guidance": "Up to 400 mg caffeine/day (~4 cups) is generally safe.",
            },
        ],
        "source": "Cornelis et al. 2006, JAMA (PMID: 16522833)",
    },
}


def generate_drug_dosing(star_alleles, lifestyle_findings=None):
    """Generate drug-specific dosing recommendations from pharmacogenomic profile.

    Parameters
    ----------
    star_alleles : dict
        Star allele calling results from call_star_alleles().
    lifestyle_findings : list, optional
        Lifestyle findings from analyze_lifestyle_health().

    Returns
    -------
    dict with keys:
        recommendations: list of {drug, category, genes, action, dose_guidance, source}
        warnings: list of critical warnings (poor/no function findings)
        summary: text summary
    """
    if not star_alleles:
        return {
            "recommendations": [],
            "warnings": [],
            "summary": "No pharmacogenomic data available for drug dosing.",
        }

    findings = lifestyle_findings or []
    recommendations = []
    warnings = []

    for drug, rule_set in DRUG_DOSING_RULES.items():
        # Check if we have data for any relevant gene
        relevant_genes = rule_set["genes"]
        has_data = any(g in star_alleles for g in relevant_genes)
        if not has_data:
            continue

        for rule in rule_set["rules"]:
            try:
                if rule["condition"](star_alleles, findings):
                    rec = {
                        "drug": drug.replace("_", " ").title(),
                        "category": rule_set["category"],
                        "genes": relevant_genes,
                        "action": rule["action"],
                        "dose_guidance": rule["dose_guidance"],
                        "source": rule_set["source"],
                    }
                    recommendations.append(rec)

                    # Flag critical warnings
                    if any(w in rule["action"].upper() for w in
                           ["CONTRAINDICATED", "FATAL", "AVOID", "LIFE-THREATENING"]):
                        warnings.append(rec)
            except (KeyError, TypeError):
                continue

    n_recs = len(recommendations)
    n_warn = len(warnings)
    if n_warn:
        summary = f"{n_recs} drug dosing recommendations, including {n_warn} critical warning(s)."
    elif n_recs:
        summary = f"{n_recs} drug dosing recommendations based on your pharmacogenomic profile."
    else:
        summary = "No drug dosing adjustments needed based on available pharmacogenomic data."

    return {
        "recommendations": recommendations,
        "warnings": warnings,
        "summary": summary,
    }
