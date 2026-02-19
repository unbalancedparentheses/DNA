"""Gene-gene interaction (epistasis) analysis.

Evaluates known multi-gene interactions where the combined effect of
variants differs from each individual effect — e.g., MTHFR + COMT
combined methylation impact.

Only well-characterized, published interactions are included.
"""


# Each model defines:
#   name: human-readable interaction name
#   genes: set of genes involved
#   conditions: list of condition dicts, each with:
#     required: dict mapping gene -> list of qualifying statuses
#     effect: description of the combined effect
#     risk_level: "low", "moderate", "high"
#     mechanism: brief mechanism explanation
#     actions: list of recommended actions

EPISTASIS_MODELS = {
    "mthfr_comt_methylation": {
        "name": "MTHFR + COMT: Methylation-Catecholamine Interaction",
        "genes": {"MTHFR", "COMT"},
        "conditions": [
            {
                "required": {
                    "MTHFR": ["reduced", "severely_reduced"],
                    "COMT": ["slow"],
                },
                "effect": (
                    "Reduced folate conversion (MTHFR) combined with slow "
                    "catecholamine clearance (COMT) creates a dual burden: "
                    "impaired methylation AND catecholamine accumulation. "
                    "This may amplify anxiety, stress sensitivity, and "
                    "supplement sensitivity (especially methyl donors)."
                ),
                "risk_level": "high",
                "mechanism": (
                    "MTHFR reduces 5-MTHF production, the methyl donor for "
                    "homocysteine-to-methionine conversion. COMT requires SAMe "
                    "(from methionine) to degrade catecholamines. When both "
                    "are impaired, methyl group supply is low AND demand is "
                    "backed up."
                ),
                "actions": [
                    "Start methylfolate at LOW dose (200-400mcg) — slow COMT "
                    "individuals may over-respond to methyl donors",
                    "Prioritize magnesium glycinate (300-400mg) for COMT support",
                    "Avoid combining high-dose methylfolate + methylB12 initially",
                    "Stress management is critical — meditation, breathwork, sleep",
                    "Monitor for irritability/anxiety when starting methyl supplements",
                ],
            },
        ],
    },
    "mthfr_mtrr_methylation": {
        "name": "MTHFR + MTRR: Double Methylation Cycle Impairment",
        "genes": {"MTHFR", "MTRR"},
        "conditions": [
            {
                "required": {
                    "MTHFR": ["reduced", "severely_reduced"],
                    "MTRR": ["reduced", "severely_reduced"],
                },
                "effect": (
                    "Both arms of the methylation cycle are impaired: "
                    "MTHFR reduces folate activation, MTRR impairs B12 "
                    "recycling. Homocysteine may accumulate from two "
                    "independent bottlenecks."
                ),
                "risk_level": "high",
                "mechanism": (
                    "MTHFR converts folate to 5-MTHF. MTRR regenerates "
                    "methylcobalamin, the cofactor for methionine synthase. "
                    "When both are impaired, the methionine synthase reaction "
                    "is doubly compromised."
                ),
                "actions": [
                    "Methylfolate (400-800mcg) + methylcobalamin (1000-5000mcg)",
                    "Add riboflavin (B2) 25-50mg — cofactor for MTHFR enzyme",
                    "Monitor homocysteine — target <10 \u03bcmol/L",
                    "Emphasize leafy greens, legumes, liver for dietary folate",
                    "Consider TMG (trimethylglycine) as alternate methyl donor",
                ],
            },
        ],
    },
    "caffeine_sensitivity": {
        "name": "CYP1A2 + COMT + ADORA2A: Caffeine Sensitivity Triad",
        "genes": {"CYP1A2", "COMT", "ADORA2A"},
        "conditions": [
            {
                "required": {
                    "CYP1A2": ["slow", "intermediate"],
                    "COMT": ["slow"],
                },
                "effect": (
                    "Slow caffeine metabolism (CYP1A2) combined with slow "
                    "catecholamine clearance (COMT) means caffeine both stays "
                    "in the body longer AND its stimulatory effects are amplified "
                    "by delayed norepinephrine/dopamine breakdown."
                ),
                "risk_level": "moderate",
                "mechanism": (
                    "CYP1A2 metabolizes caffeine in the liver. COMT degrades "
                    "the catecholamines that caffeine stimulates. Slow variants "
                    "in both create a double hit: more caffeine exposure AND "
                    "more neurotransmitter accumulation."
                ),
                "actions": [
                    "Limit caffeine to <100mg/day (one small coffee) before 10am",
                    "Consider switching to green tea (L-theanine counterbalances)",
                    "Avoid caffeine entirely if experiencing anxiety or insomnia",
                    "Allow 10+ hours caffeine clearance before sleep",
                ],
            },
            {
                "required": {
                    "CYP1A2": ["slow", "intermediate"],
                    "ADORA2A": ["anxiety_prone"],
                },
                "effect": (
                    "Slow caffeine metabolism with anxiety-prone adenosine "
                    "receptors: caffeine stays active longer AND triggers "
                    "stronger anxiogenic response."
                ),
                "risk_level": "moderate",
                "mechanism": (
                    "CYP1A2 prolongs caffeine half-life. ADORA2A variant "
                    "increases anxiety response to adenosine receptor blockade "
                    "by caffeine."
                ),
                "actions": [
                    "Strong caffeine sensitivity likely — consider elimination",
                    "If consuming caffeine, maximum 50-100mg before 9am",
                    "L-theanine (200mg) may help if caffeine is consumed",
                ],
            },
        ],
    },
    "blood_pressure_polygenic": {
        "name": "ACE + AGT + AGTR1: Blood Pressure Polygenic Risk",
        "genes": {"ACE", "AGT", "AGTR1", "GNB3"},
        "conditions": [
            {
                "required": {
                    "ACE": ["elevated", "high"],
                    "AGT": ["elevated", "high"],
                },
                "effect": (
                    "Elevated ACE activity combined with elevated angiotensinogen "
                    "production amplifies the renin-angiotensin system — "
                    "substantially increases hypertension risk beyond either "
                    "variant alone."
                ),
                "risk_level": "high",
                "mechanism": (
                    "AGT produces angiotensinogen (the substrate). ACE converts "
                    "it to angiotensin II (the vasoconstrictor). When both are "
                    "elevated, more substrate AND more conversion = higher "
                    "angiotensin II levels."
                ),
                "actions": [
                    "Home blood pressure monitoring (target <130/80)",
                    "DASH diet — high potassium, low sodium (<2300mg/day)",
                    "150+ minutes/week aerobic exercise",
                    "Discuss ACE inhibitor suitability with physician if needed",
                    "Limit alcohol, manage stress",
                ],
            },
            {
                "required": {
                    "ACE": ["elevated", "high"],
                    "AGTR1": ["elevated", "high"],
                },
                "effect": (
                    "Elevated ACE activity combined with more sensitive "
                    "angiotensin II receptor — more angiotensin II produced "
                    "AND stronger vascular response to it."
                ),
                "risk_level": "high",
                "mechanism": (
                    "ACE increases angiotensin II production. AGTR1 variant "
                    "enhances receptor response. Combined effect on "
                    "vasoconstriction is multiplicative."
                ),
                "actions": [
                    "Blood pressure monitoring essential",
                    "Sodium restriction particularly important",
                    "ARBs (angiotensin receptor blockers) may be especially effective",
                    "Regular aerobic exercise — 30+ min most days",
                ],
            },
        ],
    },
    "apoe_mthfr_cardiovascular": {
        "name": "APOE e4 + MTHFR: Cardiovascular Risk Interaction",
        "genes": {"APOE", "MTHFR"},
        "conditions": [
            {
                "required": {
                    "APOE": ["e3_e4", "e4_e4"],
                    "MTHFR": ["reduced", "severely_reduced"],
                },
                "effect": (
                    "APOE e4 increases LDL cholesterol and cardiovascular risk. "
                    "Combined with MTHFR-driven elevated homocysteine, both "
                    "lipid and homocysteine cardiovascular risk pathways are "
                    "activated simultaneously."
                ),
                "risk_level": "high",
                "mechanism": (
                    "APOE e4 impairs LDL receptor recycling (higher LDL). "
                    "MTHFR reduces homocysteine clearance. Both are independent "
                    "cardiovascular risk factors with potential synergy."
                ),
                "actions": [
                    "Regular lipid panel (focus on LDL-C and apoB)",
                    "Monitor homocysteine — target <10 \u03bcmol/L",
                    "Mediterranean diet (reduce saturated fat, increase omega-3)",
                    "Exercise regularly — particularly beneficial for APOE e4",
                    "Consider early statin discussion with physician if LDL elevated",
                ],
            },
        ],
    },
    "metabolic_risk": {
        "name": "FTO + TCF7L2: Metabolic Syndrome Risk",
        "genes": {"FTO", "TCF7L2"},
        "conditions": [
            {
                "required": {
                    "FTO": ["risk", "high_risk"],
                    "TCF7L2": ["risk", "high_risk"],
                },
                "effect": (
                    "FTO increases obesity susceptibility, TCF7L2 increases "
                    "diabetes risk. Combined, the metabolic syndrome risk "
                    "is amplified: weight gain tendency feeds directly into "
                    "insulin resistance pathway."
                ),
                "risk_level": "moderate",
                "mechanism": (
                    "FTO variants increase appetite/reduced satiety signaling. "
                    "TCF7L2 impairs insulin secretion from beta cells. "
                    "Excess weight (FTO) worsens insulin resistance, which "
                    "compounds the TCF7L2-driven beta cell dysfunction."
                ),
                "actions": [
                    "Weight management is particularly important",
                    "Monitor fasting glucose and HbA1c annually",
                    "Regular physical activity (improves insulin sensitivity)",
                    "Emphasize low-glycemic-index foods",
                    "Consider time-restricted eating if appropriate",
                ],
            },
        ],
    },
    "bdnf_comt_stress": {
        "name": "BDNF + COMT: Stress Response Interaction",
        "genes": {"BDNF", "COMT"},
        "conditions": [
            {
                "required": {
                    "BDNF": ["reduced", "met_carrier"],
                    "COMT": ["slow"],
                },
                "effect": (
                    "Reduced BDNF production (Met carrier) combined with slow "
                    "COMT creates vulnerability to stress-related cognitive "
                    "effects: less neuroplasticity reserve AND more "
                    "catecholamine accumulation under stress."
                ),
                "risk_level": "moderate",
                "mechanism": (
                    "BDNF Val66Met reduces activity-dependent BDNF secretion. "
                    "COMT Val158Met (slow) reduces dopamine/norepinephrine "
                    "degradation. Under stress, reduced BDNF impairs "
                    "hippocampal function while excess catecholamines impair "
                    "prefrontal function."
                ),
                "actions": [
                    "Regular aerobic exercise is critical (strongest BDNF booster)",
                    "Stress management: meditation, breathwork, adequate sleep",
                    "Avoid chronic high-stress environments when possible",
                    "Social connection supports BDNF signaling",
                    "Consider mindfulness-based stress reduction (MBSR)",
                ],
            },
        ],
    },
    "iron_inflammation": {
        "name": "HFE + IL6: Iron-Inflammation Interaction",
        "genes": {"HFE", "IL6"},
        "conditions": [
            {
                "required": {
                    "HFE": ["carrier", "homozygous"],
                    "IL6": ["high"],
                },
                "effect": (
                    "HFE variant (iron accumulation tendency) combined with "
                    "elevated IL-6 (inflammatory) creates a cycle: excess iron "
                    "promotes inflammation, and inflammation alters iron "
                    "distribution (hepcidin axis)."
                ),
                "risk_level": "moderate",
                "mechanism": (
                    "HFE mutations impair hepcidin response, allowing excess iron "
                    "absorption. IL-6 variant increases baseline inflammation. "
                    "Excess iron catalyzes reactive oxygen species (Fenton reaction), "
                    "amplifying inflammatory damage."
                ),
                "actions": [
                    "Monitor ferritin and CRP together (iron + inflammation)",
                    "Consider regular blood donation if ferritin runs high",
                    "Anti-inflammatory diet: omega-3, colorful vegetables",
                    "Avoid iron supplements unless deficiency confirmed",
                    "Vitamin C with meals increases iron absorption — moderate intake",
                ],
            },
        ],
    },
}


# Severity weights for dosage-aware epistasis scoring.
# Higher = more severe status for the gene.
_STATUS_SEVERITY = {
    # MTHFR
    "severely_reduced": 1.0, "reduced": 0.5,
    # COMT
    "slow": 1.0, "intermediate": 0.5,
    # CYP1A2
    "slow": 1.0,
    # ADORA2A
    "anxiety_prone": 0.8,
    # ACE/AGT/AGTR1
    "high": 1.0, "elevated": 0.7,
    # APOE
    "e4_e4": 1.0, "e3_e4": 0.5,
    # FTO/TCF7L2
    "high_risk": 1.0, "risk": 0.5,
    # BDNF
    "met_carrier": 0.6,
    # HFE
    "homozygous": 1.0, "carrier": 0.5,
    # IL6
    "high": 1.0,
}


def evaluate_epistasis(findings: list) -> list:
    """Evaluate gene-gene interactions from lifestyle/health findings.

    Uses dosage-aware matching: each matched gene contributes a severity
    score (0-1), and the interaction's effective strength is the product
    of gene severities. This prevents binary matching that treats mild
    and severe variants identically.

    Args:
        findings: List of finding dicts from analyze_lifestyle_health(),
                  each with 'gene', 'status', and 'magnitude' keys.

    Returns:
        List of interaction dicts with keys: id, name, genes_involved,
        effect, risk_level, mechanism, actions, severity_score.
    """
    # Build gene -> set of statuses + max magnitude lookup
    gene_status = {}
    gene_magnitude = {}
    for f in findings:
        gene = f.get('gene', '')
        status = f.get('status', '')
        mag = f.get('magnitude', 0)
        if gene and status:
            if gene not in gene_status:
                gene_status[gene] = set()
            gene_status[gene].add(status)
            gene_magnitude[gene] = max(gene_magnitude.get(gene, 0), mag)

    interactions = []

    for model_id, model in EPISTASIS_MODELS.items():
        for condition in model['conditions']:
            # Check if all required genes have qualifying statuses
            matched_genes = {}
            gene_severities = {}
            match = True
            for gene, qualifying_statuses in condition['required'].items():
                if gene not in gene_status:
                    match = False
                    break
                overlap = gene_status[gene] & set(qualifying_statuses)
                if not overlap:
                    match = False
                    break
                matched_genes[gene] = sorted(overlap)
                # Severity: max of status severity weights for matched statuses
                severity = max(
                    _STATUS_SEVERITY.get(s, 0.5) for s in overlap
                )
                # Also factor in magnitude (0-6 scale, normalized to 0-1)
                mag_factor = min(gene_magnitude.get(gene, 3) / 6.0, 1.0)
                gene_severities[gene] = severity * (0.5 + 0.5 * mag_factor)

            if match:
                # Overall interaction severity: geometric mean of gene severities
                severity_product = 1.0
                for s in gene_severities.values():
                    severity_product *= s
                severity_score = severity_product ** (1.0 / len(gene_severities))

                # Adjust risk level based on severity score
                base_risk = condition['risk_level']
                if severity_score < 0.3 and base_risk == 'high':
                    effective_risk = 'moderate'
                elif severity_score >= 0.7 and base_risk == 'moderate':
                    effective_risk = 'high'
                else:
                    effective_risk = base_risk

                interactions.append({
                    'id': model_id,
                    'name': model['name'],
                    'genes_involved': matched_genes,
                    'effect': condition['effect'],
                    'risk_level': effective_risk,
                    'mechanism': condition['mechanism'],
                    'actions': condition['actions'],
                    'severity_score': round(severity_score, 2),
                })

    # Sort by risk level (high first), then by severity score
    risk_order = {'high': 0, 'moderate': 1, 'low': 2}
    interactions.sort(key=lambda x: (risk_order.get(x['risk_level'], 3), -x['severity_score']))

    return interactions
