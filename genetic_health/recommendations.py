"""Personalized recommendations synthesis.

Analyzes ALL results together to detect convergent risk -- when multiple
independent signals (gene variants, PRS, pathogenic variants, epistasis)
point to the same condition.  Produces prioritized, actionable output
that a genetic counselor would give.
"""

from collections import defaultdict


# =========================================================================
# RISK GROUP DEFINITIONS
# =========================================================================
# Each group maps a health area to the genes, PRS conditions, and disease
# keywords that feed into it.  When multiple signals converge on the same
# group, priority is elevated.

RISK_GROUPS = {
    "blood_pressure": {
        "title": "Blood Pressure & Hypertension",
        "genes": {"AGTR1", "AGT", "ACE", "ADD1", "GNB3", "ADRB1"},
        "prs_condition": "hypertension",
        "disease_keywords": ["salt-sensitive", "hypertension", "blood pressure"],
        "actions": [
            "Home blood pressure monitoring (target <130/80 mmHg)",
            "Sodium intake <2,000 mg/day (DASH diet pattern)",
            "150+ min/week aerobic exercise (brisk walking, cycling, swimming)",
            "Potassium-rich foods: bananas, sweet potatoes, spinach",
            "Limit alcohol to <= 1 drink/day",
        ],
        "doctor_note": (
            "Multiple genetic signals converge on the renin-angiotensin-aldosterone "
            "system.  ACE inhibitors or ARBs may be especially effective given the "
            "RAAS pathway variants.  Consider earlier pharmacotherapy if lifestyle "
            "measures insufficient."
        ),
        "monitoring": [
            {"test": "Blood pressure", "frequency": "Weekly (home)", "reason": "Multiple BP gene variants"},
            {"test": "Renal function (BMP)", "frequency": "Annually", "reason": "Hypertension screening"},
        ],
    },
    "metabolic_diabetes": {
        "title": "Metabolic Health & Type 2 Diabetes",
        "genes": {"TCF7L2", "IRS2", "IGF2BP2", "CAPN10", "SOD2", "FTO"},
        "prs_condition": "t2d",
        "disease_keywords": ["diabetes", "insulin", "glucose"],
        "actions": [
            "Monitor fasting glucose and HbA1c annually",
            "Maintain healthy weight (BMI <25 or waist <94 cm men / <80 cm women)",
            "Low-glycaemic-index diet; emphasise whole grains, legumes, vegetables",
            "150+ min/week moderate exercise (improves insulin sensitivity)",
            "Consider time-restricted eating (e.g., 16:8) if appropriate",
        ],
        "doctor_note": (
            "TCF7L2 risk allele impairs beta-cell insulin secretion.  Combined with "
            "additional metabolic variants and/or elevated PRS, early metformin "
            "discussion is warranted if pre-diabetic (HbA1c 5.7-6.4%)."
        ),
        "monitoring": [
            {"test": "Fasting glucose / HbA1c", "frequency": "Annually", "reason": "Diabetes risk variants"},
            {"test": "Lipid panel", "frequency": "Annually", "reason": "Metabolic syndrome screening"},
        ],
    },
    "lung_health": {
        "title": "Lung & Respiratory Health",
        "genes": set(),
        "prs_condition": None,
        "disease_keywords": ["alpha-1 antitrypsin", "SERPINA1", "emphysema", "COPD"],
        "actions": [
            "Absolutely avoid smoking and second-hand smoke",
            "Request AAT (alpha-1 antitrypsin) level blood test",
            "Annual spirometry to track lung function (FEV1/FVC)",
            "Minimise exposure to air pollution, dust, chemical fumes",
            "Flu and pneumonia vaccinations up to date",
        ],
        "doctor_note": (
            "SERPINA1 pathogenic variant detected.  Measure serum AAT level.  If low, "
            "consider referral to pulmonology for baseline assessment and monitoring "
            "plan.  Augmentation therapy may be relevant for severe deficiency."
        ),
        "monitoring": [
            {"test": "Serum AAT level", "frequency": "Once (baseline)", "reason": "SERPINA1 pathogenic variant"},
            {"test": "Spirometry (FEV1/FVC)", "frequency": "Annually", "reason": "Lung function tracking"},
        ],
    },
    "caffeine": {
        "title": "Caffeine Sensitivity",
        "genes": {"CYP1A2", "ADORA2A", "COMT"},
        "prs_condition": None,
        "disease_keywords": [],
        "actions": [
            "Limit caffeine to <100 mg/day (1 small coffee) before 10 am",
            "Consider switching to green tea (L-theanine counterbalances jitters)",
            "Allow 10+ hours caffeine clearance before sleep",
            "Avoid caffeine entirely if experiencing anxiety or insomnia",
        ],
        "doctor_note": (
            "Slow CYP1A2 metabolism prolongs caffeine half-life.  Combined with "
            "COMT/ADORA2A variants, caffeine may trigger anxiety and raise "
            "cardiovascular stress.  Counsel patient on caffeine reduction."
        ),
        "monitoring": [],
    },
    "methylation": {
        "title": "Methylation & B-Vitamin Metabolism",
        "genes": {"MTHFR", "COMT", "MTRR", "PEMT"},
        "prs_condition": None,
        "disease_keywords": [],
        "actions": [
            "Methylfolate (L-5-MTHF) 400-800 mcg/day instead of folic acid",
            "Methylcobalamin (B12) 1,000 mcg sublingual",
            "Emphasise folate-rich foods: leafy greens, legumes, liver",
            "Avoid synthetic folic acid in fortified foods when possible",
            "Riboflavin (B2) 25-50 mg/day as MTHFR cofactor",
        ],
        "doctor_note": (
            "Multiple methylation cycle variants detected.  Monitor homocysteine "
            "(target <10 umol/L).  If slow COMT, start methyl donors at low dose "
            "and titrate up to avoid overmethylation symptoms."
        ),
        "monitoring": [
            {"test": "Homocysteine", "frequency": "Annually", "reason": "MTHFR / methylation variants"},
            {"test": "B12 + MMA", "frequency": "Annually", "reason": "MTRR / B12 recycling"},
        ],
    },
    "cardiovascular": {
        "title": "Cardiovascular & Clotting Risk",
        "genes": {"APOE", "F5", "F2", "LRP8", "CETP"},
        "prs_condition": "cad",
        "disease_keywords": ["myocardial infarction", "thrombosis", "thromboembolism",
                             "coronary", "heart disease"],
        "actions": [
            "Regular lipid panel (LDL-C, apoB, Lp(a) baseline)",
            "Mediterranean diet: olive oil, fish, nuts, vegetables",
            "150+ min/week aerobic exercise",
            "Know signs of DVT/PE if clotting variants present",
            "Stay hydrated and move frequently on long flights",
        ],
        "doctor_note": (
            "Convergent cardiovascular signals.  APOE e4 carriers benefit especially "
            "from exercise and saturated fat reduction.  If F5 Leiden or F2 present, "
            "discuss thromboprophylaxis for high-risk situations (surgery, immobility)."
        ),
        "monitoring": [
            {"test": "Lipid panel (LDL, apoB)", "frequency": "Annually", "reason": "Cardiovascular risk variants"},
            {"test": "Lp(a)", "frequency": "Once (baseline)", "reason": "Independent CV risk marker"},
        ],
    },
    "iron": {
        "title": "Iron Metabolism (Hemochromatosis Risk)",
        "genes": {"HFE"},
        "prs_condition": None,
        "disease_keywords": ["hemochromatosis", "iron overload"],
        "actions": [
            "Do NOT supplement iron unless deficiency is confirmed by blood test",
            "Consider regular blood donation to keep ferritin in range",
            "Moderate vitamin C intake with meals (enhances iron absorption)",
            "Annual ferritin and iron panel monitoring",
        ],
        "doctor_note": (
            "HFE variant detected.  Monitor ferritin and transferrin saturation.  "
            "If ferritin >300 ng/mL (men) or >200 ng/mL (women), consider "
            "therapeutic phlebotomy.  Screen family members."
        ),
        "monitoring": [
            {"test": "Ferritin + iron panel", "frequency": "Annually", "reason": "HFE variant"},
        ],
    },
    "skin": {
        "title": "Skin & UV Protection",
        "genes": {"MC1R"},
        "prs_condition": None,
        "disease_keywords": [],
        "actions": [
            "Daily SPF 30+ sunscreen on exposed skin (even on cloudy days)",
            "Wear protective clothing and hats during peak UV hours (10am-4pm)",
            "Annual full-body skin exam by dermatologist",
            "Topical retinoids for skin health (consult dermatologist)",
        ],
        "doctor_note": (
            "MC1R variant associated with reduced melanin, increased melanoma risk, "
            "and accelerated photoaging.  Annual dermatology screening recommended."
        ),
        "monitoring": [
            {"test": "Full-body skin exam", "frequency": "Annually", "reason": "MC1R variant"},
        ],
    },
    "cancer_screening": {
        "title": "Cancer Screening Awareness",
        "genes": {"TP53", "XRCC3"},
        "prs_condition": "breast_cancer",
        "disease_keywords": ["cancer", "carcinoma", "neoplasm", "tumor"],
        "actions": [
            "Follow age-appropriate cancer screening guidelines rigorously",
            "Discuss enhanced screening schedule with oncologist if high-risk",
            "Minimise known carcinogen exposure (tobacco, excess alcohol, UV)",
            "Maintain healthy weight (obesity is a modifiable cancer risk factor)",
        ],
        "doctor_note": (
            "Genetic signals suggest elevated cancer vigilance is warranted.  "
            "Consider referral to genetic counselor for comprehensive cancer "
            "risk assessment.  Enhanced screening may include earlier or more "
            "frequent imaging."
        ),
        "monitoring": [
            {"test": "Cancer screening (age-appropriate)", "frequency": "Per guidelines", "reason": "Cancer risk signals"},
        ],
    },
    "fitness": {
        "title": "Fitness & Exercise Response",
        "genes": {"ACTN3", "PPARGC1A", "COL5A1"},
        "prs_condition": None,
        "disease_keywords": [],
        "actions": [
            "Tailor training to your muscle fiber type (see ACTN3 result)",
            "Include both strength and endurance training regardless of genetics",
            "Pay attention to connective tissue health if COL5A1 variant present",
            "Progressive overload with adequate recovery between sessions",
        ],
        "doctor_note": "",
        "monitoring": [],
    },
}


# =========================================================================
# PRIORITY SCORING
# =========================================================================

def _compute_priority(signals):
    """Determine priority level from a list of signal types.

    signals: list of strings like 'pathogenic', 'prs_high', 'prs_elevated',
             'gene_signal', 'disease_keyword', 'epistasis'

    Returns 'high', 'moderate', or 'low'.
    """
    has_pathogenic = "pathogenic" in signals
    has_prs_high = "prs_high" in signals
    has_prs_elevated = "prs_elevated" in signals
    gene_count = signals.count("gene_signal")
    has_epistasis = "epistasis" in signals
    has_high_mag = "high_magnitude" in signals

    if has_pathogenic or has_prs_high or gene_count >= 3:
        return "high"
    if has_prs_elevated or gene_count >= 2 or has_high_mag or has_epistasis:
        return "moderate"
    return "low"


def _build_why(signals_detail):
    """Build a human-readable 'why' explanation from signal details."""
    parts = []
    for detail in signals_detail:
        parts.append(detail)
    return "; ".join(parts) if parts else "Single genetic signal detected."


# =========================================================================
# DRUG CARD
# =========================================================================

def _build_drug_card(findings, disease_findings):
    """Build a consolidated drug-gene interaction card.

    Groups CYP/drug-metabolism findings + ClinVar drug response by gene.
    """
    card = {}

    drug_genes = {"CYP2C9", "CYP2C19", "CYP1A2", "CYP2D6", "CYP3A4",
                  "VKORC1", "DPYD", "TPMT", "UGT1A1", "SLCO1B1"}

    for f in findings:
        gene = f.get("gene", "")
        if gene in drug_genes:
            if gene not in card:
                card[gene] = {"gene": gene, "entries": []}
            card[gene]["entries"].append({
                "rsid": f.get("rsid", ""),
                "genotype": f.get("genotype", ""),
                "status": f.get("status", ""),
                "description": f.get("description", ""),
                "source": "Lifestyle/SNP DB",
            })

    if disease_findings:
        for f in disease_findings.get("drug_response", []):
            gene = f.get("gene", "")
            if not gene:
                continue
            if gene not in card:
                card[gene] = {"gene": gene, "entries": []}
            card[gene]["entries"].append({
                "rsid": f.get("rsid", ""),
                "genotype": f.get("user_genotype", ""),
                "status": "",
                "description": (f.get("traits", "") or "Drug response variant"),
                "source": "ClinVar",
            })

    return list(card.values())


# =========================================================================
# GOOD NEWS
# =========================================================================

def _build_good_news(findings, disease_findings):
    """Collect protective and favorable findings."""
    good = []

    for f in findings:
        status = f.get("status", "")
        if status in ("protective", "longevity", "optimal", "fast"):
            good.append({
                "gene": f.get("gene", ""),
                "description": f.get("description", ""),
            })

    if disease_findings:
        for f in disease_findings.get("protective", []):
            good.append({
                "gene": f.get("gene", "Unknown"),
                "description": f.get("traits", "") or "Protective variant",
            })

    return good


# =========================================================================
# MAIN ENTRY POINT
# =========================================================================

def generate_recommendations(findings, disease_findings=None,
                             ancestry_results=None, prs_results=None,
                             epistasis_results=None):
    """Analyze ALL results to generate prioritized, convergent recommendations.

    Args:
        findings: list of lifestyle/health finding dicts
        disease_findings: dict with 'pathogenic', 'risk_factor', etc.
        ancestry_results: ancestry estimation dict (unused for now, reserved)
        prs_results: dict of PRS condition_id -> result dict
        epistasis_results: list of epistasis interaction dicts

    Returns:
        dict with keys: priorities, drug_card, monitoring_schedule, good_news
    """
    if findings is None:
        findings = []

    # Build quick lookups
    gene_findings = {}
    for f in findings:
        gene = f.get("gene", "")
        if gene:
            if gene not in gene_findings or f.get("magnitude", 0) > gene_findings[gene].get("magnitude", 0):
                gene_findings[gene] = f

    # Scan each risk group for convergent signals
    priorities = []

    for group_id, group in RISK_GROUPS.items():
        signals = []
        signals_detail = []

        # 1. Check gene signals
        matched_genes = []
        for gene in group["genes"]:
            if gene in gene_findings:
                gf = gene_findings[gene]
                mag = gf.get("magnitude", 0)
                status = gf.get("status", "")
                # Skip 'normal' / reference statuses
                if status in ("normal", "reference", "typical", ""):
                    continue
                signals.append("gene_signal")
                if mag >= 3:
                    signals.append("high_magnitude")
                matched_genes.append(f"{gene} ({status})")

        if matched_genes:
            signals_detail.append(f"Gene variants: {', '.join(matched_genes)}")

        # 2. Check PRS
        if prs_results and group["prs_condition"]:
            for cid, r in prs_results.items():
                if cid == group["prs_condition"]:
                    pct = r.get("percentile", 50)
                    cat = r.get("risk_category", "average")
                    if cat == "high":
                        signals.append("prs_high")
                        signals_detail.append(
                            f"PRS {r['name']}: {pct:.0f}th percentile (high)")
                    elif cat == "elevated":
                        signals.append("prs_elevated")
                        signals_detail.append(
                            f"PRS {r['name']}: {pct:.0f}th percentile (elevated)")

        # 3. Check disease findings (pathogenic/risk_factor)
        if disease_findings:
            for category in ("pathogenic", "likely_pathogenic"):
                for f in disease_findings.get(category, []):
                    gene_name = (f.get("gene", "") or "").upper()
                    traits = (f.get("traits", "") or "").lower()
                    # Match by gene in group genes
                    if gene_name in {g.upper() for g in group["genes"]}:
                        signals.append("pathogenic")
                        signals_detail.append(
                            f"Pathogenic variant: {f.get('gene', '')} "
                            f"({f.get('traits', 'unknown')})")
                    # Match by disease keywords
                    elif any(kw in traits for kw in group["disease_keywords"]):
                        signals.append("pathogenic")
                        signals_detail.append(
                            f"Pathogenic variant: {f.get('gene', '')} ({traits[:60]})")

            for f in disease_findings.get("risk_factor", []):
                traits = (f.get("traits", "") or "").lower()
                if any(kw in traits for kw in group["disease_keywords"]):
                    signals.append("disease_keyword")
                    # Only add detail for the first few
                    if signals_detail.count("disease_keyword") == 0:
                        signals_detail.append(
                            f"ClinVar risk factors related to "
                            f"{group['title'].lower()}")

        # 4. Check epistasis
        if epistasis_results:
            for interaction in epistasis_results:
                interaction_genes = set(interaction.get("genes_involved", {}).keys())
                if interaction_genes & group["genes"]:
                    signals.append("epistasis")
                    signals_detail.append(
                        f"Epistasis: {interaction['name']} "
                        f"({interaction['risk_level']})")

        # Only include groups with at least one real signal
        if not signals:
            continue

        priority = _compute_priority(signals)
        why = _build_why(signals_detail)

        priorities.append({
            "id": group_id,
            "title": group["title"],
            "priority": priority,
            "why": why,
            "actions": group["actions"],
            "doctor_note": group["doctor_note"],
            "monitoring": group["monitoring"],
            "signal_count": len(signals),
        })

    # Sort: high > moderate > low, then by signal count descending
    priority_order = {"high": 0, "moderate": 1, "low": 2}
    priorities.sort(key=lambda p: (priority_order.get(p["priority"], 3),
                                   -p["signal_count"]))

    # Build drug card
    drug_card = _build_drug_card(findings, disease_findings)

    # Build deduplicated monitoring schedule
    monitoring_schedule = _deduplicate_monitoring(priorities)

    # Build good news
    good_news = _build_good_news(findings, disease_findings)

    return {
        "priorities": priorities,
        "drug_card": drug_card,
        "monitoring_schedule": monitoring_schedule,
        "good_news": good_news,
    }


def _deduplicate_monitoring(priorities):
    """Collect monitoring items from all priorities, deduplicate by test name."""
    seen = {}
    freq_order = {
        "Weekly (home)": 0,
        "Monthly": 1,
        "Quarterly": 2,
        "Annually": 3,
        "Once (baseline)": 4,
        "Per guidelines": 5,
    }

    for p in priorities:
        for m in p.get("monitoring", []):
            test = m["test"]
            if test not in seen:
                seen[test] = m
            else:
                # Keep the more frequent one
                existing_freq = freq_order.get(seen[test]["frequency"], 10)
                new_freq = freq_order.get(m["frequency"], 10)
                if new_freq < existing_freq:
                    seen[test] = m

    # Sort by urgency (most frequent first)
    schedule = sorted(seen.values(),
                      key=lambda m: freq_order.get(m["frequency"], 10))
    return schedule
