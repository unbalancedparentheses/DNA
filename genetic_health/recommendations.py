"""Personalized recommendations synthesis.

Analyzes ALL results together to detect convergent risk -- when multiple
independent signals (gene variants, PRS, pathogenic variants, epistasis)
point to the same condition.  Produces prioritized, actionable output
that a genetic counselor would give.
"""

from collections import defaultdict

from .clinical_context import get_clinical_context, get_related_pathways, PATHWAYS


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
        "disease_keywords": ["alpha-1 antitrypsin", "alpha-1-antitrypsin",
                            "emphysema", "copd"],
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
        "genes": {"XRCC3"},
        "prs_condition": "breast_cancer",
        "disease_keywords": ["cancer", "carcinoma", "malignant neoplasm"],
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
# SPECIALIST REFERRAL TRIGGERS
# =========================================================================

SPECIALIST_REFERRALS = {
    "pulmonologist": {
        "title": "Pulmonologist",
        "triggers": [
            {"source": "acmg", "genes": ["SERPINA1"]},
            {"source": "disease", "keywords": ["alpha-1 antitrypsin", "COPD", "emphysema"]},
        ],
        "reason_template": "{gene}: {condition} — baseline pulmonary function assessment",
    },
    "genetic_counselor": {
        "title": "Genetic Counselor",
        "triggers": [
            {"source": "acmg", "any": True},
            {"source": "pathogenic_count", "min": 1},
        ],
        "reason_template": "Medically actionable genetic finding(s) detected",
    },
    "cardiologist": {
        "title": "Cardiologist",
        "triggers": [
            {"source": "priority", "group": "blood_pressure", "min_priority": "high"},
            {"source": "priority", "group": "cardiovascular", "min_priority": "high"},
        ],
        "reason_template": "Converging genetic signals in cardiovascular/BP pathways",
    },
    "hematologist": {
        "title": "Hematologist",
        "triggers": [
            {"source": "gene_status", "gene": "HFE",
             "statuses": ["carrier", "at_risk", "high_risk"]},
        ],
        "reason_template": "HFE variant — monitor iron studies, refer if ferritin elevated",
    },
    "pharmacist_pgx": {
        "title": "Pharmacist (PGx Review)",
        "triggers": [
            {"source": "star_alleles_actionable", "any": True},
        ],
        "reason_template": "Actionable pharmacogenomic variants: {genes}",
    },
    "dermatologist": {
        "title": "Dermatologist",
        "triggers": [
            {"source": "priority", "group": "skin", "min_priority": "low"},
        ],
        "reason_template": "MC1R variant — annual skin screening recommended",
    },
}


def _build_specialist_referrals(priorities, acmg, disease_findings,
                                star_alleles, gene_findings):
    """Evaluate specialist referral triggers and return matched referrals."""
    referrals = []
    priority_map = {p["id"]: p for p in priorities}
    priority_rank = {"high": 0, "moderate": 1, "low": 2}

    pathogenic_count = 0
    if disease_findings:
        pathogenic_count = (len(disease_findings.get("pathogenic", []))
                           + len(disease_findings.get("likely_pathogenic", [])))

    acmg_findings = (acmg or {}).get("acmg_findings", [])
    acmg_genes = {f.get("gene", "") for f in acmg_findings}

    actionable_star = []
    if star_alleles:
        for gene, r in star_alleles.items():
            if r.get("phenotype") in ("poor", "intermediate", "rapid", "ultrarapid"):
                actionable_star.append(gene)

    for ref_id, ref in SPECIALIST_REFERRALS.items():
        triggered = False
        reason = ref["reason_template"]

        for trigger in ref["triggers"]:
            src = trigger["source"]

            if src == "acmg":
                if trigger.get("any") and acmg_findings:
                    triggered = True
                elif trigger.get("genes"):
                    matched = acmg_genes & set(trigger["genes"])
                    if matched:
                        triggered = True
                        gene = next(iter(matched))
                        condition = next(
                            (f.get("traits", "variant detected")
                             for f in acmg_findings if f.get("gene") == gene),
                            "variant detected",
                        )
                        reason = reason.format(gene=gene, condition=condition)

            elif src == "pathogenic_count":
                if pathogenic_count >= trigger.get("min", 1):
                    triggered = True

            elif src == "disease":
                if disease_findings:
                    for cat in ("pathogenic", "likely_pathogenic"):
                        for f in disease_findings.get(cat, []):
                            traits = (f.get("traits", "") or "").lower()
                            if any(kw in traits for kw in trigger["keywords"]):
                                triggered = True
                                break
                        if triggered:
                            break

            elif src == "priority":
                group_id = trigger["group"]
                min_pri = trigger["min_priority"]
                if group_id in priority_map:
                    p = priority_map[group_id]
                    if priority_rank.get(p["priority"], 3) <= priority_rank.get(min_pri, 3):
                        triggered = True

            elif src == "gene_status":
                gene = trigger["gene"]
                if gene in gene_findings:
                    status = gene_findings[gene].get("status", "")
                    if status in trigger["statuses"]:
                        triggered = True

            elif src == "star_alleles_actionable":
                if trigger.get("any") and actionable_star:
                    triggered = True
                    reason = reason.format(genes=", ".join(actionable_star))

            if triggered:
                break

        if triggered:
            urgency = "routine"
            if ref_id in ("genetic_counselor", "pulmonologist"):
                urgency = "soon"
            referrals.append({
                "specialist": ref["title"],
                "reason": reason,
                "urgency": urgency,
            })

    return referrals


# =========================================================================
# CLINICAL CONTEXT INTEGRATION
# =========================================================================

def _build_clinical_insights(gene_findings):
    """Extract clinical context insights for findings with context entries.

    Returns a list of per-gene insight dicts (mechanism, implications, actions)
    for all findings that have clinical context data.
    """
    insights = []
    for gene, gf in gene_findings.items():
        status = gf.get("status", "")
        mag = gf.get("magnitude", 0)
        if mag < 1:
            continue
        ctx = get_clinical_context(gene, status)
        if not ctx:
            continue
        insights.append({
            "gene": gene,
            "status": status,
            "magnitude": mag,
            "mechanism": ctx.get("mechanism", ""),
            "implications": ctx.get("implications", []),
            "actions": ctx.get("actions", []),
            "interactions": ctx.get("interactions", []),
            "pathways": get_related_pathways(gene),
        })
    insights.sort(key=lambda x: -x["magnitude"])
    return insights


def _overlay_clinical_actions(priorities, gene_findings):
    """Augment priority group actions with specific clinical_context actions.

    For each priority group, find matching gene findings that have clinical
    context and append their specific actions (deduplicating against existing).
    """
    for p in priorities:
        group_id = p["id"]
        group_def = RISK_GROUPS.get(group_id)
        if not group_def:
            continue

        extra_actions = []
        for gene in group_def["genes"]:
            if gene not in gene_findings:
                continue
            gf = gene_findings[gene]
            status = gf.get("status", "")
            ctx = get_clinical_context(gene, status)
            if not ctx or not ctx.get("actions"):
                continue
            for action in ctx["actions"]:
                if action not in p["actions"] and action not in extra_actions:
                    extra_actions.append(action)

        if extra_actions:
            p["clinical_actions"] = extra_actions


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
                             epistasis_results=None,
                             star_alleles=None, acmg=None):
    """Analyze ALL results to generate prioritized, convergent recommendations.

    Args:
        findings: list of lifestyle/health finding dicts
        disease_findings: dict with 'pathogenic', 'risk_factor', etc.
        ancestry_results: ancestry estimation dict
        prs_results: dict of PRS condition_id -> result dict
        epistasis_results: list of epistasis interaction dicts
        star_alleles: dict of gene -> star allele result
        acmg: dict with ACMG secondary findings

    Returns:
        dict with keys: priorities, drug_card, monitoring_schedule,
        good_news, specialist_referrals, clinical_insights
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
                # Skip normal/reference statuses and low-magnitude polymorphisms
                if status in ("normal", "reference", "typical", "") or mag < 2:
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
                    # Use only the primary conditions (first 3 pipe-
                    # separated fields) to avoid matching on loosely
                    # associated somatic conditions that ClinVar appends
                    # later in the traits field (COSMIC annotations etc.)
                    full_traits = (f.get("traits", "") or "").lower()
                    parts = full_traits.split("|")
                    primary_trait = "|".join(parts[:3]).strip()
                    # Match by gene in group genes
                    if gene_name in {g.upper() for g in group["genes"]}:
                        signals.append("pathogenic")
                        signals_detail.append(
                            f"Pathogenic variant: {f.get('gene', '')} "
                            f"({f.get('traits', 'unknown')})")
                    # Match by disease keywords on primary condition only
                    # (require gold_stars >= 2 to avoid low-evidence entries)
                    elif f.get('gold_stars', 0) >= 2 and any(
                        kw in primary_trait for kw in group["disease_keywords"]
                    ):
                        signals.append("pathogenic")
                        signals_detail.append(
                            f"Pathogenic variant: {f.get('gene', '')} "
                            f"({primary_trait[:60]})")

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

    # Overlay clinical_context actions onto matching priority groups
    _overlay_clinical_actions(priorities, gene_findings)

    # Build drug card
    drug_card = _build_drug_card(findings, disease_findings)

    # Build deduplicated monitoring schedule
    monitoring_schedule = _deduplicate_monitoring(priorities)

    # Build good news
    good_news = _build_good_news(findings, disease_findings)

    # Build specialist referrals
    specialist_referrals = _build_specialist_referrals(
        priorities, acmg, disease_findings, star_alleles, gene_findings)

    # Build clinical insights
    clinical_insights = _build_clinical_insights(gene_findings)

    return {
        "priorities": priorities,
        "drug_card": drug_card,
        "monitoring_schedule": monitoring_schedule,
        "good_news": good_news,
        "specialist_referrals": specialist_referrals,
        "clinical_insights": clinical_insights,
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
