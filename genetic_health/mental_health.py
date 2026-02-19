"""Mental health susceptibility profiling.

Synthesizes neurotransmitter, stress response, and mood-related gene
variants into a mental health genetic profile with treatment matching.
"""


# Mental health-associated SNPs organized by domain
MENTAL_HEALTH_SNPS = {
    # --- Depression / Mood ---
    "rs6295": {
        "gene": "HTR1A",
        "domain": "depression",
        "risk_allele": "G",
        "description": "5-HT1A receptor (C-1019G) — reduced serotonin autoreceptor feedback",
        "effect": "GG genotype associated with increased depression risk and SSRI resistance",
        "reference": "Lemonde et al. 2003 (J Neurosci)",
    },
    "rs1800532": {
        "gene": "TPH1",
        "domain": "depression",
        "risk_allele": "A",
        "description": "Tryptophan hydroxylase 1 — serotonin synthesis rate-limiting enzyme",
        "effect": "A allele associated with reduced serotonin synthesis and depression risk",
        "reference": "Zill et al. 2004 (Mol Psychiatry)",
    },
    # --- Anxiety / Stress ---
    "rs110402": {
        "gene": "CRHR1",
        "domain": "anxiety",
        "risk_allele": "A",
        "description": "Corticotropin-releasing hormone receptor 1 — HPA axis regulation",
        "effect": "Modulates cortisol stress response; A allele may be protective post-trauma",
        "reference": "Bradley et al. 2008 (Arch Gen Psychiatry)",
    },
    "rs9296158": {
        "gene": "FKBP5",
        "domain": "anxiety",
        "risk_allele": "A",
        "description": "FK506-binding protein 5 — glucocorticoid receptor sensitivity",
        "effect": "A allele associated with impaired cortisol recovery after stress",
        "reference": "Binder et al. 2008 (JAMA)",
    },
    # --- Addiction ---
    "rs279858": {
        "gene": "GABRA2",
        "domain": "addiction",
        "risk_allele": "G",
        "description": "GABA-A receptor alpha-2 subunit — inhibitory neurotransmission",
        "effect": "G allele associated with alcohol dependence and externalizing behavior",
        "reference": "Edenberg et al. 2004 (Am J Hum Genet)",
    },
    "rs16969968": {
        "gene": "CHRNA5",
        "domain": "addiction",
        "risk_allele": "A",
        "description": "Nicotinic receptor alpha-5 — nicotine sensitivity",
        "effect": "A allele increases nicotine dependence risk and cigarettes/day",
        "reference": "Saccone et al. 2007 (Hum Mol Genet)",
    },
    # --- ADHD ---
    "rs1611115": {
        "gene": "DBH",
        "domain": "adhd",
        "risk_allele": "T",
        "description": "Dopamine beta-hydroxylase — converts dopamine to norepinephrine",
        "effect": "T allele reduces DBH activity; associated with ADHD traits",
        "reference": "Cubells et al. 2011 (Am J Med Genet B)",
    },
}

# Genes already in the lifestyle SNP database that contribute to mental health
_EXISTING_MENTAL_GENES = {
    "depression": {
        "BDNF": {"risk_statuses": ["reduced", "met_carrier"],
                  "impact": "Reduced neuroplasticity and BDNF secretion"},
        "SLC6A4": {"risk_statuses": ["short", "reduced"],
                    "impact": "Reduced serotonin reuptake transporter efficiency"},
        "MTHFR": {"risk_statuses": ["reduced", "severely_reduced"],
                   "impact": "Impaired methylation affects neurotransmitter synthesis"},
        "TNF": {"risk_statuses": ["high", "elevated"],
                "impact": "Neuroinflammation linked to treatment-resistant depression"},
        "IL6": {"risk_statuses": ["high", "elevated"],
                "impact": "Inflammatory cytokine elevation drives sickness behavior/anhedonia"},
    },
    "anxiety": {
        "COMT": {"risk_statuses": ["slow"],
                 "impact": "Elevated prefrontal catecholamines increase anxiety proneness"},
        "ADORA2A": {"risk_statuses": ["anxiety_prone"],
                     "impact": "Caffeine-induced anxiety sensitivity"},
    },
    "addiction": {
        "ANKK1": {"risk_statuses": ["risk", "reduced"],
                   "impact": "Reduced DRD2 density — reward deficiency, addiction vulnerability"},
        "OPRM1": {"risk_statuses": ["reduced", "altered"],
                   "impact": "Altered opioid receptor sensitivity"},
        "ALDH2": {"risk_statuses": ["reduced", "non_functional"],
                   "impact": "Alcohol flush reaction — protective against alcoholism"},
        "ADH1B": {"risk_statuses": ["fast"],
                   "impact": "Fast alcohol metabolism — modestly protective"},
    },
    "cognitive": {
        "COMT": {"risk_statuses": ["slow"],
                 "impact": "Higher prefrontal dopamine — better working memory but stress sensitivity"},
        "BDNF": {"risk_statuses": ["reduced", "met_carrier"],
                 "impact": "Reduced learning-dependent BDNF release"},
    },
}

# Treatment matching based on genetic profile
_TREATMENT_MATCHING = {
    "ssri_response": {
        "favorable": ["SLC6A4:normal", "CYP2C19:normal", "CYP2D6:normal", "HTR1A:CC"],
        "caution": ["SLC6A4:short", "CYP2C19:poor", "COMT:slow", "HTR1A:GG"],
        "note": "Poor CYP2C19 + slow COMT may cause SSRI hypersensitivity",
    },
    "exercise_therapy": {
        "favorable": ["BDNF:any", "COMT:slow", "PPARGC1A:any"],
        "note": "Exercise is the strongest BDNF booster; particularly effective for BDNF Met carriers",
    },
    "mindfulness": {
        "favorable": ["COMT:slow", "FKBP5:risk"],
        "note": "Slow COMT responders benefit most from stress reduction techniques",
    },
}


def profile_mental_health(genome_by_rsid, lifestyle_findings=None, star_alleles=None):
    """Build a mental health genetic profile.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict.
    lifestyle_findings : list or None
        Findings from analyze_lifestyle_health().
    star_alleles : dict or None
        Star allele results.

    Returns
    -------
    dict with keys: domains, new_markers, treatment_notes, resilience_factors,
                    risk_factors, recommendations, summary.
    """
    if lifestyle_findings is None:
        lifestyle_findings = []
    if star_alleles is None:
        star_alleles = {}

    # Build gene -> status from lifestyle
    gene_status = {}
    for f in lifestyle_findings:
        gene = f.get("gene", "")
        status = f.get("status", "")
        if gene and status:
            gene_status[gene] = status

    # --- Score new mental health SNPs ---
    new_markers = []
    for rsid, info in MENTAL_HEALTH_SNPS.items():
        entry = genome_by_rsid.get(rsid)
        if not entry:
            continue
        genotype = entry.get("genotype", "")
        if not genotype:
            continue
        copies = sum(1 for a in genotype if a == info["risk_allele"])
        new_markers.append({
            "rsid": rsid,
            "gene": info["gene"],
            "domain": info["domain"],
            "genotype": genotype,
            "risk_copies": copies,
            "description": info["description"],
            "effect": info["effect"],
            "reference": info["reference"],
        })

    # --- Score each domain ---
    domains = {}
    for domain_name, domain_genes in _EXISTING_MENTAL_GENES.items():
        risk_score = 0.0
        protective_score = 0.0
        signals = []

        for gene, gene_info in domain_genes.items():
            status = gene_status.get(gene, "")
            if status in gene_info["risk_statuses"]:
                risk_score += 1.0
                signals.append(f"{gene}: {gene_info['impact']}")
            elif status and status not in ("normal", "reference", "typical"):
                # Neutral signal
                pass
            else:
                protective_score += 0.3

        # Add new marker signals
        for m in new_markers:
            if m["domain"] == domain_name and m["risk_copies"] > 0:
                risk_score += m["risk_copies"] * 0.5
                signals.append(f"{m['gene']}: {m['effect']}")

        # Normalize
        total = risk_score + protective_score
        if total > 0:
            risk_pct = round(risk_score / (risk_score + max(protective_score, 1)) * 100)
        else:
            risk_pct = 30  # baseline

        if risk_pct >= 60:
            level = "elevated"
        elif risk_pct >= 40:
            level = "moderate"
        else:
            level = "low"

        domains[domain_name] = {
            "risk_level": level,
            "risk_score": risk_pct,
            "signals": signals[:5],
        }

    # --- Risk and resilience factors ---
    risk_factors = []
    resilience_factors = []

    if gene_status.get("COMT") == "slow":
        risk_factors.append("Slow COMT: higher anxiety proneness, but better working memory")
        resilience_factors.append("Slow COMT: enhanced cognitive performance under low stress")
    if gene_status.get("BDNF") in ("reduced", "met_carrier"):
        risk_factors.append("BDNF Met carrier: reduced neuroplasticity under stress")
        resilience_factors.append("BDNF Met: responds strongly to exercise-induced BDNF boost")
    if gene_status.get("ALDH2") in ("reduced", "non_functional"):
        resilience_factors.append("ALDH2 variant: natural alcohol aversion (protective)")
    if gene_status.get("SLC6A4") in ("short", "reduced"):
        risk_factors.append("Short serotonin transporter: increased environmental sensitivity")
        resilience_factors.append("Short 5-HTTLPR: also more responsive to positive environments")

    # Add new marker findings
    for m in new_markers:
        if m["risk_copies"] >= 2:
            risk_factors.append(f"{m['gene']}: Homozygous risk — {m['effect']}")
        elif m["risk_copies"] == 1 and m["domain"] != "addiction":
            risk_factors.append(f"{m['gene']}: Heterozygous — moderate {m['domain']} signal")

    # --- Treatment notes ---
    treatment_notes = []

    # SSRI matching
    cyp2c19_pheno = star_alleles.get("CYP2C19", {}).get("phenotype", "normal")
    if cyp2c19_pheno in ("poor", "intermediate"):
        treatment_notes.append(
            "CYP2C19 reduced: Start SSRIs (citalopram, escitalopram, sertraline) "
            "at 50% dose. Consider alternatives less dependent on CYP2C19."
        )
    elif cyp2c19_pheno in ("ultrarapid", "rapid"):
        treatment_notes.append(
            "CYP2C19 ultrarapid: Standard SSRI doses may be insufficient. "
            "Consider higher doses or non-CYP2C19 alternatives."
        )

    if gene_status.get("COMT") == "slow":
        treatment_notes.append(
            "Slow COMT: May be sensitive to stimulant medications and methyl donors. "
            "Mindfulness-based therapy particularly effective."
        )

    if gene_status.get("BDNF") in ("reduced", "met_carrier"):
        treatment_notes.append(
            "BDNF Met carrier: Exercise prescription is evidence-based — "
            "aerobic exercise 150+ min/week is as effective as SSRIs for mild-moderate depression."
        )

    # --- Recommendations ---
    recommendations = [
        "Regular aerobic exercise (strongest natural antidepressant — boosts BDNF, serotonin, dopamine)",
        "Sleep optimization (circadian disruption worsens all mental health domains)",
        "Social connection (protective across all genetic risk profiles)",
    ]

    if domains.get("anxiety", {}).get("risk_level") == "elevated":
        recommendations.append(
            "Stress management: meditation, breathwork, or CBT-based anxiety protocol"
        )
    if domains.get("addiction", {}).get("risk_level") == "elevated":
        recommendations.append(
            "Awareness of addiction vulnerability — set boundaries with alcohol/substances"
        )
    if any(m["gene"] == "FKBP5" and m["risk_copies"] > 0 for m in new_markers):
        recommendations.append(
            "Trauma-informed care: FKBP5 variant suggests altered stress recovery — "
            "consider EMDR or trauma-focused CBT if relevant"
        )

    # Summary
    elevated = [d for d, info in domains.items() if info["risk_level"] == "elevated"]
    if elevated:
        summary = f"Elevated genetic signals in: {', '.join(elevated)}. Proactive strategies recommended."
    else:
        summary = "No strongly elevated mental health genetic signals. Maintain wellness practices."

    return {
        "domains": domains,
        "new_markers": new_markers,
        "treatment_notes": treatment_notes,
        "resilience_factors": resilience_factors,
        "risk_factors": risk_factors[:8],
        "recommendations": recommendations,
        "summary": summary,
    }
