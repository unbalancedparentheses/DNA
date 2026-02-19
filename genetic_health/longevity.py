"""Longevity and healthspan profiling.

Synthesizes longevity-associated genes, disease PRS, and lifestyle
variant data into a consolidated aging/healthspan profile.
"""


# Longevity-associated SNPs with published evidence
LONGEVITY_SNPS = {
    "rs2802292": {
        "gene": "FOXO3",
        "name": "Forkhead box O3 — autophagy & longevity",
        "longevity_allele": "G",
        "effect": "protective",
        "or_longevity": 1.17,
        "reference": "Willcox et al. 2008 (PNAS)",
    },
    "rs5882": {
        "gene": "CETP",
        "name": "Cholesteryl ester transfer protein — HDL & longevity",
        "longevity_allele": "G",
        "effect": "protective",
        "or_longevity": 1.14,
        "reference": "Barzilai et al. 2003 (JAMA)",
    },
    "rs429358": {
        "gene": "APOE",
        "name": "APOE e4 allele — Alzheimer's & cardiovascular risk",
        "longevity_allele": "T",  # T = e2/e3, C = e4
        "effect": "risk_if_C",
        "or_longevity": 0.73,  # e4 carriers live shorter on average
        "reference": "Deelen et al. 2011 (Aging Cell)",
    },
    "rs1042522": {
        "gene": "TP53",
        "name": "Tumor protein p53 — cancer surveillance",
        "longevity_allele": "G",
        "effect": "context_dependent",
        "or_longevity": 1.0,
        "reference": "Bojesen & Nordestgaard 2008 (Cell)",
    },
    "rs4880": {
        "gene": "SOD2",
        "name": "Superoxide dismutase — oxidative stress defense",
        "longevity_allele": "T",  # Ala variant — targeted to mitochondria
        "effect": "protective",
        "or_longevity": 1.10,
        "reference": "Soerensen et al. 2012 (Mech Ageing Dev)",
    },
    "rs2536": {
        "gene": "MTOR",
        "name": "Mechanistic target of rapamycin — growth/aging pathway",
        "longevity_allele": "C",
        "effect": "protective",
        "or_longevity": 1.08,
        "reference": "Johnson et al. 2013 (Nature)",
    },
    "rs7069102": {
        "gene": "SIRT1",
        "name": "Sirtuin 1 — NAD+ deacetylase, caloric restriction pathway",
        "longevity_allele": "G",
        "effect": "protective",
        "or_longevity": 1.12,
        "reference": "Figarska et al. 2013 (Age)",
    },
    "rs9536314": {
        "gene": "KLOTHO",
        "name": "Klotho — anti-aging hormone, kidney & brain protection",
        "longevity_allele": "T",
        "effect": "protective",
        "or_longevity": 1.15,
        "reference": "Arking et al. 2002 (PNAS)",
    },
    "rs2229765": {
        "gene": "IGF1R",
        "name": "IGF-1 receptor — growth hormone/insulin pathway",
        "longevity_allele": "A",
        "effect": "protective",
        "or_longevity": 1.09,
        "reference": "Suh et al. 2008 (PNAS)",
    },
    "rs10936599": {
        "gene": "TERC",
        "name": "Telomerase RNA component — telomere maintenance",
        "longevity_allele": "C",
        "effect": "protective",
        "or_longevity": 1.06,
        "reference": "Codd et al. 2013 (Nature Genetics)",
    },
    "rs2736100": {
        "gene": "TERT",
        "name": "Telomerase reverse transcriptase — telomere length",
        "longevity_allele": "C",
        "effect": "protective",
        "or_longevity": 1.07,
        "reference": "Codd et al. 2013 (Nature Genetics)",
    },
}

# Healthspan domains: map lifestyle gene statuses to domain scores
_DOMAIN_GENES = {
    "cardiovascular": {
        "genes": ["APOE", "CETP", "AGT", "ACE", "AGTR1", "F5", "F2", "MTHFR"],
        "protective_statuses": {"protective", "longevity", "optimal", "normal", "reference"},
        "risk_statuses": {"e3_e4", "e4_e4", "elevated", "high", "high_risk",
                          "heterozygous", "homozygous", "reduced", "severely_reduced"},
    },
    "metabolic": {
        "genes": ["FTO", "TCF7L2", "PPARG", "PPARA", "PPARGC1A", "MTNR1B"],
        "protective_statuses": {"normal", "reference", "protective"},
        "risk_statuses": {"risk", "high_risk", "slow", "reduced"},
    },
    "neurological": {
        "genes": ["APOE", "BDNF", "COMT", "ANKK1"],
        "protective_statuses": {"normal", "fast", "optimal", "protective"},
        "risk_statuses": {"e3_e4", "e4_e4", "reduced", "met_carrier", "slow", "risk"},
    },
    "inflammation": {
        "genes": ["TNF", "IL6", "SOD2", "GSTP1"],
        "protective_statuses": {"normal", "reference", "low", "protective"},
        "risk_statuses": {"high", "elevated", "reduced", "slow"},
    },
    "cancer_defense": {
        "genes": ["TP53", "SOD2", "NAT2", "GSTP1"],
        "protective_statuses": {"normal", "fast", "protective"},
        "risk_statuses": {"slow", "reduced", "risk"},
    },
}


def profile_longevity(genome_by_rsid, lifestyle_findings=None, apoe=None,
                      prs_results=None):
    """Build a comprehensive longevity and healthspan profile.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict.
    lifestyle_findings : list or None
        Findings from analyze_lifestyle_health().
    apoe : dict or None
        APOE haplotype result.
    prs_results : dict or None
        PRS results dict.

    Returns
    -------
    dict with keys: longevity_score, longevity_alleles, healthspan_domains,
                    top_risks, top_protective, interventions, summary.
    """
    if lifestyle_findings is None:
        lifestyle_findings = []

    # --- Score longevity alleles ---
    longevity_alleles = []
    protective_count = 0
    risk_count = 0
    total_checked = 0

    for rsid, info in LONGEVITY_SNPS.items():
        entry = genome_by_rsid.get(rsid)
        if not entry:
            continue
        genotype = entry.get("genotype", "")
        if not genotype:
            continue

        total_checked += 1
        allele = info["longevity_allele"]
        copies = sum(1 for a in genotype if a == allele)

        if copies == 2:
            status = "homozygous_protective"
            protective_count += 2
        elif copies == 1:
            status = "heterozygous"
            protective_count += 1
        else:
            status = "absent"
            risk_count += 1

        longevity_alleles.append({
            "rsid": rsid,
            "gene": info["gene"],
            "name": info["name"],
            "genotype": genotype,
            "longevity_allele": allele,
            "copies": copies,
            "status": status,
            "reference": info["reference"],
        })

    # Longevity score: 0-100 scale based on protective allele fraction
    max_possible = total_checked * 2 if total_checked > 0 else 1
    raw_score = protective_count / max_possible
    longevity_score = round(raw_score * 100, 1)

    # --- Healthspan domain scoring ---
    gene_status_map = {}
    for f in lifestyle_findings:
        gene = f.get("gene", "")
        status = f.get("status", "")
        if gene and status:
            gene_status_map[gene] = status

    # Factor in APOE
    if apoe and apoe.get("apoe_type") != "Unknown":
        apoe_type = apoe["apoe_type"]
        if "e4" in apoe_type:
            gene_status_map.setdefault("APOE", apoe_type.replace("/", "_"))

    healthspan_domains = {}
    for domain, config in _DOMAIN_GENES.items():
        domain_score = 50  # neutral baseline
        genes_found = 0
        for gene in config["genes"]:
            status = gene_status_map.get(gene, "")
            if not status:
                continue
            genes_found += 1
            if status in config["protective_statuses"]:
                domain_score += 8
            elif status in config["risk_statuses"]:
                domain_score -= 10
        domain_score = max(10, min(90, domain_score))
        healthspan_domains[domain] = {
            "score": domain_score,
            "genes_found": genes_found,
            "rating": "good" if domain_score >= 60 else "average" if domain_score >= 40 else "attention",
        }

    # --- Top risks and protective factors ---
    top_risks = []
    top_protective = []

    if prs_results:
        for cid, r in prs_results.items():
            cat = r.get("risk_category", "average")
            if cat in ("elevated", "high"):
                top_risks.append(f"Elevated PRS for {r['name']} ({r['percentile']:.0f}th percentile)")
            elif cat == "low":
                top_protective.append(f"Low PRS for {r['name']} ({r['percentile']:.0f}th percentile)")

    if apoe:
        risk = apoe.get("risk_level", "")
        if risk in ("elevated", "high"):
            top_risks.append(f"APOE {apoe['apoe_type']}: {apoe['description']}")
        elif risk == "reduced":
            top_protective.append(f"APOE {apoe['apoe_type']}: Protective for Alzheimer's")

    for a in longevity_alleles:
        if a["copies"] == 2:
            top_protective.append(f"{a['gene']}: Homozygous for longevity allele")
        elif a["copies"] == 0:
            top_risks.append(f"{a['gene']}: No longevity allele copies")

    # --- Interventions ranked by genetic support ---
    interventions = []
    if gene_status_map.get("PPARGC1A") or gene_status_map.get("BDNF"):
        interventions.append({
            "intervention": "Regular aerobic exercise (150+ min/week)",
            "genetic_support": "high",
            "why": "Strongest BDNF booster + mitochondrial biogenesis via PGC-1α",
        })
    interventions.append({
        "intervention": "Mediterranean dietary pattern",
        "genetic_support": "high",
        "why": "Reduces inflammation (TNF/IL6), supports APOE-linked cardiovascular health",
    })
    if any(a["gene"] == "FOXO3" and a["copies"] > 0 for a in longevity_alleles):
        interventions.append({
            "intervention": "Time-restricted eating (12-16hr overnight fast)",
            "genetic_support": "moderate",
            "why": "FOXO3 protective allele activates autophagy pathways enhanced by fasting",
        })
    if any(a["gene"] == "SIRT1" for a in longevity_alleles):
        interventions.append({
            "intervention": "Resveratrol / NAD+ precursors (NMN, NR)",
            "genetic_support": "moderate",
            "why": "SIRT1 pathway responds to NAD+ boosting compounds",
        })
    interventions.append({
        "intervention": "Sleep optimization (7-9 hrs, consistent schedule)",
        "genetic_support": "high",
        "why": "Telomere maintenance (TERC/TERT) requires adequate sleep",
    })
    interventions.append({
        "intervention": "Stress management (meditation, social connection)",
        "genetic_support": "moderate",
        "why": "Reduces cortisol-driven telomere shortening and inflammation",
    })

    # Summary
    if longevity_score >= 65:
        summary = "Your longevity genetic profile is favorable."
    elif longevity_score >= 40:
        summary = "Your longevity profile is average — lifestyle optimization has high impact."
    else:
        summary = "Your longevity profile shows some risk factors — proactive interventions recommended."

    return {
        "longevity_score": longevity_score,
        "longevity_alleles": longevity_alleles,
        "alleles_checked": total_checked,
        "healthspan_domains": healthspan_domains,
        "top_risks": top_risks[:8],
        "top_protective": top_protective[:8],
        "interventions": interventions,
        "summary": summary,
    }
