"""Sleep and circadian profiling from chronotype-associated SNPs.

Combines clock gene variants to predict chronotype (morning/evening),
sleep architecture tendencies, and caffeine timing recommendations.
"""


# Chronotype-associated SNPs
CHRONOTYPE_SNPS = {
    "rs1801260": {
        "gene": "CLOCK",
        "trait": "chronotype",
        "evening_allele": "C",
        "weight": 1.0,
        "description": "CLOCK 3111T>C — evening preference",
        "reference": "Katzenberg et al. 1998 (Sleep)",
    },
    "rs2304672": {
        "gene": "CLOCK",
        "trait": "chronotype",
        "evening_allele": "C",
        "weight": 0.7,
        "description": "CLOCK — additional chronotype variant",
        "reference": "Mishima et al. 2005 (BMC Med Genet)",
    },
    "rs934945": {
        "gene": "PER2",
        "trait": "chronotype",
        "evening_allele": "G",
        "weight": 0.9,
        "description": "PER2 — period circadian protein 2",
        "reference": "Carpen et al. 2006 (Sleep)",
    },
    "rs228697": {
        "gene": "PER3",
        "trait": "chronotype",
        "evening_allele": "C",
        "weight": 0.8,
        "description": "PER3 — short sleep tolerance, evening preference",
        "reference": "Archer et al. 2003 (Sleep)",
    },
    "rs2278749": {
        "gene": "ARNTL",
        "trait": "chronotype",
        "evening_allele": "T",
        "weight": 0.6,
        "description": "BMAL1/ARNTL — core circadian oscillator",
        "reference": "Woon et al. 2007 (Chronobiol Int)",
    },
    "rs2305160": {
        "gene": "NPAS2",
        "trait": "chronotype",
        "evening_allele": "G",
        "weight": 0.7,
        "description": "NPAS2 — neuronal PAS domain protein 2",
        "reference": "Johansson et al. 2003 (J Biol Rhythms)",
    },
    "rs10830963": {
        "gene": "MTNR1B",
        "trait": "melatonin",
        "evening_allele": "G",
        "weight": 0.8,
        "description": "Melatonin receptor 1B — delayed melatonin onset",
        "reference": "Bonnefond et al. 2012 (Nature Genetics)",
    },
    "rs6295": {
        "gene": "HTR1A",
        "trait": "sleep_quality",
        "evening_allele": "G",
        "weight": 0.5,
        "description": "Serotonin 1A receptor — REM sleep regulation",
        "reference": "Huang et al. 2004 (Neuropsychopharmacology)",
    },
    "rs73598374": {
        "gene": "ADA",
        "trait": "sleep_depth",
        "evening_allele": "T",
        "weight": 0.7,
        "description": "Adenosine deaminase — deep sleep pressure",
        "reference": "Rétey et al. 2005 (PNAS)",
    },
}


def profile_sleep(genome_by_rsid, lifestyle_findings=None):
    """Build a sleep and circadian profile from genome data.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict.
    lifestyle_findings : list or None
        Lifestyle findings for caffeine/adenosine context.

    Returns
    -------
    dict with keys: chronotype, chronotype_score, sleep_markers,
                    optimal_sleep_window, caffeine_cutoff,
                    recommendations, confidence.
    """
    if lifestyle_findings is None:
        lifestyle_findings = []

    evening_score = 0.0
    total_weight = 0.0
    markers_found = 0
    sleep_markers = []

    for rsid, info in CHRONOTYPE_SNPS.items():
        entry = genome_by_rsid.get(rsid)
        if not entry:
            continue
        genotype = entry.get("genotype", "")
        if not genotype:
            continue

        markers_found += 1
        allele = info["evening_allele"]
        copies = sum(1 for a in genotype if a == allele)
        contribution = copies * info["weight"]
        evening_score += contribution
        total_weight += 2 * info["weight"]  # max possible for this SNP

        sleep_markers.append({
            "rsid": rsid,
            "gene": info["gene"],
            "trait": info["trait"],
            "genotype": genotype,
            "evening_allele_copies": copies,
            "description": info["description"],
        })

    # Normalize to 0-100 scale (0=extreme morning, 100=extreme evening)
    if total_weight > 0:
        chronotype_score = round((evening_score / total_weight) * 100, 1)
    else:
        chronotype_score = 50.0

    # Classify chronotype
    if chronotype_score >= 70:
        chronotype = "Definite Evening (Night Owl)"
        optimal_window = "12:00 AM – 8:00 AM"
        caffeine_cutoff = "12:00 PM (noon)"
        peak_alertness = "10:00 AM – 2:00 PM and 7:00 PM – 11:00 PM"
    elif chronotype_score >= 55:
        chronotype = "Moderate Evening"
        optimal_window = "11:30 PM – 7:30 AM"
        caffeine_cutoff = "1:00 PM"
        peak_alertness = "10:00 AM – 1:00 PM and 6:00 PM – 10:00 PM"
    elif chronotype_score >= 45:
        chronotype = "Intermediate (Neither)"
        optimal_window = "11:00 PM – 7:00 AM"
        caffeine_cutoff = "2:00 PM"
        peak_alertness = "9:00 AM – 12:00 PM and 3:00 PM – 7:00 PM"
    elif chronotype_score >= 30:
        chronotype = "Moderate Morning"
        optimal_window = "10:00 PM – 6:00 AM"
        caffeine_cutoff = "2:00 PM"
        peak_alertness = "8:00 AM – 12:00 PM"
    else:
        chronotype = "Definite Morning (Early Bird)"
        optimal_window = "9:30 PM – 5:30 AM"
        caffeine_cutoff = "12:00 PM (noon)"
        peak_alertness = "6:00 AM – 11:00 AM"

    # Check for caffeine sensitivity from lifestyle findings
    caffeine_sensitive = False
    for f in lifestyle_findings:
        if f.get("gene") == "CYP1A2" and f.get("status") in ("slow", "intermediate"):
            caffeine_sensitive = True
        if f.get("gene") == "ADORA2A" and f.get("status") == "anxiety_prone":
            caffeine_sensitive = True

    if caffeine_sensitive:
        caffeine_cutoff = "10:00 AM or avoid entirely"

    # Confidence
    if markers_found >= 6:
        confidence = "high"
    elif markers_found >= 3:
        confidence = "moderate"
    else:
        confidence = "low"

    # Recommendations
    recommendations = []
    recommendations.append(f"Target sleep window: {optimal_window}")
    recommendations.append(f"Last caffeine by: {caffeine_cutoff}")
    recommendations.append(f"Peak alertness hours: {peak_alertness}")

    if chronotype_score >= 60:
        recommendations.append(
            "Bright light exposure within 30 min of waking helps shift circadian rhythm earlier"
        )
        recommendations.append(
            "Blue light filter (f.lux/Night Shift) starting 2 hours before target bedtime"
        )
        recommendations.append(
            "Melatonin (0.3-1mg) 2-3 hours before target bedtime may help if shifting schedule"
        )
    elif chronotype_score <= 40:
        recommendations.append(
            "Avoid bright light in the evening to preserve early chronotype"
        )
        recommendations.append(
            "Early morning exercise reinforces morning circadian rhythm"
        )

    if caffeine_sensitive:
        recommendations.append(
            "Your caffeine metabolism genes suggest high sensitivity — "
            "consider switching to green tea or eliminating caffeine"
        )

    # Check for deep sleep markers
    deep_sleep_note = ""
    for m in sleep_markers:
        if m["gene"] == "ADA" and m["evening_allele_copies"] > 0:
            deep_sleep_note = (
                "ADA variant detected: you likely build sleep pressure faster, "
                "meaning better deep sleep but higher sensitivity to sleep deprivation."
            )

    return {
        "chronotype": chronotype,
        "chronotype_score": chronotype_score,
        "sleep_markers": sleep_markers,
        "markers_found": markers_found,
        "optimal_sleep_window": optimal_window,
        "caffeine_cutoff": caffeine_cutoff,
        "peak_alertness": peak_alertness,
        "caffeine_sensitive": caffeine_sensitive,
        "deep_sleep_note": deep_sleep_note,
        "recommendations": recommendations,
        "confidence": confidence,
    }
