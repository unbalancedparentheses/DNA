"""Estimate maternal mitochondrial haplogroup from 23andMe MT SNP data."""

# Decision tree of mitochondrial haplogroup-defining SNPs.
# Each entry: rsID -> {allele, haplogroup, description, parent}
# Walk from root (most specific match wins).

MT_HAPLOGROUP_TREE = [
    # --- African L lineages ---
    {"rsid": "rs2853499", "allele": "G", "haplogroup": "L0",
     "description": "Ancient African maternal lineage (Khoisan, Southern African)"},
    {"rsid": "rs3088309", "allele": "G", "haplogroup": "L1",
     "description": "West/Central African maternal lineage"},
    {"rsid": "rs28358571", "allele": "T", "haplogroup": "L2",
     "description": "West African maternal lineage (~25% of sub-Saharan Africans)"},
    {"rsid": "rs2854128", "allele": "T", "haplogroup": "L3",
     "description": "Out-of-Africa ancestral lineage (~70,000 years ago)"},

    # --- Out-of-Africa: M and N branches ---
    {"rsid": "rs28358579", "allele": "C", "haplogroup": "M",
     "description": "Southern route out of Africa (common in South/East Asia, Oceania)"},
    {"rsid": "rs28358587", "allele": "T", "haplogroup": "N",
     "description": "Northern route out of Africa (ancestral to most European/Asian haplogroups)"},

    # --- European haplogroups (all derive from N via R) ---
    {"rsid": "rs2032658", "allele": "G", "haplogroup": "H",
     "description": "Most common European haplogroup (~40% of Europeans)"},
    {"rsid": "rs2853825", "allele": "T", "haplogroup": "H1",
     "description": "H1 subclade (~25% of Europeans, common in Iberia and W. Europe)"},
    {"rsid": "rs28359168", "allele": "A", "haplogroup": "H2",
     "description": "H2 subclade (associated with Neolithic expansion in Europe)"},
    {"rsid": "rs3928306", "allele": "A", "haplogroup": "J",
     "description": "European/Near Eastern lineage (~12% Europeans, associated with farming expansion)"},
    {"rsid": "rs28358270", "allele": "T", "haplogroup": "J1",
     "description": "J1 subclade (Neolithic expansion from Near East)"},
    {"rsid": "rs2853826", "allele": "T", "haplogroup": "T",
     "description": "European lineage (~9% of Europeans, sister clade to J)"},
    {"rsid": "rs41456348", "allele": "C", "haplogroup": "T2",
     "description": "T2 subclade (most common T subclade in Europe)"},
    {"rsid": "rs2853518", "allele": "T", "haplogroup": "K",
     "description": "European lineage (~10% of Europeans, Ashkenazi Jewish ~32%)"},
    {"rsid": "rs3135028", "allele": "G", "haplogroup": "U",
     "description": "Ancient European/West Asian lineage (~11% of Europeans)"},
    {"rsid": "rs28358280", "allele": "T", "haplogroup": "U5",
     "description": "U5: oldest European-specific lineage (~10,000 years, Mesolithic hunter-gatherers)"},
    {"rsid": "rs35219750", "allele": "C", "haplogroup": "V",
     "description": "European lineage (~4% of Europeans, common in Scandinavia and Iberia)"},
    {"rsid": "rs2853515", "allele": "G", "haplogroup": "W",
     "description": "European/West Asian lineage (~3% of Europeans)"},
    {"rsid": "rs3134015", "allele": "A", "haplogroup": "X",
     "description": "Found in Europe, W. Asia, and Native Americans (~2% of Europeans)"},

    # --- East Asian and Native American haplogroups (derive from M or N) ---
    {"rsid": "rs28977561", "allele": "A", "haplogroup": "A",
     "description": "East Asian and Native American maternal lineage"},
    {"rsid": "rs8896", "allele": "T", "haplogroup": "B",
     "description": "East Asian, Polynesian, and Native American maternal lineage"},
    {"rsid": "rs41323649", "allele": "C", "haplogroup": "C",
     "description": "East Asian, Central Asian, and Native American maternal lineage"},
    {"rsid": "rs28357968", "allele": "T", "haplogroup": "D",
     "description": "East Asian and some Native American maternal lineage"},
    {"rsid": "rs41419549", "allele": "G", "haplogroup": "G",
     "description": "East Asian maternal lineage (common in northern East Asia)"},

    # --- South/West Asian ---
    {"rsid": "rs28359172", "allele": "G", "haplogroup": "R",
     "description": "South and West Asian / European macro-haplogroup"},
]

# Lineage region mapping
_LINEAGE_MAP = {
    "L0": "African", "L1": "African", "L2": "African", "L3": "African",
    "M": "South/East Asian", "N": "Eurasian",
    "H": "European", "H1": "European", "H2": "European",
    "J": "European/Near Eastern", "J1": "European/Near Eastern",
    "T": "European", "T2": "European",
    "K": "European", "U": "European/West Asian", "U5": "European",
    "V": "European", "W": "European/West Asian", "X": "Eurasian/Native American",
    "A": "East Asian/Native American", "B": "East Asian/Polynesian/Native American",
    "C": "East Asian/Native American", "D": "East Asian/Native American",
    "G": "East Asian", "R": "South/West Asian",
}


def estimate_mt_haplogroup(genome_by_rsid):
    """Estimate maternal mitochondrial haplogroup from genome data.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict with keys: haplogroup, description, confidence, markers_found,
        markers_tested, lineage, details.
    """
    matches = []
    markers_found = 0
    markers_tested = len(MT_HAPLOGROUP_TREE)

    for marker in MT_HAPLOGROUP_TREE:
        rsid = marker["rsid"]
        entry = genome_by_rsid.get(rsid)
        if entry is None:
            continue

        genotype = entry.get("genotype", "")
        if not genotype:
            continue

        markers_found += 1
        target_allele = marker["allele"]

        # MT DNA is haploid: genotype is typically a single character,
        # but 23andMe may report it as two identical characters
        if target_allele in genotype:
            matches.append(marker)

    # Most specific match wins (later entries in tree are more specific)
    if matches:
        best = matches[-1]
        haplogroup = best["haplogroup"]
        description = best["description"]
    else:
        haplogroup = "Unknown"
        description = "No defining mitochondrial SNPs matched"

    # Confidence based on markers found
    if markers_found >= 15:
        confidence = "high"
    elif markers_found >= 5:
        confidence = "moderate"
    elif markers_found >= 1:
        confidence = "low"
    else:
        confidence = "none"

    lineage = _LINEAGE_MAP.get(haplogroup, "Unknown") + " maternal"

    return {
        "haplogroup": haplogroup,
        "description": description,
        "confidence": confidence,
        "markers_found": markers_found,
        "markers_tested": markers_tested,
        "lineage": lineage,
        "details": [
            {
                "rsid": m["rsid"],
                "haplogroup": m["haplogroup"],
                "allele": m["allele"],
            }
            for m in matches
        ],
    }
