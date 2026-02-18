"""CPIC-style star-allele calling for key pharmacogenes.

Calls diplotypes and maps to metabolizer phenotypes for CYP2C19, CYP2C9,
and CYP2D6 based on defining SNPs available in 23andMe data.
"""

# Star allele definitions: {gene: {star_allele: {rsid: variant_allele}}}
# *1 is always the reference (wild-type) — no defining variants.

STAR_ALLELE_DEFINITIONS = {
    "CYP2C19": {
        "function_map": {
            "*1": "normal",
            "*2": "no_function",
            "*3": "no_function",
            "*17": "increased",
        },
        "alleles": {
            "*2": {"rs4244285": "A"},
            "*3": {"rs4986893": "A"},
            "*17": {"rs12248560": "T"},
        },
        "snps": ["rs4244285", "rs4986893", "rs12248560"],
    },
    "CYP2C9": {
        "function_map": {
            "*1": "normal",
            "*2": "decreased",
            "*3": "decreased",
        },
        "alleles": {
            "*2": {"rs1799853": "T"},
            "*3": {"rs1057910": "C"},
        },
        "snps": ["rs1799853", "rs1057910"],
    },
    "CYP2D6": {
        "function_map": {
            "*1": "normal",
            "*4": "no_function",
            "*10": "decreased",
        },
        "alleles": {
            "*4": {"rs3892097": "A"},
            "*10": {"rs1065852": "T"},
        },
        "snps": ["rs3892097", "rs1065852"],
    },
}

# Phenotype mapping from diplotype function pairs
# (sorted pair of functions) -> metabolizer phenotype
_PHENOTYPE_MAP = {
    ("increased", "increased"): "ultrarapid",
    ("increased", "normal"): "rapid",
    ("normal", "normal"): "normal",
    ("decreased", "normal"): "intermediate",
    ("decreased", "increased"): "normal",  # one up, one down ≈ normal
    ("decreased", "decreased"): "intermediate",
    ("increased", "no_function"): "intermediate",
    ("no_function", "normal"): "intermediate",
    ("decreased", "no_function"): "poor",
    ("no_function", "no_function"): "poor",
}


def _call_gene(gene, definitions, genome_by_rsid):
    """Call star alleles for a single gene.

    Returns dict with gene, diplotype, phenotype, snps_found, clinical_note.
    """
    alleles_def = definitions["alleles"]
    function_map = definitions["function_map"]
    all_snps = definitions["snps"]

    # Check which SNPs are present
    snps_found = []
    snps_missing = []
    for rsid in all_snps:
        if rsid in genome_by_rsid and genome_by_rsid[rsid].get("genotype"):
            snps_found.append(rsid)
        else:
            snps_missing.append(rsid)

    if not snps_found:
        return {
            "gene": gene,
            "diplotype": "Unknown",
            "phenotype": "Unknown",
            "snps_found": 0,
            "snps_total": len(all_snps),
            "clinical_note": f"No defining SNPs found for {gene} in genome data.",
        }

    # Count variant alleles for each star allele on each chromosome copy
    # We identify alleles present in the genotype
    called_alleles = []  # list of (star_allele, count_of_copies)

    for star_name, defining_snps in alleles_def.items():
        # For each defining SNP, check if variant allele is present
        # Since most star alleles are defined by a single SNP in 23andMe,
        # we check allele count directly
        total_copies = 0
        all_defined = True

        for rsid, variant_allele in defining_snps.items():
            entry = genome_by_rsid.get(rsid)
            if entry is None or not entry.get("genotype"):
                all_defined = False
                continue

            genotype = entry["genotype"]
            copies = sum(1 for a in genotype if a == variant_allele)
            total_copies += copies

        if all_defined or total_copies > 0:
            called_alleles.append((star_name, total_copies))

    # Build diplotype: two allele calls
    # Start with *1/*1, then replace based on variant alleles found
    allele1 = "*1"
    allele2 = "*1"

    for star_name, copies in sorted(called_alleles, key=lambda x: -x[1]):
        if copies >= 2:
            allele1 = star_name
            allele2 = star_name
        elif copies == 1:
            if allele1 == "*1":
                allele1 = star_name
            elif allele2 == "*1":
                allele2 = star_name

    # Sort alleles for consistent representation
    diplotype_alleles = sorted([allele1, allele2])
    diplotype = f"{diplotype_alleles[0]}/{diplotype_alleles[1]}"

    # Map to phenotype
    func1 = function_map.get(diplotype_alleles[0], "normal")
    func2 = function_map.get(diplotype_alleles[1], "normal")
    func_pair = tuple(sorted([func1, func2]))
    phenotype = _PHENOTYPE_MAP.get(func_pair)
    if phenotype is None:
        phenotype = "Indeterminate"
        notes.append(
            f"Unexpected function combination ({func1} + {func2}). "
            "Clinical interpretation needed."
        )

    # Clinical notes
    notes = []
    if gene == "CYP2D6":
        notes.append(
            "23andMe cannot detect CYP2D6 gene deletions (*5) or "
            "duplications (*1xN). This result may be incomplete."
        )
    if snps_missing:
        notes.append(
            f"Missing SNPs: {', '.join(snps_missing)}. "
            "Result based on available data only."
        )

    clinical_note = " ".join(notes) if notes else "All defining SNPs found."

    return {
        "gene": gene,
        "diplotype": diplotype,
        "phenotype": phenotype,
        "snps_found": len(snps_found),
        "snps_total": len(all_snps),
        "clinical_note": clinical_note,
    }


def call_star_alleles(genome_by_rsid):
    """Call star alleles for all defined pharmacogenes.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict mapping gene name to result dict with keys:
        gene, diplotype, phenotype, snps_found, snps_total, clinical_note.
    """
    results = {}
    for gene, definitions in STAR_ALLELE_DEFINITIONS.items():
        results[gene] = _call_gene(gene, definitions, genome_by_rsid)
    return results
