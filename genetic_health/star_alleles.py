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
            "*4": "no_function",
            "*5": "no_function",
            "*6": "no_function",
            "*7": "no_function",
            "*8": "no_function",
            "*17": "increased",
        },
        "alleles": {
            "*2": {"rs4244285": "A"},
            "*3": {"rs4986893": "A"},
            "*4": {"rs28399504": "A"},
            "*5": {"rs56337013": "T"},
            "*6": {"rs72552267": "A"},
            "*7": {"rs72558186": "T"},
            "*8": {"rs41291556": "T"},
            "*17": {"rs12248560": "T"},
        },
        "snps": ["rs4244285", "rs4986893", "rs28399504", "rs56337013",
                 "rs72552267", "rs72558186", "rs41291556", "rs12248560"],
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
            "*2": "normal",
            "*4": "no_function",
            "*5": "no_function",
            "*6": "no_function",
            "*9": "decreased",
            "*10": "decreased",
            "*41": "decreased",
        },
        "alleles": {
            "*2": {"rs16947": "A"},
            "*4": {"rs3892097": "A"},
            "*6": {"rs5030655": "T"},
            "*9": {"rs5030656": "A"},
            "*10": {"rs1065852": "T"},
            "*41": {"rs28371725": "T"},
        },
        "snps": ["rs16947", "rs3892097", "rs5030655", "rs5030656",
                 "rs1065852", "rs28371725"],
        "clinical_note": (
            "23andMe cannot detect CYP2D6 gene deletions (*5) or "
            "duplications (*1xN/*2xN). Copy number variants account "
            "for ~20-30% of poor/ultrarapid metabolizers. This result "
            "may underestimate phenotype severity."
        ),
    },
    "DPYD": {
        "function_map": {
            "*1": "normal",
            "*2A": "no_function",
            "*13": "decreased",
        },
        "alleles": {
            "*2A": {"rs3918290": "A"},
            "*13": {"rs55886062": "A"},
        },
        "snps": ["rs3918290", "rs55886062"],
        "clinical_note": "DPYD deficiency can cause fatal 5-FU/capecitabine toxicity. Pre-treatment testing recommended.",
    },
    "TPMT": {
        "function_map": {
            "*1": "normal",
            "*2": "no_function",
            "*3A": "no_function",
            "*3B": "no_function",
            "*3C": "decreased",
            "*4": "no_function",
        },
        "alleles": {
            "*2": {"rs1800462": "G"},
            "*3A": {"rs1800460": "A", "rs1142345": "C"},
            "*3B": {"rs1800460": "A"},
            "*3C": {"rs1142345": "C"},
            "*4": {"rs1800584": "A"},
        },
        "snps": ["rs1800462", "rs1800460", "rs1142345", "rs1800584"],
        "clinical_note": "TPMT deficiency can cause life-threatening myelosuppression with azathioprine/6-MP.",
    },
    "UGT1A1": {
        "function_map": {
            "*1": "normal",
            "*28": "decreased",
            "*6": "decreased",
        },
        "alleles": {
            "*28": {"rs8175347": "A"},
            "*6": {"rs4148323": "A"},
        },
        "snps": ["rs8175347", "rs4148323"],
        "clinical_note": "UGT1A1*28 causes Gilbert syndrome and irinotecan toxicity. UGT1A1*6 common in East Asian populations.",
    },
}

# Phenotype mapping from diplotype function pairs
# (sorted pair of functions) -> metabolizer phenotype
# All 10 combinations of {normal, increased, decreased, no_function} covered.
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
# Verify completeness at import time
_FUNCTIONS = ("normal", "increased", "decreased", "no_function")
for _f1 in _FUNCTIONS:
    for _f2 in _FUNCTIONS:
        _pair = tuple(sorted((_f1, _f2)))
        assert _pair in _PHENOTYPE_MAP, f"Missing phenotype mapping for {_pair}"


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
            "coverage": 0.0,
            "confidence": "low",
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

    # Clinical notes
    notes = []

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

    # Gene-specific clinical note from definitions
    gene_note = definitions.get("clinical_note")
    if gene_note:
        notes.append(gene_note)
    if snps_missing:
        notes.append(
            f"Missing SNPs: {', '.join(snps_missing)}. "
            "Result based on available data only."
        )

    clinical_note = " ".join(notes) if notes else "All defining SNPs found."

    coverage = len(snps_found) / len(all_snps) if all_snps else 0
    confidence = "high" if coverage >= 0.8 else "moderate" if coverage >= 0.5 else "low"

    return {
        "gene": gene,
        "diplotype": diplotype,
        "phenotype": phenotype,
        "snps_found": len(snps_found),
        "snps_total": len(all_snps),
        "coverage": round(coverage, 2),
        "confidence": confidence,
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
