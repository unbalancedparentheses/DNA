"""APOE haplotype calling from rs429358 + rs7412.

Combines two SNPs into APOE epsilon type (e2/e2 through e4/e4)
with Alzheimer's disease risk context.
"""

# APOE haplotype lookup: (rs429358_genotype, rs7412_genotype) -> haplotype
# rs429358: T=e2/e3 allele, C=e4 allele
# rs7412:   T=e2 allele, C=e3/e4 allele
_HAPLOTYPE_TABLE = {
    ("TT", "TT"): "e2/e2",
    ("TT", "CT"): "e2/e3",
    ("TT", "CC"): "e3/e3",
    ("CT", "CC"): "e3/e4",
    ("CT", "CT"): "e2/e4",
    ("CC", "CC"): "e4/e4",
    # Handle reversed allele order
    ("TT", "TC"): "e2/e3",
    ("TC", "CC"): "e3/e4",
    ("TC", "CT"): "e2/e4",
    ("TC", "TC"): "e2/e4",
}

_RISK_INFO = {
    "e2/e2": {"risk_level": "reduced", "alzheimer_or": 0.6,
              "description": "Lowest Alzheimer's risk. May increase cardiovascular risk (type III hyperlipoproteinemia)."},
    "e2/e3": {"risk_level": "reduced", "alzheimer_or": 0.6,
              "description": "Below-average Alzheimer's risk. One protective e2 allele."},
    "e3/e3": {"risk_level": "average", "alzheimer_or": 1.0,
              "description": "Most common genotype (~60% of population). Average Alzheimer's risk."},
    "e3/e4": {"risk_level": "elevated", "alzheimer_or": 2.8,
              "description": "One e4 allele increases Alzheimer's risk ~2-3x. Exercise, sleep, and cardiovascular health may mitigate risk."},
    "e2/e4": {"risk_level": "moderate", "alzheimer_or": 1.2,
              "description": "Mixed: one risk allele (e4) and one protective allele (e2). Net risk near average."},
    "e4/e4": {"risk_level": "high", "alzheimer_or": 12.0,
              "description": "Highest genetic risk for late-onset Alzheimer's (~12x). Aggressive cardiovascular and lifestyle optimization recommended."},
}


def call_apoe_haplotype(genome_by_rsid):
    """Determine APOE epsilon haplotype from rs429358 and rs7412.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict with keys: apoe_type, risk_level, alzheimer_or, description,
                    confidence, details.
    """
    rs429358 = genome_by_rsid.get("rs429358", {}).get("genotype", "")
    rs7412 = genome_by_rsid.get("rs7412", {}).get("genotype", "")

    details = []
    if rs429358:
        details.append(f"rs429358: {rs429358}")
    if rs7412:
        details.append(f"rs7412: {rs7412}")

    if not rs429358 or not rs7412:
        return {
            "apoe_type": "Unknown",
            "risk_level": "unknown",
            "alzheimer_or": None,
            "description": "Insufficient data — both rs429358 and rs7412 are needed.",
            "confidence": "low",
            "details": details,
        }

    # Normalize genotype order (sort alleles)
    key = (rs429358, rs7412)
    haplotype = _HAPLOTYPE_TABLE.get(key)

    if haplotype is None:
        # Try sorting each genotype's alleles
        sorted_key = ("".join(sorted(rs429358)), "".join(sorted(rs7412)))
        haplotype = _HAPLOTYPE_TABLE.get(sorted_key)

    if haplotype is None:
        return {
            "apoe_type": "Unknown",
            "risk_level": "unknown",
            "alzheimer_or": None,
            "description": f"Unexpected genotype combination: rs429358={rs429358}, rs7412={rs7412}.",
            "confidence": "low",
            "details": details,
        }

    info = _RISK_INFO[haplotype]
    return {
        "apoe_type": haplotype,
        "risk_level": info["risk_level"],
        "alzheimer_or": info["alzheimer_or"],
        "description": info["description"],
        "confidence": "high",
        "details": details,
    }
