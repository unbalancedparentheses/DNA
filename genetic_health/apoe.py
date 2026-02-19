"""APOE haplotype calling from rs429358 + rs7412.

Combines two SNPs into APOE epsilon type (e2/e2 through e4/e4)
with Alzheimer's disease risk context.

APOE allele encoding:
  rs429358  rs7412   epsilon
  T         T        e2
  T         C        e3
  C         C        e4
  C         T        impossible (not observed in nature)

Diploid genotypes are decoded allele-by-allele.
"""


def _decode_haplotype(rs429358_genotype, rs7412_genotype):
    """Deterministic APOE haplotype from rs429358 + rs7412 genotypes.

    Decodes each allele pair into an epsilon allele, then sorts
    for canonical representation.
    """
    if len(rs429358_genotype) != 2 or len(rs7412_genotype) != 2:
        return None

    # Map each allele position to an epsilon allele
    epsilon_alleles = []
    for i in range(2):
        a358 = rs429358_genotype[i]
        a7412 = rs7412_genotype[i]
        if a358 == "T" and a7412 == "T":
            epsilon_alleles.append("e2")
        elif a358 == "T" and a7412 == "C":
            epsilon_alleles.append("e3")
        elif a358 == "C" and a7412 == "C":
            epsilon_alleles.append("e4")
        else:
            # C+T is not a natural APOE haplotype; try other phasing
            return None

    if len(epsilon_alleles) == 2:
        epsilon_alleles.sort()
        return f"{epsilon_alleles[0]}/{epsilon_alleles[1]}"
    return None


def _try_all_phasings(rs429358, rs7412):
    """Try all possible allele phasings to find a valid APOE haplotype."""
    # 23andMe genotypes are unphased, so TC could be T|C or C|T
    alleles_358 = [rs429358[0] + rs429358[1], rs429358[1] + rs429358[0]]
    alleles_7412 = [rs7412[0] + rs7412[1], rs7412[1] + rs7412[0]]

    for a358 in alleles_358:
        for a7412 in alleles_7412:
            result = _decode_haplotype(a358, a7412)
            if result:
                return result
    return None

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

    # Deterministic decoding — try all allele phasings
    haplotype = _try_all_phasings(rs429358, rs7412)

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
