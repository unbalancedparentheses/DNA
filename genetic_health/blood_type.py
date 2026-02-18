"""Predict ABO blood group and Rh factor from 23andMe SNP data."""


def predict_blood_type(genome_by_rsid):
    """Predict blood type from ABO and Rh proxy SNPs.

    Uses:
    - rs505922 (C/T): proxy for ABO O allele (T = O)
    - rs8176746 (C/T): B antigen (T = B)
    - rs590787 (C/T): proxy for RhD status (T = Rh-)

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.

    Returns
    -------
    dict with keys: blood_type, abo, rh, confidence, details.
    """
    details = []

    # --- ABO determination ---
    abo_proxy = genome_by_rsid.get("rs505922", {}).get("genotype", "")
    b_allele = genome_by_rsid.get("rs8176746", {}).get("genotype", "")

    abo = None
    abo_confidence = 0

    if abo_proxy:
        details.append(f"rs505922 (ABO proxy): {abo_proxy}")
        abo_confidence += 1
    if b_allele:
        details.append(f"rs8176746 (B allele): {b_allele}")
        abo_confidence += 1

    # Determine ABO from proxy SNPs
    has_b = b_allele and "T" in b_allele

    if abo_proxy == "TT":
        # Strong O signal
        if has_b:
            abo = "B"  # Rare: B allele + O proxy — B wins
        else:
            abo = "O"
    elif abo_proxy in ("CT", "TC"):
        # One O allele, one A/B allele
        if has_b:
            abo = "B"
        else:
            abo = "A"
    elif abo_proxy == "CC":
        # No O alleles
        if b_allele == "TT":
            abo = "B"
        elif has_b:
            abo = "AB"
        else:
            abo = "A"
    else:
        # No ABO proxy data
        if has_b:
            abo = "B"
        else:
            abo = None

    # --- Rh determination ---
    rh_snp = genome_by_rsid.get("rs590787", {}).get("genotype", "")
    rh = None
    rh_confidence = 0

    if rh_snp:
        details.append(f"rs590787 (Rh proxy): {rh_snp}")
        rh_confidence = 1
        if rh_snp == "TT":
            rh = "-"
        else:
            rh = "+"

    # --- Assemble result ---
    if abo and rh:
        blood_type = f"{abo}{rh}"
    elif abo:
        blood_type = f"{abo}?"
    elif rh:
        blood_type = f"?{rh}"
    else:
        blood_type = "Unknown"

    total_confidence = abo_confidence + rh_confidence
    if total_confidence >= 3:
        confidence = "high"
    elif total_confidence >= 1:
        confidence = "moderate"
    else:
        confidence = "low"

    return {
        "blood_type": blood_type,
        "abo": abo or "Unknown",
        "rh": rh or "Unknown",
        "confidence": confidence,
        "details": details,
    }
