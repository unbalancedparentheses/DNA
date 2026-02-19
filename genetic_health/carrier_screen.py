"""Carrier screening organizer.

Reorganizes heterozygous recessive ClinVar findings into a structured
carrier screening report grouped by disease system.
"""

from .analysis import classify_zygosity

# Disease system classification by gene
_GENE_SYSTEMS = {
    # Hematologic
    "HBB": "Hematologic", "HBA1": "Hematologic", "HBA2": "Hematologic",
    "G6PD": "Hematologic", "F5": "Hematologic", "F2": "Hematologic",
    # Metabolic
    "CFTR": "Metabolic/Pulmonary", "PAH": "Metabolic", "GALT": "Metabolic",
    "GBA": "Metabolic", "HEXA": "Metabolic", "HEXB": "Metabolic",
    "SMPD1": "Metabolic", "GAA": "Metabolic", "IDUA": "Metabolic",
    "BTD": "Metabolic", "ASL": "Metabolic", "ASPA": "Metabolic",
    "ACADM": "Metabolic", "MCCC2": "Metabolic",
    # Neurologic
    "SMN1": "Neurologic", "DYNC1H1": "Neurologic",
    "FMR1": "Neurologic", "MECP2": "Neurologic",
    # Cancer predisposition
    "BRCA1": "Cancer Predisposition", "BRCA2": "Cancer Predisposition",
    "ATM": "Cancer Predisposition", "PALB2": "Cancer Predisposition",
    "CHEK2": "Cancer Predisposition", "BAP1": "Cancer Predisposition",
    "MUTYH": "Cancer Predisposition", "APC": "Cancer Predisposition",
    "MLH1": "Cancer Predisposition", "MSH2": "Cancer Predisposition",
    "MSH6": "Cancer Predisposition", "PMS2": "Cancer Predisposition",
    "TP53": "Cancer Predisposition", "PTEN": "Cancer Predisposition",
    "STK11": "Cancer Predisposition", "CDH1": "Cancer Predisposition",
    "BMPR1A": "Cancer Predisposition", "SMAD4": "Cancer Predisposition",
    # Connective tissue
    "FBN1": "Connective Tissue", "COL3A1": "Connective Tissue",
    "COL1A1": "Connective Tissue", "COL1A2": "Connective Tissue",
    "ELN": "Connective Tissue",
    # Cardiac
    "RYR1": "Cardiac/Muscular", "RYR2": "Cardiac",
    "CACNA1S": "Cardiac/Muscular", "SCN5A": "Cardiac",
    "KCNQ1": "Cardiac", "KCNH2": "Cardiac",
    "MYBPC3": "Cardiac", "MYH7": "Cardiac",
    "TNNT2": "Cardiac", "LMNA": "Cardiac",
    "DSP": "Cardiac", "PKP2": "Cardiac",
    # Renal
    "PKD1": "Renal", "PKD2": "Renal", "ATP7B": "Renal/Hepatic",
    "SLC12A3": "Renal",
    # Immune
    "SERPINA1": "Immune/Pulmonary", "AIRE": "Immune", "FOXP3": "Immune",
    # Sensory
    "GJB2": "Sensory (Hearing)", "GJB6": "Sensory (Hearing)",
    "SLC26A4": "Sensory (Hearing)", "RPE65": "Sensory (Vision)",
    "ABCA4": "Sensory (Vision)",
    # Endocrine
    "OTC": "Endocrine/Metabolic", "CYP21A2": "Endocrine",
    "HFE": "Iron Metabolism",
    # Vascular
    "ACVRL1": "Vascular", "ENG": "Vascular",
    "TGFBR1": "Vascular", "TGFBR2": "Vascular",
    "FBN1": "Connective Tissue/Vascular",
}

# Inheritance patterns
_GENE_INHERITANCE = {
    "CFTR": "autosomal recessive",
    "HBB": "autosomal recessive",
    "HEXA": "autosomal recessive",
    "HEXB": "autosomal recessive",
    "GBA": "autosomal recessive",
    "PAH": "autosomal recessive",
    "SMN1": "autosomal recessive",
    "GJB2": "autosomal recessive",
    "SERPINA1": "autosomal recessive (codominant)",
    "HFE": "autosomal recessive",
    "GAA": "autosomal recessive",
    "SMPD1": "autosomal recessive",
    "GALT": "autosomal recessive",
    "IDUA": "autosomal recessive",
    "BTD": "autosomal recessive",
    "ASL": "autosomal recessive",
    "ASPA": "autosomal recessive",
    "ACADM": "autosomal recessive",
    "ABCA4": "autosomal recessive",
    "SLC26A4": "autosomal recessive",
    "CYP21A2": "autosomal recessive",
    "ATP7B": "autosomal recessive",
    "G6PD": "X-linked",
    "OTC": "X-linked",
    "FMR1": "X-linked",
    "F5": "autosomal dominant (incomplete penetrance)",
    "F2": "autosomal dominant (incomplete penetrance)",
    "BRCA1": "autosomal dominant",
    "BRCA2": "autosomal dominant",
    "MLH1": "autosomal dominant",
    "MSH2": "autosomal dominant",
    "MSH6": "autosomal dominant",
    "PMS2": "autosomal dominant",
    "APC": "autosomal dominant",
    "MUTYH": "autosomal recessive",
    "ATM": "autosomal dominant (moderate penetrance)",
    "PALB2": "autosomal dominant (moderate penetrance)",
    "CHEK2": "autosomal dominant (moderate penetrance)",
    "TP53": "autosomal dominant",
    "RYR1": "autosomal dominant / autosomal recessive",
    "CACNA1S": "autosomal dominant",
}

# Conditions commonly screened in couples (recessive + cancer predisposition)
_COUPLES_CONDITIONS = frozenset({
    "CFTR", "HBB", "HEXA", "SMN1", "GBA", "PAH", "GALT",
    "GJB2", "SMPD1", "ASPA", "ACADM", "BTD", "GAA",
    "BRCA1", "BRCA2", "ATM", "PALB2", "CHEK2", "MUTYH",
    "MLH1", "MSH2", "MSH6", "PMS2",
})


def organize_carrier_findings(disease_findings):
    """Organize carrier findings from ClinVar into a structured report.

    Parameters
    ----------
    disease_findings : dict
        ClinVar findings dict with pathogenic, likely_pathogenic, etc.

    Returns
    -------
    dict with keys: carriers, total_carriers, by_system, couples_relevant.
    """
    if not disease_findings:
        return {
            "carriers": [],
            "total_carriers": 0,
            "by_system": {},
            "couples_relevant": [],
        }

    carriers = []

    for category in ("pathogenic", "likely_pathogenic"):
        for finding in disease_findings.get(category, []):
            status, desc = classify_zygosity(finding)
            if status != "CARRIER":
                continue

            gene = finding.get("gene", "Unknown")
            gene_upper = gene.upper()
            system = _GENE_SYSTEMS.get(gene_upper, "Other")
            inheritance = _GENE_INHERITANCE.get(
                gene_upper,
                (finding.get("inheritance") or "unknown").lower()
            )

            condition = finding.get("traits", "").split(";")[0].strip() or "Unknown condition"
            is_couples = gene_upper in _COUPLES_CONDITIONS

            reproductive_note = ""
            if "recessive" in inheritance:
                reproductive_note = ("If partner is also a carrier, each child has "
                                     "a 25% chance of being affected.")
            elif "x-linked" in inheritance:
                reproductive_note = ("X-linked: carrier females may pass to sons "
                                     "(50% chance affected) and daughters (50% chance carrier).")

            carriers.append({
                "gene": gene,
                "condition": condition,
                "inheritance": inheritance,
                "system": system,
                "reproductive_note": reproductive_note,
                "couples_relevant": is_couples,
                "rsid": finding.get("rsid", ""),
                "genotype": finding.get("user_genotype", ""),
                "gold_stars": finding.get("gold_stars", 0),
            })

    # Group by system
    by_system = {}
    for c in carriers:
        sys_name = c["system"]
        by_system.setdefault(sys_name, []).append(c)

    couples_relevant = [c for c in carriers if c["couples_relevant"]]

    return {
        "carriers": carriers,
        "total_carriers": len(carriers),
        "by_system": by_system,
        "couples_relevant": couples_relevant,
    }
