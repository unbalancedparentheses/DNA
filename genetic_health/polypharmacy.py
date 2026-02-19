"""Polypharmacy risk assessment from combined pharmacogenomic profiles.

Flags dangerous drug combinations based on multiple pharmacogene variants
acting on the same metabolic pathway or drug class.
"""


# Each rule defines a multi-gene interaction that creates polypharmacy risk.
# genes: dict of gene -> list of qualifying phenotypes/statuses
# drugs_affected: list of drug classes at risk
# severity: "high", "moderate", "low"
# warning: user-facing explanation
# action: clinical recommendation

POLYPHARMACY_RULES = [
    {
        "id": "warfarin_compound",
        "name": "Warfarin Compound Sensitivity",
        "genes": {
            "CYP2C9": ["poor", "intermediate"],
            "VKORC1": ["sensitive", "highly_sensitive"],
        },
        "gene_source": {"CYP2C9": "star_alleles", "VKORC1": "lifestyle"},
        "drugs_affected": ["warfarin", "acenocoumarol"],
        "severity": "high",
        "warning": (
            "CYP2C9 reduced metabolism + VKORC1 increased sensitivity = "
            "compounded warfarin sensitivity. Standard doses may cause "
            "dangerous bleeding. Require 50-80% dose reduction."
        ),
        "action": (
            "Use FDA-approved warfarin dosing algorithm incorporating both "
            "CYP2C9 and VKORC1 genotypes. Consider direct oral anticoagulants "
            "(DOACs) as alternatives. INR monitoring critical."
        ),
    },
    {
        "id": "opioid_sensitivity",
        "name": "Opioid Response Complexity",
        "genes": {
            "CYP2D6": ["poor", "intermediate"],
            "OPRM1": ["reduced", "altered"],
        },
        "gene_source": {"CYP2D6": "star_alleles", "OPRM1": "lifestyle"},
        "drugs_affected": ["codeine", "tramadol", "oxycodone", "hydrocodone"],
        "severity": "high",
        "warning": (
            "Poor CYP2D6 metabolism prevents conversion of codeine/tramadol "
            "to active metabolites. Combined with altered opioid receptor "
            "sensitivity, pain relief may be minimal or unpredictable."
        ),
        "action": (
            "Avoid codeine and tramadol (prodrugs requiring CYP2D6). "
            "Use morphine or non-opioid alternatives. If opioids necessary, "
            "titrate carefully with pain specialist guidance."
        ),
    },
    {
        "id": "ssri_sensitivity",
        "name": "SSRI/Antidepressant Sensitivity Triad",
        "genes": {
            "CYP2C19": ["poor", "intermediate"],
            "COMT": ["slow"],
        },
        "gene_source": {"CYP2C19": "star_alleles", "COMT": "lifestyle"},
        "drugs_affected": ["citalopram", "escitalopram", "sertraline", "amitriptyline"],
        "severity": "moderate",
        "warning": (
            "Slow CYP2C19 metabolism increases SSRI blood levels. Combined "
            "with slow COMT (elevated catecholamines), serotonin syndrome "
            "risk is elevated. Activation side effects more likely."
        ),
        "action": (
            "Start SSRIs at 50% standard dose. Monitor for activation, "
            "agitation, or serotonergic symptoms. CYP2C19-metabolized SSRIs "
            "(citalopram, escitalopram) especially affected. Consider "
            "alternatives like bupropion or mirtazapine."
        ),
    },
    {
        "id": "statin_myopathy",
        "name": "Statin Myopathy Risk",
        "genes": {
            "SLCO1B1": ["reduced", "poor_transporter"],
        },
        "gene_source": {"SLCO1B1": "lifestyle"},
        "drugs_affected": ["simvastatin", "atorvastatin", "rosuvastatin"],
        "severity": "moderate",
        "warning": (
            "Reduced SLCO1B1 transport increases statin blood levels, "
            "raising myopathy risk 5-17x with simvastatin. Higher risk "
            "when combined with CYP3A4 inhibitors (clarithromycin, "
            "itraconazole, grapefruit)."
        ),
        "action": (
            "Avoid simvastatin >20mg. Prefer pravastatin or rosuvastatin "
            "(lower myopathy risk). Avoid concurrent CYP3A4 inhibitors. "
            "Monitor CK levels if muscle symptoms appear."
        ),
    },
    {
        "id": "chemo_toxicity",
        "name": "Chemotherapy Toxicity Risk",
        "genes": {
            "DPYD": ["poor", "intermediate"],
        },
        "gene_source": {"DPYD": "star_alleles"},
        "drugs_affected": ["5-fluorouracil", "capecitabine", "tegafur"],
        "severity": "high",
        "warning": (
            "DPYD deficiency can cause fatal fluoropyrimidine toxicity. "
            "Even partial deficiency (intermediate metabolizer) requires "
            "50% dose reduction. This is a life-threatening interaction."
        ),
        "action": (
            "MANDATORY pre-treatment DPYD testing before any fluoropyrimidine. "
            "Intermediate metabolizers: reduce dose by 50%. Poor metabolizers: "
            "contraindicated. Consider alternative chemotherapy regimens."
        ),
    },
    {
        "id": "thiopurine_toxicity",
        "name": "Thiopurine Myelosuppression Risk",
        "genes": {
            "TPMT": ["poor", "intermediate"],
        },
        "gene_source": {"TPMT": "star_alleles"},
        "drugs_affected": ["azathioprine", "6-mercaptopurine", "thioguanine"],
        "severity": "high",
        "warning": (
            "TPMT deficiency causes accumulation of cytotoxic thioguanine "
            "nucleotides, leading to life-threatening myelosuppression. "
            "Even intermediate metabolizers need dose adjustments."
        ),
        "action": (
            "Intermediate metabolizers: reduce thiopurine dose by 30-50%. "
            "Poor metabolizers: reduce by 90% or use alternative. "
            "Monitor CBC weekly for first 8 weeks of therapy."
        ),
    },
    {
        "id": "cyp2d6_ultrarapid_codeine",
        "name": "Ultrarapid CYP2D6 + Codeine Toxicity",
        "genes": {
            "CYP2D6": ["ultrarapid", "rapid"],
        },
        "gene_source": {"CYP2D6": "star_alleles"},
        "drugs_affected": ["codeine", "tramadol"],
        "severity": "high",
        "warning": (
            "Ultrarapid CYP2D6 metabolism converts codeine to morphine "
            "too quickly, causing respiratory depression. Fatal cases "
            "reported, especially in children and breastfeeding mothers."
        ),
        "action": (
            "AVOID codeine completely. Use non-CYP2D6 alternatives: "
            "morphine (already active), NSAIDs, or acetaminophen. "
            "If breastfeeding, codeine is contraindicated."
        ),
    },
    {
        "id": "ppi_reduced_efficacy",
        "name": "PPI Reduced Efficacy (Ultrarapid CYP2C19)",
        "genes": {
            "CYP2C19": ["ultrarapid", "rapid"],
        },
        "gene_source": {"CYP2C19": "star_alleles"},
        "drugs_affected": ["omeprazole", "lansoprazole", "pantoprazole"],
        "severity": "moderate",
        "warning": (
            "Ultrarapid CYP2C19 clears PPIs too fast, reducing acid "
            "suppression. H. pylori eradication rates drop significantly. "
            "Standard PPI doses may be insufficient."
        ),
        "action": (
            "Increase PPI dose or switch to rabeprazole (less CYP2C19 "
            "dependent). For H. pylori: use higher-dose PPI regimens. "
            "Consider vonoprazan if available."
        ),
    },
    {
        "id": "clopidogrel_resistance",
        "name": "Clopidogrel Resistance (Poor CYP2C19)",
        "genes": {
            "CYP2C19": ["poor", "intermediate"],
        },
        "gene_source": {"CYP2C19": "star_alleles"},
        "drugs_affected": ["clopidogrel"],
        "severity": "high",
        "warning": (
            "CYP2C19 is required to activate clopidogrel. Poor/intermediate "
            "metabolizers have reduced antiplatelet effect, increasing risk "
            "of stent thrombosis and cardiovascular events."
        ),
        "action": (
            "Use alternative antiplatelet: prasugrel or ticagrelor "
            "(not CYP2C19 dependent). FDA black box warning on clopidogrel "
            "for poor CYP2C19 metabolizers."
        ),
    },
]


def assess_polypharmacy(genome_by_rsid, star_alleles=None, lifestyle_findings=None):
    """Assess polypharmacy risks from combined pharmacogenomic profiles.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict.
    star_alleles : dict or None
        Star allele results from call_star_alleles().
    lifestyle_findings : list or None
        Lifestyle/health findings from analyze_lifestyle_health().

    Returns
    -------
    dict with keys: warnings, total_warnings, by_severity, drug_card_additions.
    """
    if lifestyle_findings is None:
        lifestyle_findings = []
    if star_alleles is None:
        star_alleles = {}

    # Build gene -> status/phenotype lookup
    gene_phenotypes = {}

    # From star alleles
    for gene, result in star_alleles.items():
        phenotype = result.get("phenotype", "normal")
        if phenotype != "Unknown":
            gene_phenotypes[gene] = phenotype

    # From lifestyle findings
    for f in lifestyle_findings:
        gene = f.get("gene", "")
        status = f.get("status", "")
        if gene and status:
            gene_phenotypes[gene] = status

    warnings = []

    for rule in POLYPHARMACY_RULES:
        match = True
        matched_genes = {}

        for gene, qualifying in rule["genes"].items():
            phenotype = gene_phenotypes.get(gene, "")
            if phenotype and phenotype in qualifying:
                matched_genes[gene] = phenotype
            else:
                match = False
                break

        if match:
            warnings.append({
                "id": rule["id"],
                "name": rule["name"],
                "severity": rule["severity"],
                "matched_genes": matched_genes,
                "drugs_affected": rule["drugs_affected"],
                "warning": rule["warning"],
                "action": rule["action"],
            })

    # Sort by severity
    severity_order = {"high": 0, "moderate": 1, "low": 2}
    warnings.sort(key=lambda w: severity_order.get(w["severity"], 3))

    by_severity = {}
    for w in warnings:
        by_severity.setdefault(w["severity"], []).append(w)

    return {
        "warnings": warnings,
        "total_warnings": len(warnings),
        "by_severity": by_severity,
    }
