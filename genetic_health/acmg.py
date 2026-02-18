"""ACMG Secondary Findings v3.2 — flag medically actionable genes.

Filters existing ClinVar pathogenic/likely_pathogenic findings against the
ACMG SF v3.2 list of 81 genes recommended for reporting as secondary findings.
"""

# ACMG SF v3.2 gene list (81 genes)
ACMG_GENES = frozenset({
    # Hereditary cancer
    "BRCA1", "BRCA2", "TP53", "MLH1", "MSH2", "MSH6", "PMS2", "APC",
    "MUTYH", "RB1", "MEN1", "RET", "VHL", "SDHB", "SDHD", "PTEN",
    "STK11", "BMPR1A", "SMAD4", "CDH1", "PALB2", "ATM", "CHEK2",
    "BAP1", "DICER1", "HOXB13", "FLCN", "MET", "MITF", "NTHL1",
    "MAX", "SDHA", "SDHAF2", "SDHC", "TMEM127", "NF2",
    # Cardiovascular
    "LDLR", "APOB", "PCSK9", "MYBPC3", "MYH7", "SCN5A", "KCNQ1",
    "KCNH2", "LMNA", "RYR1", "RYR2", "CACNA1S", "PKP2", "DSP",
    "DSG2", "DSC2", "TMEM43", "TTN", "ACTA2", "FBN1", "TGFBR1",
    "TGFBR2", "MYH11", "COL3A1", "SMAD3", "ACVRL1", "ENG",
    "TRDN", "CASQ2",
    # Metabolic / other
    "GLA", "BTD", "OTC", "ATP7B", "HFE", "SERPINA1", "HBB",
    "GAA", "HEXA", "HEXB", "SMPD1", "RPE65", "HNF1A",
    # Tumor predisposition
    "TSC1", "TSC2", "WT1",
})

_ACTIONABILITY = {
    "BRCA1": "Hereditary breast/ovarian cancer — increased surveillance, risk-reducing surgery options",
    "BRCA2": "Hereditary breast/ovarian cancer — increased surveillance, risk-reducing surgery options",
    "TP53": "Li-Fraumeni syndrome — whole-body MRI, multi-cancer surveillance",
    "MLH1": "Lynch syndrome — colonoscopy every 1-2 years from age 25",
    "MSH2": "Lynch syndrome — colonoscopy every 1-2 years from age 25",
    "MSH6": "Lynch syndrome — colonoscopy every 1-2 years from age 25",
    "PMS2": "Lynch syndrome — colonoscopy every 1-2 years from age 25",
    "LDLR": "Familial hypercholesterolemia — early statin therapy, lipid monitoring",
    "APOB": "Familial hypercholesterolemia — early statin therapy, lipid monitoring",
    "PCSK9": "Familial hypercholesterolemia — early statin therapy, lipid monitoring",
    "MYBPC3": "Hypertrophic cardiomyopathy — echocardiography, activity modification",
    "MYH7": "Hypertrophic cardiomyopathy — echocardiography, activity modification",
    "SCN5A": "Brugada/Long QT — ECG monitoring, avoid triggering drugs",
    "KCNQ1": "Long QT syndrome — ECG monitoring, avoid QT-prolonging drugs",
    "KCNH2": "Long QT syndrome — ECG monitoring, avoid QT-prolonging drugs",
    "RYR1": "Malignant hyperthermia susceptibility — avoid triggering anesthetics",
    "FBN1": "Marfan syndrome — aortic imaging, activity modification",
    "HFE": "Hereditary hemochromatosis — iron/ferritin monitoring, phlebotomy",
    "SERPINA1": "Alpha-1 antitrypsin deficiency — avoid smoking, pulmonary function monitoring",
    "HBB": "Sickle cell/beta-thalassemia — hematology follow-up",
    "GLA": "Fabry disease — enzyme replacement therapy available",
    "ATP7B": "Wilson disease — copper studies, chelation therapy",
    "APC": "Familial adenomatous polyposis — colonoscopy from early teens",
    "RET": "MEN2 — thyroid cancer screening, prophylactic thyroidectomy consideration",
    "VHL": "Von Hippel-Lindau — multi-organ surveillance",
    "RB1": "Retinoblastoma — ophthalmologic screening",
    "PTEN": "PTEN hamartoma tumor syndrome — multi-cancer surveillance",
    "PKP2": "Arrhythmogenic cardiomyopathy — cardiac imaging, activity restriction",
    "LMNA": "Dilated cardiomyopathy/muscular dystrophy — cardiac monitoring",
    "BTD": "Biotinidase deficiency — biotin supplementation",
    "OTC": "Ornithine transcarbamylase deficiency — dietary management, emergency protocol",
    "GAA": "Pompe disease — enzyme replacement therapy",
    "RPE65": "Retinal dystrophy — gene therapy (Luxturna) available",
}


def flag_acmg_findings(disease_findings):
    """Filter ClinVar findings for ACMG SF v3.2 medically actionable genes.

    Parameters
    ----------
    disease_findings : dict
        ClinVar findings dict with keys: pathogenic, likely_pathogenic,
        risk_factor, drug_response, protective.

    Returns
    -------
    dict with keys: acmg_findings, genes_screened, genes_with_variants, summary.
    """
    if not disease_findings:
        return {
            "acmg_findings": [],
            "genes_screened": len(ACMG_GENES),
            "genes_with_variants": 0,
            "summary": "No ClinVar data available for ACMG screening.",
        }

    acmg_findings = []
    genes_found = set()

    for category in ("pathogenic", "likely_pathogenic"):
        for finding in disease_findings.get(category, []):
            gene = (finding.get("gene") or "").upper()
            if gene in ACMG_GENES:
                genes_found.add(gene)
                actionability = _ACTIONABILITY.get(gene, "Medically actionable — consult genetic counselor")
                acmg_findings.append({
                    **finding,
                    "acmg_category": category,
                    "acmg_actionability": actionability,
                })

    acmg_findings.sort(key=lambda x: (-x.get("gold_stars", 0), x.get("gene", "")))

    if acmg_findings:
        summary = (f"{len(acmg_findings)} variant(s) found in {len(genes_found)} "
                   f"ACMG-recommended gene(s). Genetic counseling recommended.")
    else:
        summary = "No pathogenic/likely pathogenic variants in ACMG SF v3.2 genes."

    return {
        "acmg_findings": acmg_findings,
        "genes_screened": len(ACMG_GENES),
        "genes_with_variants": len(genes_found),
        "summary": summary,
    }
