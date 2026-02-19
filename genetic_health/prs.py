"""Polygenic Risk Score (PRS) calculation for common conditions.

Implements simplified PRS models using published GWAS effect sizes
(log odds ratios) for 8 well-studied conditions. Each model uses
~20-30 top-associated SNPs.

Limitations:
  - No LD clumping: SNPs in the same LD block are deduplicated by rsID but
    not pruned by r². Nearby SNPs may be correlated, inflating variance.
  - EUR-calibrated: effect sizes and frequencies from European-ancestry GWAS.
    Percentiles are adjusted by a continuous ancestry factor for non-EUR.
  - Confidence intervals are approximate (assumes independent SNPs).

References per condition are embedded in PRS_MODELS.
"""

import math

# =============================================================================
# PRS MODELS
# =============================================================================
# Each condition has:
#   name: display name
#   reference: primary GWAS publication
#   snps: list of {rsid, risk_allele, log_or, gene, eur_freq}
#
# log_or values are natural log of odds ratio from published GWAS.
# eur_freq is European risk allele frequency (for population mean/var).

PRS_MODELS = {
    "type2_diabetes": {
        "name": "Type 2 Diabetes",
        "reference": "Mahajan et al. 2018 (Nature Genetics)",
        "snps": [
            {"rsid": "rs7903146", "risk_allele": "T", "log_or": 0.322, "gene": "TCF7L2", "eur_freq": 0.290},
            {"rsid": "rs13266634", "risk_allele": "C", "log_or": 0.131, "gene": "SLC30A8", "eur_freq": 0.700},
            {"rsid": "rs1111875", "risk_allele": "C", "log_or": 0.113, "gene": "HHEX", "eur_freq": 0.560},
            {"rsid": "rs10811661", "risk_allele": "T", "log_or": 0.182, "gene": "CDKN2A/B", "eur_freq": 0.830},
            {"rsid": "rs5219", "risk_allele": "T", "log_or": 0.148, "gene": "KCNJ11", "eur_freq": 0.370},
            {"rsid": "rs4402960", "risk_allele": "T", "log_or": 0.113, "gene": "IGF2BP2", "eur_freq": 0.310},
            {"rsid": "rs12255372", "risk_allele": "T", "log_or": 0.262, "gene": "TCF7L2", "eur_freq": 0.260},
            {"rsid": "rs9939609", "risk_allele": "A", "log_or": 0.131, "gene": "FTO", "eur_freq": 0.415},
            {"rsid": "rs7756992", "risk_allele": "G", "log_or": 0.157, "gene": "CDKAL1", "eur_freq": 0.310},
            {"rsid": "rs10946398", "risk_allele": "C", "log_or": 0.148, "gene": "CDKAL1", "eur_freq": 0.310},
            {"rsid": "rs4607103", "risk_allele": "C", "log_or": 0.105, "gene": "ADAMTS9", "eur_freq": 0.760},
            {"rsid": "rs2943641", "risk_allele": "C", "log_or": 0.105, "gene": "IRS1", "eur_freq": 0.630},
            {"rsid": "rs1801282", "risk_allele": "C", "log_or": 0.148, "gene": "PPARG", "eur_freq": 0.880},
            {"rsid": "rs8050136", "risk_allele": "A", "log_or": 0.131, "gene": "FTO", "eur_freq": 0.410},
            {"rsid": "rs10830963", "risk_allele": "G", "log_or": 0.095, "gene": "MTNR1B", "eur_freq": 0.280},
            {"rsid": "rs7578597", "risk_allele": "T", "log_or": 0.148, "gene": "THADA", "eur_freq": 0.900},
            {"rsid": "rs780094", "risk_allele": "C", "log_or": 0.078, "gene": "GCKR", "eur_freq": 0.390},
            {"rsid": "rs3802177", "risk_allele": "G", "log_or": 0.131, "gene": "SLC30A8", "eur_freq": 0.700},
            {"rsid": "rs231362", "risk_allele": "G", "log_or": 0.095, "gene": "KCNQ1", "eur_freq": 0.520},
            {"rsid": "rs2237895", "risk_allele": "C", "log_or": 0.105, "gene": "KCNQ1", "eur_freq": 0.390},
            {"rsid": "rs163184", "risk_allele": "G", "log_or": 0.095, "gene": "KCNQ1", "eur_freq": 0.490},
            {"rsid": "rs11708067", "risk_allele": "A", "log_or": 0.095, "gene": "ADCY5", "eur_freq": 0.760},
            {"rsid": "rs11605924", "risk_allele": "A", "log_or": 0.078, "gene": "CRY2", "eur_freq": 0.490},
            {"rsid": "rs7961581", "risk_allele": "C", "log_or": 0.095, "gene": "TSPAN8/LGR5", "eur_freq": 0.270},
            {"rsid": "rs6795735", "risk_allele": "C", "log_or": 0.078, "gene": "ADAMTS9", "eur_freq": 0.560},
        ],
    },

    "coronary_artery_disease": {
        "name": "Coronary Artery Disease",
        "reference": "Nikpay et al. 2015 (Nature Genetics); Aragam et al. 2022",
        "snps": [
            {"rsid": "rs4977574", "risk_allele": "G", "log_or": 0.223, "gene": "CDKN2A/B (9p21)", "eur_freq": 0.490},
            {"rsid": "rs1333049", "risk_allele": "C", "log_or": 0.223, "gene": "9p21.3", "eur_freq": 0.470},
            {"rsid": "rs12526453", "risk_allele": "C", "log_or": 0.095, "gene": "PHACTR1", "eur_freq": 0.660},
            {"rsid": "rs6725887", "risk_allele": "C", "log_or": 0.131, "gene": "WDR12", "eur_freq": 0.150},
            {"rsid": "rs9982601", "risk_allele": "T", "log_or": 0.148, "gene": "SLC22A3-LPAL2-LPA", "eur_freq": 0.130},
            {"rsid": "rs2259816", "risk_allele": "T", "log_or": 0.078, "gene": "HNF1A", "eur_freq": 0.340},
            {"rsid": "rs501120", "risk_allele": "T", "log_or": 0.095, "gene": "CXCL12", "eur_freq": 0.820},
            {"rsid": "rs17465637", "risk_allele": "C", "log_or": 0.095, "gene": "MIA3", "eur_freq": 0.720},
            {"rsid": "rs6922269", "risk_allele": "A", "log_or": 0.095, "gene": "MTHFD1L", "eur_freq": 0.260},
            {"rsid": "rs2505083", "risk_allele": "C", "log_or": 0.065, "gene": "KIAA1462", "eur_freq": 0.420},
            {"rsid": "rs3184504", "risk_allele": "T", "log_or": 0.078, "gene": "SH2B3", "eur_freq": 0.480},
            {"rsid": "rs11206510", "risk_allele": "T", "log_or": 0.095, "gene": "PCSK9", "eur_freq": 0.820},
            {"rsid": "rs515135", "risk_allele": "G", "log_or": 0.095, "gene": "APOB", "eur_freq": 0.790},
            {"rsid": "rs46522", "risk_allele": "T", "log_or": 0.065, "gene": "UBE2Z", "eur_freq": 0.520},
            {"rsid": "rs2954029", "risk_allele": "A", "log_or": 0.078, "gene": "TRIB1", "eur_freq": 0.540},
            {"rsid": "rs1746048", "risk_allele": "C", "log_or": 0.095, "gene": "CXCL12", "eur_freq": 0.870},
            {"rsid": "rs12190287", "risk_allele": "C", "log_or": 0.078, "gene": "TCF21", "eur_freq": 0.620},
            {"rsid": "rs7025486", "risk_allele": "A", "log_or": 0.065, "gene": "DAB2IP", "eur_freq": 0.260},
            {"rsid": "rs12413409", "risk_allele": "G", "log_or": 0.095, "gene": "CYP17A1-NT5C2", "eur_freq": 0.890},
            {"rsid": "rs964184", "risk_allele": "G", "log_or": 0.105, "gene": "ZPR1/ZNF259", "eur_freq": 0.130},
            {"rsid": "rs4773144", "risk_allele": "G", "log_or": 0.065, "gene": "COL4A1/A2", "eur_freq": 0.440},
            {"rsid": "rs1122608", "risk_allele": "G", "log_or": 0.095, "gene": "LDLR", "eur_freq": 0.770},
            {"rsid": "rs6544713", "risk_allele": "T", "log_or": 0.078, "gene": "SWAP70", "eur_freq": 0.310},
            {"rsid": "rs2246833", "risk_allele": "T", "log_or": 0.065, "gene": "LIPA", "eur_freq": 0.380},
            {"rsid": "rs264", "risk_allele": "G", "log_or": 0.078, "gene": "LPL", "eur_freq": 0.860},
        ],
    },

    "hypertension": {
        "name": "Hypertension",
        "reference": "Evangelou et al. 2018 (Nature Genetics)",
        "snps": [
            {"rsid": "rs699", "risk_allele": "G", "log_or": 0.039, "gene": "AGT", "eur_freq": 0.410},
            {"rsid": "rs5186", "risk_allele": "C", "log_or": 0.049, "gene": "AGTR1", "eur_freq": 0.285},
            {"rsid": "rs1799998", "risk_allele": "T", "log_or": 0.044, "gene": "CYP11B2", "eur_freq": 0.440},
            {"rsid": "rs4961", "risk_allele": "G", "log_or": 0.039, "gene": "ADD1", "eur_freq": 0.200},
            {"rsid": "rs1378942", "risk_allele": "C", "log_or": 0.049, "gene": "CSK-ULK3", "eur_freq": 0.650},
            {"rsid": "rs11191548", "risk_allele": "T", "log_or": 0.059, "gene": "NT5C2", "eur_freq": 0.920},
            {"rsid": "rs17249754", "risk_allele": "G", "log_or": 0.054, "gene": "ATP2B1", "eur_freq": 0.840},
            {"rsid": "rs12946454", "risk_allele": "T", "log_or": 0.034, "gene": "PLCD3", "eur_freq": 0.260},
            {"rsid": "rs16998073", "risk_allele": "T", "log_or": 0.039, "gene": "FGF5", "eur_freq": 0.290},
            {"rsid": "rs1004467", "risk_allele": "A", "log_or": 0.044, "gene": "CYP17A1", "eur_freq": 0.900},
            {"rsid": "rs3184504", "risk_allele": "T", "log_or": 0.039, "gene": "SH2B3", "eur_freq": 0.480},
            {"rsid": "rs11556924", "risk_allele": "T", "log_or": 0.034, "gene": "ZC3HC1", "eur_freq": 0.620},
            {"rsid": "rs4373814", "risk_allele": "G", "log_or": 0.034, "gene": "CACNB2", "eur_freq": 0.550},
            {"rsid": "rs12940887", "risk_allele": "T", "log_or": 0.034, "gene": "ZNF652", "eur_freq": 0.390},
            {"rsid": "rs1530440", "risk_allele": "T", "log_or": 0.039, "gene": "c10orf107", "eur_freq": 0.180},
            {"rsid": "rs13082711", "risk_allele": "T", "log_or": 0.034, "gene": "SLC4A7", "eur_freq": 0.220},
            {"rsid": "rs381815", "risk_allele": "T", "log_or": 0.029, "gene": "PLEKHA7", "eur_freq": 0.260},
            {"rsid": "rs805303", "risk_allele": "G", "log_or": 0.029, "gene": "BAT2/BAT5", "eur_freq": 0.620},
            {"rsid": "rs2681472", "risk_allele": "A", "log_or": 0.054, "gene": "ATP2B1", "eur_freq": 0.790},
            {"rsid": "rs13139571", "risk_allele": "A", "log_or": 0.034, "gene": "CACNA1C", "eur_freq": 0.370},
            {"rsid": "rs932764", "risk_allele": "G", "log_or": 0.029, "gene": "PLCE1", "eur_freq": 0.450},
            {"rsid": "rs2521501", "risk_allele": "T", "log_or": 0.039, "gene": "FES", "eur_freq": 0.300},
        ],
    },

    "breast_cancer": {
        "name": "Breast Cancer",
        "reference": "Michailidou et al. 2017 (Nature)",
        "snps": [
            {"rsid": "rs2981582", "risk_allele": "A", "log_or": 0.262, "gene": "FGFR2", "eur_freq": 0.380},
            {"rsid": "rs3803662", "risk_allele": "A", "log_or": 0.198, "gene": "TOX3/TNRC9", "eur_freq": 0.260},
            {"rsid": "rs889312", "risk_allele": "C", "log_or": 0.113, "gene": "MAP3K1", "eur_freq": 0.280},
            {"rsid": "rs13281615", "risk_allele": "G", "log_or": 0.078, "gene": "8q24", "eur_freq": 0.400},
            {"rsid": "rs3817198", "risk_allele": "C", "log_or": 0.078, "gene": "LSP1", "eur_freq": 0.310},
            {"rsid": "rs4973768", "risk_allele": "T", "log_or": 0.078, "gene": "SLC4A7/NEK10", "eur_freq": 0.470},
            {"rsid": "rs10941679", "risk_allele": "G", "log_or": 0.095, "gene": "5p12", "eur_freq": 0.260},
            {"rsid": "rs2046210", "risk_allele": "A", "log_or": 0.095, "gene": "ESR1", "eur_freq": 0.340},
            {"rsid": "rs2981578", "risk_allele": "G", "log_or": 0.198, "gene": "FGFR2", "eur_freq": 0.400},
            {"rsid": "rs1011970", "risk_allele": "T", "log_or": 0.065, "gene": "CDKN2A/B", "eur_freq": 0.170},
            {"rsid": "rs704010", "risk_allele": "A", "log_or": 0.065, "gene": "ZMIZ1", "eur_freq": 0.380},
            {"rsid": "rs614367", "risk_allele": "T", "log_or": 0.148, "gene": "11q13", "eur_freq": 0.150},
            {"rsid": "rs999737", "risk_allele": "C", "log_or": 0.095, "gene": "RAD51L1", "eur_freq": 0.760},
            {"rsid": "rs11249433", "risk_allele": "G", "log_or": 0.078, "gene": "1p11.2", "eur_freq": 0.400},
            {"rsid": "rs865686", "risk_allele": "T", "log_or": 0.095, "gene": "9q31.2", "eur_freq": 0.600},
            {"rsid": "rs10995190", "risk_allele": "G", "log_or": 0.131, "gene": "ZNF365", "eur_freq": 0.850},
            {"rsid": "rs1562430", "risk_allele": "T", "log_or": 0.078, "gene": "8q24", "eur_freq": 0.580},
            {"rsid": "rs6001930", "risk_allele": "T", "log_or": 0.095, "gene": "MKL1/MYBL1", "eur_freq": 0.230},
            {"rsid": "rs17530068", "risk_allele": "A", "log_or": 0.065, "gene": "6q14", "eur_freq": 0.250},
            {"rsid": "rs4808801", "risk_allele": "A", "log_or": 0.065, "gene": "EBF1/RNF146", "eur_freq": 0.350},
            {"rsid": "rs11820646", "risk_allele": "T", "log_or": 0.054, "gene": "2q35", "eur_freq": 0.400},
            {"rsid": "rs941764", "risk_allele": "G", "log_or": 0.054, "gene": "CCND1", "eur_freq": 0.330},
        ],
    },

    "age_related_macular_degeneration": {
        "name": "Age-Related Macular Degeneration",
        "reference": "Fritsche et al. 2016 (Nature Genetics)",
        "snps": [
            {"rsid": "rs10490924", "risk_allele": "T", "log_or": 0.693, "gene": "ARMS2/HTRA1", "eur_freq": 0.210},
            {"rsid": "rs1061170", "risk_allele": "C", "log_or": 0.470, "gene": "CFH", "eur_freq": 0.350},
            {"rsid": "rs2230199", "risk_allele": "G", "log_or": 0.262, "gene": "C3", "eur_freq": 0.200},
            {"rsid": "rs10033900", "risk_allele": "T", "log_or": 0.131, "gene": "CFI", "eur_freq": 0.500},
            {"rsid": "rs429608", "risk_allele": "G", "log_or": 0.336, "gene": "C2/CFB", "eur_freq": 0.050},
            {"rsid": "rs943080", "risk_allele": "T", "log_or": 0.148, "gene": "VEGFA", "eur_freq": 0.540},
            {"rsid": "rs3764261", "risk_allele": "A", "log_or": 0.131, "gene": "CETP", "eur_freq": 0.330},
            {"rsid": "rs13278062", "risk_allele": "T", "log_or": 0.131, "gene": "TNFRSF10A", "eur_freq": 0.490},
            {"rsid": "rs4698775", "risk_allele": "G", "log_or": 0.148, "gene": "CFI", "eur_freq": 0.470},
            {"rsid": "rs920299", "risk_allele": "C", "log_or": 0.095, "gene": "LIPC", "eur_freq": 0.470},
            {"rsid": "rs5749482", "risk_allele": "C", "log_or": 0.182, "gene": "TIMP3", "eur_freq": 0.950},
            {"rsid": "rs3812111", "risk_allele": "T", "log_or": 0.095, "gene": "COL10A1", "eur_freq": 0.430},
            {"rsid": "rs8135665", "risk_allele": "T", "log_or": 0.095, "gene": "SLC16A8", "eur_freq": 0.220},
            {"rsid": "rs1999930", "risk_allele": "G", "log_or": 0.095, "gene": "FRK/COL10A1", "eur_freq": 0.740},
            {"rsid": "rs6795735", "risk_allele": "C", "log_or": 0.078, "gene": "ADAMTS9", "eur_freq": 0.560},
            {"rsid": "rs334353", "risk_allele": "T", "log_or": 0.078, "gene": "TGFBR1", "eur_freq": 0.220},
        ],
    },
    "prostate_cancer": {
        "name": "Prostate Cancer",
        "reference": "Schumacher et al. 2018 (Nature Genetics)",
        "snps": [
            {"rsid": "rs10993994", "risk_allele": "T", "log_or": 0.198, "gene": "MSMB", "eur_freq": 0.380},
            {"rsid": "rs1447295", "risk_allele": "A", "log_or": 0.157, "gene": "8q24", "eur_freq": 0.120},
            {"rsid": "rs6983267", "risk_allele": "G", "log_or": 0.131, "gene": "8q24", "eur_freq": 0.500},
            {"rsid": "rs16901979", "risk_allele": "A", "log_or": 0.336, "gene": "8q24", "eur_freq": 0.030},
            {"rsid": "rs4430796", "risk_allele": "A", "log_or": 0.148, "gene": "HNF1B", "eur_freq": 0.490},
            {"rsid": "rs1859962", "risk_allele": "G", "log_or": 0.131, "gene": "17q24.3", "eur_freq": 0.460},
            {"rsid": "rs2735839", "risk_allele": "G", "log_or": 0.148, "gene": "KLK3 (PSA)", "eur_freq": 0.850},
            {"rsid": "rs17632542", "risk_allele": "T", "log_or": 0.182, "gene": "KLK3", "eur_freq": 0.920},
            {"rsid": "rs12621278", "risk_allele": "A", "log_or": 0.223, "gene": "ITGA6", "eur_freq": 0.060},
            {"rsid": "rs10486567", "risk_allele": "G", "log_or": 0.131, "gene": "JAZF1", "eur_freq": 0.770},
            {"rsid": "rs2660753", "risk_allele": "T", "log_or": 0.148, "gene": "3p12", "eur_freq": 0.100},
            {"rsid": "rs7679673", "risk_allele": "C", "log_or": 0.113, "gene": "TET2", "eur_freq": 0.550},
            {"rsid": "rs11649743", "risk_allele": "G", "log_or": 0.113, "gene": "HNF1B", "eur_freq": 0.740},
            {"rsid": "rs721048", "risk_allele": "A", "log_or": 0.148, "gene": "EHBP1", "eur_freq": 0.190},
            {"rsid": "rs9364554", "risk_allele": "T", "log_or": 0.095, "gene": "SLC22A3", "eur_freq": 0.280},
            {"rsid": "rs12418451", "risk_allele": "A", "log_or": 0.113, "gene": "11q13", "eur_freq": 0.270},
            {"rsid": "rs8102476", "risk_allele": "C", "log_or": 0.095, "gene": "19q13", "eur_freq": 0.530},
            {"rsid": "rs11568818", "risk_allele": "G", "log_or": 0.078, "gene": "MMP7", "eur_freq": 0.470},
            {"rsid": "rs7931342", "risk_allele": "G", "log_or": 0.095, "gene": "10q11", "eur_freq": 0.510},
            {"rsid": "rs2242652", "risk_allele": "A", "log_or": 0.113, "gene": "TERT", "eur_freq": 0.200},
        ],
    },

    "ischemic_stroke": {
        "name": "Ischemic Stroke",
        "reference": "Malik et al. 2018 (Nature Genetics)",
        "snps": [
            {"rsid": "rs12122341", "risk_allele": "T", "log_or": 0.095, "gene": "SH2B3", "eur_freq": 0.360},
            {"rsid": "rs879324", "risk_allele": "A", "log_or": 0.078, "gene": "PITX2", "eur_freq": 0.290},
            {"rsid": "rs2107595", "risk_allele": "A", "log_or": 0.095, "gene": "HDAC9", "eur_freq": 0.190},
            {"rsid": "rs12932445", "risk_allele": "A", "log_or": 0.078, "gene": "ZFHX3", "eur_freq": 0.210},
            {"rsid": "rs7859727", "risk_allele": "T", "log_or": 0.065, "gene": "TSPAN2", "eur_freq": 0.380},
            {"rsid": "rs4977574", "risk_allele": "G", "log_or": 0.095, "gene": "CDKN2A/B (9p21)", "eur_freq": 0.490},
            {"rsid": "rs11984041", "risk_allele": "T", "log_or": 0.182, "gene": "ABO", "eur_freq": 0.070},
            {"rsid": "rs6825454", "risk_allele": "C", "log_or": 0.065, "gene": "WNT2B", "eur_freq": 0.450},
            {"rsid": "rs1799983", "risk_allele": "T", "log_or": 0.078, "gene": "NOS3", "eur_freq": 0.320},
            {"rsid": "rs505922", "risk_allele": "C", "log_or": 0.095, "gene": "ABO", "eur_freq": 0.370},
            {"rsid": "rs3184504", "risk_allele": "T", "log_or": 0.065, "gene": "SH2B3", "eur_freq": 0.480},
            {"rsid": "rs1801020", "risk_allele": "A", "log_or": 0.078, "gene": "F12", "eur_freq": 0.230},
            {"rsid": "rs9349379", "risk_allele": "G", "log_or": 0.078, "gene": "PHACTR1", "eur_freq": 0.580},
            {"rsid": "rs635634", "risk_allele": "T", "log_or": 0.095, "gene": "ABO", "eur_freq": 0.200},
            {"rsid": "rs11556924", "risk_allele": "T", "log_or": 0.054, "gene": "ZC3HC1", "eur_freq": 0.620},
        ],
    },

    "colorectal_cancer": {
        "name": "Colorectal Cancer",
        "reference": "Huyghe et al. 2019 (Nature Genetics)",
        "snps": [
            {"rsid": "rs6983267", "risk_allele": "G", "log_or": 0.182, "gene": "8q24/MYC", "eur_freq": 0.500},
            {"rsid": "rs7014346", "risk_allele": "A", "log_or": 0.131, "gene": "8q24", "eur_freq": 0.360},
            {"rsid": "rs4939827", "risk_allele": "T", "log_or": 0.113, "gene": "SMAD7", "eur_freq": 0.530},
            {"rsid": "rs10411210", "risk_allele": "T", "log_or": 0.131, "gene": "RHPN2", "eur_freq": 0.900},
            {"rsid": "rs961253", "risk_allele": "A", "log_or": 0.113, "gene": "20p12.3", "eur_freq": 0.370},
            {"rsid": "rs6691170", "risk_allele": "T", "log_or": 0.095, "gene": "1q41", "eur_freq": 0.640},
            {"rsid": "rs4444235", "risk_allele": "C", "log_or": 0.095, "gene": "BMP4", "eur_freq": 0.470},
            {"rsid": "rs9929218", "risk_allele": "G", "log_or": 0.095, "gene": "CDH1", "eur_freq": 0.710},
            {"rsid": "rs4779584", "risk_allele": "T", "log_or": 0.131, "gene": "GREM1", "eur_freq": 0.190},
            {"rsid": "rs10795668", "risk_allele": "G", "log_or": 0.095, "gene": "10p14", "eur_freq": 0.670},
            {"rsid": "rs3802842", "risk_allele": "C", "log_or": 0.105, "gene": "11q23.1", "eur_freq": 0.270},
            {"rsid": "rs4813802", "risk_allele": "G", "log_or": 0.078, "gene": "20p12.3", "eur_freq": 0.340},
            {"rsid": "rs2427308", "risk_allele": "C", "log_or": 0.078, "gene": "20q13.33", "eur_freq": 0.470},
            {"rsid": "rs10936599", "risk_allele": "C", "log_or": 0.078, "gene": "MYNN/TERC", "eur_freq": 0.730},
            {"rsid": "rs16892766", "risk_allele": "C", "log_or": 0.198, "gene": "8q23.3 (EIF3H)", "eur_freq": 0.070},
            {"rsid": "rs3217810", "risk_allele": "T", "log_or": 0.148, "gene": "CCND2", "eur_freq": 0.100},
            {"rsid": "rs992157", "risk_allele": "C", "log_or": 0.078, "gene": "2q32.3", "eur_freq": 0.360},
            {"rsid": "rs1800469", "risk_allele": "A", "log_or": 0.065, "gene": "TGFB1", "eur_freq": 0.310},
            {"rsid": "rs7136702", "risk_allele": "T", "log_or": 0.065, "gene": "12q13.13", "eur_freq": 0.400},
            {"rsid": "rs11169552", "risk_allele": "C", "log_or": 0.095, "gene": "12q24.21", "eur_freq": 0.730},
        ],
    },
}


# =============================================================================
# SCORING
# =============================================================================

def _count_risk_allele(genotype: str, risk_allele: str) -> int:
    """Count copies of risk allele (0, 1, or 2)."""
    return sum(1 for a in genotype if a == risk_allele)


def _z_to_percentile(z: float) -> float:
    """Convert Z-score to percentile using error function (stdlib math)."""
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0))) * 100


def _categorize_percentile(percentile: float) -> str:
    """Map percentile to risk category."""
    if percentile < 20:
        return "low"
    elif percentile < 80:
        return "average"
    elif percentile < 95:
        return "elevated"
    else:
        return "high"


def _validate_models():
    """Validate PRS model definitions at import time."""
    for cid, model in PRS_MODELS.items():
        seen = set()
        for snp in model["snps"]:
            if snp["rsid"] in seen:
                raise ValueError(f"Duplicate rsID {snp['rsid']} in PRS model {cid}")
            seen.add(snp["rsid"])
            if not (-1 < snp["log_or"] < 1):
                raise ValueError(
                    f"Unexpected effect size {snp['log_or']} for {snp['rsid']} in {cid}"
                )

_validate_models()


def calculate_prs(genome_by_rsid: dict, ancestry_proportions: dict = None) -> dict:
    """Calculate polygenic risk scores for all conditions.

    Algorithm per condition:
        1. Deduplicate SNPs by rsID within each model
        2. raw_score = sum(log_or * risk_allele_count) for found SNPs
        3. population_mean = sum(log_or * 2 * eur_freq) for found SNPs
        4. population_var = sum(log_or^2 * 2 * eur_freq * (1-eur_freq)) for found SNPs
        5. z-score = (raw - mean) / sqrt(var)
        6. percentile via math.erf (no scipy needed)
        7. 95% confidence interval from standard error
        8. Continuous ancestry adjustment (scales z-score by EUR proportion)

    Args:
        genome_by_rsid: Dict mapping rsID -> {genotype, chromosome, position}
        ancestry_proportions: Optional dict of population proportions from ancestry module

    Returns:
        Dict mapping condition_id -> result dict
    """
    eur_prop = 1.0
    if ancestry_proportions:
        eur_prop = ancestry_proportions.get("EUR", 1.0)

    results = {}

    for condition_id, model in PRS_MODELS.items():
        raw_score = 0.0
        pop_mean = 0.0
        pop_var = 0.0
        snps_found = 0
        contributing = []
        seen_rsids = set()

        for snp in model["snps"]:
            rsid = snp["rsid"]
            if rsid in seen_rsids:
                continue
            seen_rsids.add(rsid)

            if rsid not in genome_by_rsid:
                continue

            genotype = genome_by_rsid[rsid]["genotype"]
            risk_count = _count_risk_allele(genotype, snp["risk_allele"])
            log_or = snp["log_or"]
            eur_freq = snp["eur_freq"]

            score_contribution = log_or * risk_count
            raw_score += score_contribution

            # Expected mean and variance under EUR freq
            pop_mean += log_or * 2 * eur_freq
            pop_var += (log_or ** 2) * 2 * eur_freq * (1 - eur_freq)

            snps_found += 1
            if risk_count > 0:
                contributing.append({
                    "rsid": rsid,
                    "gene": snp["gene"],
                    "risk_allele": snp["risk_allele"],
                    "copies": risk_count,
                    "log_or": log_or,
                    "contribution": score_contribution,
                })

        # Calculate percentile with confidence interval
        if snps_found > 0 and pop_var > 0:
            sd = math.sqrt(pop_var)
            z_score = (raw_score - pop_mean) / sd
            # Standard error of the PRS (approximation assuming independent SNPs)
            se = sd / math.sqrt(snps_found) if snps_found > 1 else sd
            z_se = se / sd  # SE in z-score units
            percentile = _z_to_percentile(z_score)
            ci_lower = _z_to_percentile(z_score - 1.96 * z_se)
            ci_upper = _z_to_percentile(z_score + 1.96 * z_se)
        else:
            z_score = 0.0
            percentile = 50.0
            ci_lower = 50.0
            ci_upper = 50.0

        risk_category = _categorize_percentile(percentile)

        # Continuous ancestry adjustment
        ancestry_applicable = True
        ancestry_warning = ""
        if eur_prop < 0.95:
            # Scale z-score by EUR proportion (linear decay model)
            adjusted_z = z_score * eur_prop
            adjusted_percentile = _z_to_percentile(adjusted_z)
            if eur_prop < 0.5:
                ancestry_applicable = False
                ancestry_warning = (
                    f"PRS calibrated on European populations. Your EUR ancestry "
                    f"is ~{eur_prop:.0%}, so these scores have reduced accuracy. "
                    f"Unadjusted percentile: {percentile:.0f}th."
                )
            else:
                ancestry_warning = (
                    f"PRS adjusted for {eur_prop:.0%} European ancestry "
                    f"(unadjusted: {percentile:.0f}th percentile)."
                )
            percentile = round(adjusted_percentile, 1)
            risk_category = _categorize_percentile(percentile)

        # Sort contributing SNPs by contribution (descending)
        contributing.sort(key=lambda x: -x["contribution"])

        # Unique SNP count in model (after dedup)
        unique_total = len(seen_rsids)

        results[condition_id] = {
            "name": model["name"],
            "raw_score": round(raw_score, 4),
            "z_score": round(z_score, 3),
            "percentile": round(percentile, 1),
            "ci_95_lower": round(ci_lower, 1),
            "ci_95_upper": round(ci_upper, 1),
            "risk_category": risk_category,
            "snps_found": snps_found,
            "snps_total": unique_total,
            "ancestry_applicable": ancestry_applicable,
            "ancestry_warning": ancestry_warning,
            "contributing_snps": contributing[:10],
            "reference": model["reference"],
        }

    return results
