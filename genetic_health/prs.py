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

    "alzheimers_disease": {
        "name": "Alzheimer's Disease",
        "reference": "Jansen et al. 2019 (Nature Genetics); Kunkle et al. 2019",
        "snps": [
            {"rsid": "rs429358", "risk_allele": "C", "log_or": 0.693, "gene": "APOE (e4)", "eur_freq": 0.153},
            {"rsid": "rs11136000", "risk_allele": "T", "log_or": 0.148, "gene": "CLU", "eur_freq": 0.620},
            {"rsid": "rs744373", "risk_allele": "G", "log_or": 0.148, "gene": "BIN1", "eur_freq": 0.290},
            {"rsid": "rs3851179", "risk_allele": "C", "log_or": 0.131, "gene": "PICALM", "eur_freq": 0.360},
            {"rsid": "rs3764650", "risk_allele": "G", "log_or": 0.182, "gene": "ABCA7", "eur_freq": 0.100},
            {"rsid": "rs11218343", "risk_allele": "C", "log_or": 0.148, "gene": "SORL1", "eur_freq": 0.960},
            {"rsid": "rs3865444", "risk_allele": "C", "log_or": 0.095, "gene": "CD33", "eur_freq": 0.680},
            {"rsid": "rs610932", "risk_allele": "G", "log_or": 0.095, "gene": "MS4A6A", "eur_freq": 0.580},
            {"rsid": "rs6656401", "risk_allele": "A", "log_or": 0.131, "gene": "CR1", "eur_freq": 0.200},
            {"rsid": "rs1532278", "risk_allele": "T", "log_or": 0.113, "gene": "CLU", "eur_freq": 0.620},
            {"rsid": "rs2718058", "risk_allele": "G", "log_or": 0.078, "gene": "NME8", "eur_freq": 0.610},
            {"rsid": "rs10948363", "risk_allele": "G", "log_or": 0.095, "gene": "CD2AP", "eur_freq": 0.270},
            {"rsid": "rs9349407", "risk_allele": "C", "log_or": 0.095, "gene": "CD2AP", "eur_freq": 0.240},
            {"rsid": "rs28834970", "risk_allele": "C", "log_or": 0.095, "gene": "PTK2B", "eur_freq": 0.370},
        ],
    },

    "atrial_fibrillation": {
        "name": "Atrial Fibrillation",
        "reference": "Roselli et al. 2018 (Nature Genetics)",
        "snps": [
            {"rsid": "rs2200733", "risk_allele": "T", "log_or": 0.336, "gene": "PITX2 (4q25)", "eur_freq": 0.120},
            {"rsid": "rs6843082", "risk_allele": "G", "log_or": 0.157, "gene": "PITX2 (4q25)", "eur_freq": 0.160},
            {"rsid": "rs10033464", "risk_allele": "T", "log_or": 0.198, "gene": "PITX2 (4q25)", "eur_freq": 0.080},
            {"rsid": "rs3807989", "risk_allele": "A", "log_or": 0.113, "gene": "CAV1", "eur_freq": 0.390},
            {"rsid": "rs7164883", "risk_allele": "G", "log_or": 0.095, "gene": "KCNN3", "eur_freq": 0.300},
            {"rsid": "rs2106261", "risk_allele": "T", "log_or": 0.131, "gene": "ZFHX3 (16q22)", "eur_freq": 0.190},
            {"rsid": "rs6490029", "risk_allele": "G", "log_or": 0.078, "gene": "SYNPO2L", "eur_freq": 0.440},
            {"rsid": "rs2040862", "risk_allele": "C", "log_or": 0.095, "gene": "NEURL1", "eur_freq": 0.200},
            {"rsid": "rs10821415", "risk_allele": "A", "log_or": 0.078, "gene": "KCNJ5", "eur_freq": 0.320},
            {"rsid": "rs1152591", "risk_allele": "A", "log_or": 0.065, "gene": "SCN10A", "eur_freq": 0.480},
            {"rsid": "rs7508", "risk_allele": "G", "log_or": 0.065, "gene": "CUX2", "eur_freq": 0.530},
            {"rsid": "rs6584555", "risk_allele": "T", "log_or": 0.078, "gene": "SYNE2", "eur_freq": 0.270},
        ],
    },

    "chronic_kidney_disease": {
        "name": "Chronic Kidney Disease",
        "reference": "Wuttke et al. 2019 (Nature Genetics)",
        "snps": [
            {"rsid": "rs13146355", "risk_allele": "C", "log_or": 0.095, "gene": "UMOD", "eur_freq": 0.190},
            {"rsid": "rs1260326", "risk_allele": "T", "log_or": 0.054, "gene": "GCKR", "eur_freq": 0.410},
            {"rsid": "rs10206899", "risk_allele": "G", "log_or": 0.078, "gene": "NAT8/ALMS1", "eur_freq": 0.240},
            {"rsid": "rs7805747", "risk_allele": "A", "log_or": 0.065, "gene": "PRKAG2", "eur_freq": 0.260},
            {"rsid": "rs2453533", "risk_allele": "C", "log_or": 0.054, "gene": "SHROOM3", "eur_freq": 0.440},
            {"rsid": "rs6795735", "risk_allele": "C", "log_or": 0.054, "gene": "ADAMTS9", "eur_freq": 0.560},
            {"rsid": "rs4014195", "risk_allele": "G", "log_or": 0.049, "gene": "CPS1", "eur_freq": 0.320},
            {"rsid": "rs267734", "risk_allele": "A", "log_or": 0.065, "gene": "TFDP2", "eur_freq": 0.200},
            {"rsid": "rs11959928", "risk_allele": "A", "log_or": 0.044, "gene": "DAB2", "eur_freq": 0.380},
            {"rsid": "rs9895661", "risk_allele": "T", "log_or": 0.049, "gene": "MPPED2", "eur_freq": 0.270},
            {"rsid": "rs4744712", "risk_allele": "A", "log_or": 0.044, "gene": "HNF4A", "eur_freq": 0.450},
            {"rsid": "rs1801239", "risk_allele": "G", "log_or": 0.078, "gene": "CUBN", "eur_freq": 0.100},
        ],
    },

    "obesity": {
        "name": "Obesity (BMI)",
        "reference": "Yengo et al. 2018 (Human Molecular Genetics); Locke et al. 2015",
        "snps": [
            {"rsid": "rs9939609", "risk_allele": "A", "log_or": 0.131, "gene": "FTO", "eur_freq": 0.415},
            {"rsid": "rs571312", "risk_allele": "A", "log_or": 0.113, "gene": "MC4R", "eur_freq": 0.240},
            {"rsid": "rs10938397", "risk_allele": "G", "log_or": 0.078, "gene": "GNPDA2", "eur_freq": 0.440},
            {"rsid": "rs2867125", "risk_allele": "C", "log_or": 0.078, "gene": "TMEM18", "eur_freq": 0.830},
            {"rsid": "rs543874", "risk_allele": "G", "log_or": 0.078, "gene": "SEC16B", "eur_freq": 0.190},
            {"rsid": "rs7138803", "risk_allele": "A", "log_or": 0.065, "gene": "BCDIN3D/FAIM2", "eur_freq": 0.380},
            {"rsid": "rs10767664", "risk_allele": "A", "log_or": 0.078, "gene": "BDNF", "eur_freq": 0.780},
            {"rsid": "rs2815752", "risk_allele": "A", "log_or": 0.065, "gene": "NEGR1", "eur_freq": 0.620},
            {"rsid": "rs10150332", "risk_allele": "C", "log_or": 0.054, "gene": "NRXN3", "eur_freq": 0.240},
            {"rsid": "rs6265", "risk_allele": "T", "log_or": 0.065, "gene": "BDNF Val66Met", "eur_freq": 0.190},
            {"rsid": "rs1558902", "risk_allele": "A", "log_or": 0.131, "gene": "FTO (lead)", "eur_freq": 0.420},
            {"rsid": "rs13021737", "risk_allele": "G", "log_or": 0.078, "gene": "TMEM18", "eur_freq": 0.830},
            {"rsid": "rs10182181", "risk_allele": "G", "log_or": 0.065, "gene": "ADCY3-POMC", "eur_freq": 0.460},
            {"rsid": "rs11030104", "risk_allele": "A", "log_or": 0.065, "gene": "BDNF", "eur_freq": 0.790},
            {"rsid": "rs3101336", "risk_allele": "C", "log_or": 0.065, "gene": "NEGR1", "eur_freq": 0.630},
        ],
    },

    "asthma": {
        "name": "Asthma",
        "reference": "Demenais et al. 2018 (Nature Genetics); GABRIEL Consortium",
        "snps": [
            {"rsid": "rs2305480", "risk_allele": "G", "log_or": 0.148, "gene": "GSDMB/ORMDL3 (17q21)", "eur_freq": 0.470},
            {"rsid": "rs8076131", "risk_allele": "T", "log_or": 0.131, "gene": "ORMDL3/GSDMB", "eur_freq": 0.450},
            {"rsid": "rs1295686", "risk_allele": "A", "log_or": 0.095, "gene": "IL13", "eur_freq": 0.200},
            {"rsid": "rs3894194", "risk_allele": "G", "log_or": 0.078, "gene": "TSLP", "eur_freq": 0.410},
            {"rsid": "rs2284033", "risk_allele": "G", "log_or": 0.078, "gene": "IL2RB", "eur_freq": 0.450},
            {"rsid": "rs1342326", "risk_allele": "A", "log_or": 0.095, "gene": "IL33", "eur_freq": 0.110},
            {"rsid": "rs20541", "risk_allele": "A", "log_or": 0.078, "gene": "IL13", "eur_freq": 0.220},
            {"rsid": "rs7216389", "risk_allele": "T", "log_or": 0.148, "gene": "ORMDL3", "eur_freq": 0.490},
            {"rsid": "rs1101999", "risk_allele": "T", "log_or": 0.065, "gene": "MHC/HLA", "eur_freq": 0.510},
            {"rsid": "rs3117098", "risk_allele": "C", "log_or": 0.065, "gene": "HLA-DQ", "eur_freq": 0.380},
            {"rsid": "rs10197862", "risk_allele": "C", "log_or": 0.078, "gene": "IL1RL1/IL18R1", "eur_freq": 0.150},
        ],
    },

    "migraine": {
        "name": "Migraine",
        "reference": "Gormley et al. 2016 (Nature Genetics); Anttila et al. 2013",
        "snps": [
            {"rsid": "rs11172113", "risk_allele": "T", "log_or": 0.095, "gene": "LRP1", "eur_freq": 0.580},
            {"rsid": "rs10915437", "risk_allele": "A", "log_or": 0.078, "gene": "AJAP1", "eur_freq": 0.220},
            {"rsid": "rs10166942", "risk_allele": "T", "log_or": 0.113, "gene": "TRPM8", "eur_freq": 0.800},
            {"rsid": "rs2651899", "risk_allele": "C", "log_or": 0.065, "gene": "PRDM16", "eur_freq": 0.590},
            {"rsid": "rs12134493", "risk_allele": "A", "log_or": 0.095, "gene": "TSPAN2", "eur_freq": 0.150},
            {"rsid": "rs1835740", "risk_allele": "A", "log_or": 0.095, "gene": "MTDH/AEG-1", "eur_freq": 0.240},
            {"rsid": "rs6790925", "risk_allele": "T", "log_or": 0.065, "gene": "TGFBR2", "eur_freq": 0.380},
            {"rsid": "rs9349379", "risk_allele": "G", "log_or": 0.078, "gene": "PHACTR1", "eur_freq": 0.580},
            {"rsid": "rs10504861", "risk_allele": "T", "log_or": 0.065, "gene": "near REST", "eur_freq": 0.260},
            {"rsid": "rs13208321", "risk_allele": "A", "log_or": 0.065, "gene": "FHL5", "eur_freq": 0.240},
        ],
    },

    "osteoporosis": {
        "name": "Osteoporosis (Bone Mineral Density)",
        "reference": "Morris et al. 2019 (Nature Genetics); Estrada et al. 2012",
        "snps": [
            {"rsid": "rs3736228", "risk_allele": "T", "log_or": 0.131, "gene": "LRP5", "eur_freq": 0.080},
            {"rsid": "rs2282679", "risk_allele": "C", "log_or": 0.095, "gene": "GC (VitD binding)", "eur_freq": 0.280},
            {"rsid": "rs4988235", "risk_allele": "G", "log_or": 0.054, "gene": "LCT (calcium)", "eur_freq": 0.258},
            {"rsid": "rs9921222", "risk_allele": "T", "log_or": 0.078, "gene": "MEPE", "eur_freq": 0.410},
            {"rsid": "rs4355801", "risk_allele": "A", "log_or": 0.065, "gene": "TNFRSF11B (OPG)", "eur_freq": 0.520},
            {"rsid": "rs2062377", "risk_allele": "T", "log_or": 0.054, "gene": "RANK/TNFRSF11A", "eur_freq": 0.310},
            {"rsid": "rs7521902", "risk_allele": "A", "log_or": 0.065, "gene": "WNT4", "eur_freq": 0.230},
            {"rsid": "rs227584", "risk_allele": "C", "log_or": 0.049, "gene": "CTNNB1 (Wnt)", "eur_freq": 0.340},
            {"rsid": "rs2189480", "risk_allele": "G", "log_or": 0.054, "gene": "FUBP3", "eur_freq": 0.420},
            {"rsid": "rs2504063", "risk_allele": "G", "log_or": 0.049, "gene": "ESR1 (estrogen)", "eur_freq": 0.620},
        ],
    },

    "parkinsons_disease": {
        "name": "Parkinson's Disease",
        "reference": "Nalls et al. 2019 (Lancet Neurology); Chang et al. 2017",
        "snps": [
            {"rsid": "rs356182", "risk_allele": "A", "log_or": 0.287, "gene": "SNCA", "eur_freq": 0.370},
            {"rsid": "rs34637584", "risk_allele": "T", "log_or": 0.405, "gene": "LRRK2 G2019S", "eur_freq": 0.001},
            {"rsid": "rs76763715", "risk_allele": "T", "log_or": 0.693, "gene": "GBA N370S", "eur_freq": 0.003},
            {"rsid": "rs11931074", "risk_allele": "T", "log_or": 0.262, "gene": "SNCA", "eur_freq": 0.070},
            {"rsid": "rs6532194", "risk_allele": "T", "log_or": 0.157, "gene": "BST1", "eur_freq": 0.420},
            {"rsid": "rs12456492", "risk_allele": "G", "log_or": 0.113, "gene": "RIT2", "eur_freq": 0.660},
            {"rsid": "rs2414739", "risk_allele": "A", "log_or": 0.131, "gene": "PM20D1", "eur_freq": 0.670},
            {"rsid": "rs11724635", "risk_allele": "A", "log_or": 0.113, "gene": "FAM47E/SCARB2", "eur_freq": 0.540},
            {"rsid": "rs12637471", "risk_allele": "A", "log_or": 0.148, "gene": "MCCC1", "eur_freq": 0.170},
            {"rsid": "rs823118", "risk_allele": "T", "log_or": 0.113, "gene": "NUCKS1/RAB7L1", "eur_freq": 0.650},
            {"rsid": "rs199347", "risk_allele": "A", "log_or": 0.105, "gene": "GPNMB", "eur_freq": 0.430},
            {"rsid": "rs1474055", "risk_allele": "T", "log_or": 0.148, "gene": "STK39", "eur_freq": 0.120},
        ],
    },

    "type1_diabetes": {
        "name": "Type 1 Diabetes",
        "reference": "Onengut-Gumuscu et al. 2015 (Nature Genetics); Barrett et al. 2009",
        "snps": [
            {"rsid": "rs2187668", "risk_allele": "T", "log_or": 0.693, "gene": "HLA-DQA1", "eur_freq": 0.130},
            {"rsid": "rs9273363", "risk_allele": "C", "log_or": 0.588, "gene": "HLA-DQB1", "eur_freq": 0.180},
            {"rsid": "rs2476601", "risk_allele": "A", "log_or": 0.336, "gene": "PTPN22", "eur_freq": 0.100},
            {"rsid": "rs2292239", "risk_allele": "T", "log_or": 0.182, "gene": "ERBB3", "eur_freq": 0.340},
            {"rsid": "rs3184504", "risk_allele": "T", "log_or": 0.148, "gene": "SH2B3", "eur_freq": 0.480},
            {"rsid": "rs6679677", "risk_allele": "A", "log_or": 0.336, "gene": "PTPN22 (proxy)", "eur_freq": 0.100},
            {"rsid": "rs11171710", "risk_allele": "A", "log_or": 0.131, "gene": "BACH2", "eur_freq": 0.350},
            {"rsid": "rs1990760", "risk_allele": "T", "log_or": 0.148, "gene": "IFIH1", "eur_freq": 0.610},
            {"rsid": "rs12708716", "risk_allele": "G", "log_or": 0.182, "gene": "CLEC16A", "eur_freq": 0.340},
            {"rsid": "rs3087243", "risk_allele": "G", "log_or": 0.182, "gene": "CTLA4", "eur_freq": 0.540},
            {"rsid": "rs6691977", "risk_allele": "T", "log_or": 0.113, "gene": "TAGAP", "eur_freq": 0.260},
            {"rsid": "rs689", "risk_allele": "T", "log_or": 0.405, "gene": "INS", "eur_freq": 0.730},
        ],
    },

    "rheumatoid_arthritis": {
        "name": "Rheumatoid Arthritis",
        "reference": "Okada et al. 2014 (Nature); Stahl et al. 2010 (Nature Genetics)",
        "snps": [
            {"rsid": "rs2476601", "risk_allele": "A", "log_or": 0.336, "gene": "PTPN22", "eur_freq": 0.100},
            {"rsid": "rs3184504", "risk_allele": "T", "log_or": 0.113, "gene": "SH2B3", "eur_freq": 0.480},
            {"rsid": "rs2104286", "risk_allele": "T", "log_or": 0.131, "gene": "IL2RA", "eur_freq": 0.250},
            {"rsid": "rs4810485", "risk_allele": "T", "log_or": 0.095, "gene": "CD40", "eur_freq": 0.260},
            {"rsid": "rs2228145", "risk_allele": "C", "log_or": 0.113, "gene": "IL6R", "eur_freq": 0.390},
            {"rsid": "rs874040", "risk_allele": "C", "log_or": 0.131, "gene": "RBPJ", "eur_freq": 0.380},
            {"rsid": "rs10488631", "risk_allele": "C", "log_or": 0.148, "gene": "IRF5", "eur_freq": 0.100},
            {"rsid": "rs3087243", "risk_allele": "G", "log_or": 0.148, "gene": "CTLA4", "eur_freq": 0.540},
            {"rsid": "rs7574865", "risk_allele": "T", "log_or": 0.182, "gene": "STAT4", "eur_freq": 0.230},
            {"rsid": "rs13031237", "risk_allele": "G", "log_or": 0.095, "gene": "REL", "eur_freq": 0.440},
            {"rsid": "rs2736340", "risk_allele": "T", "log_or": 0.105, "gene": "BLK", "eur_freq": 0.260},
            {"rsid": "rs6920220", "risk_allele": "A", "log_or": 0.131, "gene": "TNFAIP3", "eur_freq": 0.210},
        ],
    },

    "celiac_disease": {
        "name": "Celiac Disease",
        "reference": "Trynka et al. 2011 (Nature Genetics); Dubois et al. 2010",
        "snps": [
            {"rsid": "rs2187668", "risk_allele": "T", "log_or": 0.916, "gene": "HLA-DQ2.5 tag", "eur_freq": 0.130},
            {"rsid": "rs7454108", "risk_allele": "C", "log_or": 0.693, "gene": "HLA-DQ8 tag", "eur_freq": 0.090},
            {"rsid": "rs2476601", "risk_allele": "A", "log_or": 0.182, "gene": "PTPN22", "eur_freq": 0.100},
            {"rsid": "rs3184504", "risk_allele": "T", "log_or": 0.131, "gene": "SH2B3", "eur_freq": 0.480},
            {"rsid": "rs1990760", "risk_allele": "T", "log_or": 0.113, "gene": "IFIH1", "eur_freq": 0.610},
            {"rsid": "rs917997", "risk_allele": "A", "log_or": 0.131, "gene": "IL18RAP", "eur_freq": 0.240},
            {"rsid": "rs2816316", "risk_allele": "C", "log_or": 0.131, "gene": "RGS1", "eur_freq": 0.820},
            {"rsid": "rs13010713", "risk_allele": "A", "log_or": 0.113, "gene": "IL2/IL21", "eur_freq": 0.260},
            {"rsid": "rs1738074", "risk_allele": "T", "log_or": 0.095, "gene": "TAGAP", "eur_freq": 0.540},
            {"rsid": "rs3087243", "risk_allele": "G", "log_or": 0.113, "gene": "CTLA4", "eur_freq": 0.540},
            {"rsid": "rs6441961", "risk_allele": "C", "log_or": 0.148, "gene": "CCR3/CCR5", "eur_freq": 0.200},
            {"rsid": "rs802734", "risk_allele": "T", "log_or": 0.095, "gene": "CIITA/SOCS1", "eur_freq": 0.470},
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
    raw = 0.5 * (1.0 + math.erf(z / math.sqrt(2.0))) * 100
    return max(0.1, min(99.9, raw))


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
