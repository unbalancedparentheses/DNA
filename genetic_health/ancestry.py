"""Ancestry estimation from Ancestry-Informative Markers (AIMs).

Uses ~55 well-characterized AIMs with published allele frequency data
to estimate superpopulation proportions via maximum-likelihood scoring.

References:
  - Kosoy et al. 2009 (doi:10.1007/s00439-008-0585-5)
  - Kidd et al. 2014 (doi:10.1016/j.fsigen.2013.11.010)
  - 1000 Genomes Project (phase 3)
"""

import math

# =============================================================================
# ANCESTRY-INFORMATIVE MARKERS
# =============================================================================
# Superpopulations: EUR (European), AFR (African), EAS (East Asian),
# SAS (South Asian), AMR (Admixed American) — matching 1000 Genomes.
#
# Each entry:
#   allele: the informative allele to count (0/1/2 copies)
#   frequencies: population frequency of that allele

ANCESTRY_MARKERS = {
    # --- Pigmentation ---
    "rs1426654": {
        "gene": "SLC24A5",
        "description": "Skin pigmentation",
        "allele": "A",
        "frequencies": {"EUR": 0.978, "AFR": 0.073, "EAS": 0.013, "SAS": 0.600, "AMR": 0.553},
    },
    "rs16891982": {
        "gene": "SLC45A2",
        "description": "Skin/hair pigmentation",
        "allele": "G",
        "frequencies": {"EUR": 0.937, "AFR": 0.032, "EAS": 0.005, "SAS": 0.080, "AMR": 0.430},
    },
    "rs12913832": {
        "gene": "HERC2",
        "description": "Eye color (blue/brown)",
        "allele": "G",
        "frequencies": {"EUR": 0.723, "AFR": 0.033, "EAS": 0.010, "SAS": 0.050, "AMR": 0.280},
    },
    "rs1800407": {
        "gene": "OCA2",
        "description": "Eye color modifier",
        "allele": "A",
        "frequencies": {"EUR": 0.078, "AFR": 0.009, "EAS": 0.001, "SAS": 0.010, "AMR": 0.040},
    },
    "rs1393350": {
        "gene": "TYR",
        "description": "Tyrosinase — pigmentation",
        "allele": "A",
        "frequencies": {"EUR": 0.237, "AFR": 0.035, "EAS": 0.000, "SAS": 0.075, "AMR": 0.125},
    },
    "rs12821256": {
        "gene": "KITLG",
        "description": "Hair color (blond)",
        "allele": "C",
        "frequencies": {"EUR": 0.135, "AFR": 0.004, "EAS": 0.001, "SAS": 0.010, "AMR": 0.053},
    },
    "rs4959270": {
        "gene": "IRF4",
        "description": "Hair/skin pigmentation",
        "allele": "A",
        "frequencies": {"EUR": 0.480, "AFR": 0.155, "EAS": 0.045, "SAS": 0.200, "AMR": 0.320},
    },
    # --- DARC / Duffy ---
    "rs2814778": {
        "gene": "DARC",
        "description": "Duffy blood group (malaria resistance)",
        "allele": "C",
        "frequencies": {"EUR": 0.995, "AFR": 0.003, "EAS": 0.998, "SAS": 0.990, "AMR": 0.620},
    },
    # --- Morphology ---
    "rs3827760": {
        "gene": "EDAR",
        "description": "Hair thickness / shovel-shaped incisors",
        "allele": "A",
        "frequencies": {"EUR": 0.008, "AFR": 0.002, "EAS": 0.870, "SAS": 0.050, "AMR": 0.430},
    },
    # --- Alcohol metabolism ---
    "rs671": {
        "gene": "ALDH2",
        "description": "Alcohol flush reaction",
        "allele": "A",
        "frequencies": {"EUR": 0.001, "AFR": 0.001, "EAS": 0.190, "SAS": 0.008, "AMR": 0.020},
    },
    "rs1229984": {
        "gene": "ADH1B",
        "description": "Alcohol metabolism speed",
        "allele": "T",
        "frequencies": {"EUR": 0.044, "AFR": 0.020, "EAS": 0.700, "SAS": 0.150, "AMR": 0.120},
    },
    # --- Lactase persistence ---
    "rs4988235": {
        "gene": "MCM6/LCT",
        "description": "Lactase persistence",
        "allele": "A",
        "frequencies": {"EUR": 0.742, "AFR": 0.078, "EAS": 0.010, "SAS": 0.300, "AMR": 0.370},
    },
    # --- Drug metabolism ---
    "rs776746": {
        "gene": "CYP3A5",
        "description": "CYP3A5 expresser status",
        "allele": "T",
        "frequencies": {"EUR": 0.065, "AFR": 0.680, "EAS": 0.270, "SAS": 0.350, "AMR": 0.230},
    },
    "rs4149056": {
        "gene": "SLCO1B1",
        "description": "Statin transporter",
        "allele": "C",
        "frequencies": {"EUR": 0.155, "AFR": 0.020, "EAS": 0.120, "SAS": 0.060, "AMR": 0.090},
    },
    # --- Earwax / body odor ---
    "rs17822931": {
        "gene": "ABCC11",
        "description": "Earwax type / body odor",
        "allele": "T",
        "frequencies": {"EUR": 0.120, "AFR": 0.018, "EAS": 0.920, "SAS": 0.350, "AMR": 0.460},
    },
    # --- Immune / disease ---
    "rs8176719": {
        "gene": "ABO",
        "description": "Blood type O",
        "allele": "T",
        "frequencies": {"EUR": 0.630, "AFR": 0.700, "EAS": 0.550, "SAS": 0.590, "AMR": 0.680},
    },
    "rs12075": {
        "gene": "DARC",
        "description": "Duffy antigen receptor variant",
        "allele": "A",
        "frequencies": {"EUR": 0.427, "AFR": 0.870, "EAS": 0.040, "SAS": 0.270, "AMR": 0.420},
    },
    "rs2395029": {
        "gene": "HLA-B",
        "description": "HLA allele — HIV progression",
        "allele": "G",
        "frequencies": {"EUR": 0.056, "AFR": 0.003, "EAS": 0.001, "SAS": 0.010, "AMR": 0.025},
    },
    # --- Cardiovascular ---
    "rs1801133": {
        "gene": "MTHFR",
        "description": "Folate metabolism",
        "allele": "A",
        "frequencies": {"EUR": 0.340, "AFR": 0.090, "EAS": 0.290, "SAS": 0.160, "AMR": 0.450},
    },
    "rs6025": {
        "gene": "F5",
        "description": "Factor V Leiden",
        "allele": "A",
        "frequencies": {"EUR": 0.035, "AFR": 0.001, "EAS": 0.001, "SAS": 0.005, "AMR": 0.015},
    },
    "rs1799963": {
        "gene": "F2",
        "description": "Prothrombin G20210A",
        "allele": "A",
        "frequencies": {"EUR": 0.013, "AFR": 0.001, "EAS": 0.001, "SAS": 0.003, "AMR": 0.008},
    },
    # --- Taste / diet ---
    "rs713598": {
        "gene": "TAS2R38",
        "description": "Bitter taste perception (PTC)",
        "allele": "G",
        "frequencies": {"EUR": 0.460, "AFR": 0.600, "EAS": 0.350, "SAS": 0.500, "AMR": 0.480},
    },
    "rs1726866": {
        "gene": "TAS2R38",
        "description": "Bitter taste perception",
        "allele": "C",
        "frequencies": {"EUR": 0.460, "AFR": 0.600, "EAS": 0.350, "SAS": 0.500, "AMR": 0.480},
    },
    # --- Adiposity / metabolism ---
    "rs9939609": {
        "gene": "FTO",
        "description": "Obesity susceptibility",
        "allele": "A",
        "frequencies": {"EUR": 0.415, "AFR": 0.490, "EAS": 0.125, "SAS": 0.310, "AMR": 0.330},
    },
    "rs7903146": {
        "gene": "TCF7L2",
        "description": "Type 2 diabetes risk",
        "allele": "T",
        "frequencies": {"EUR": 0.290, "AFR": 0.290, "EAS": 0.035, "SAS": 0.260, "AMR": 0.270},
    },
    # --- APOE ---
    "rs429358": {
        "gene": "APOE",
        "description": "APOE e4 allele",
        "allele": "C",
        "frequencies": {"EUR": 0.153, "AFR": 0.265, "EAS": 0.091, "SAS": 0.110, "AMR": 0.120},
    },
    "rs7412": {
        "gene": "APOE",
        "description": "APOE e2 allele",
        "allele": "T",
        "frequencies": {"EUR": 0.076, "AFR": 0.108, "EAS": 0.090, "SAS": 0.050, "AMR": 0.055},
    },
    # --- Inflammatory ---
    "rs1800629": {
        "gene": "TNF",
        "description": "TNF-alpha promoter",
        "allele": "A",
        "frequencies": {"EUR": 0.135, "AFR": 0.112, "EAS": 0.020, "SAS": 0.050, "AMR": 0.065},
    },
    "rs1800896": {
        "gene": "IL10",
        "description": "IL-10 anti-inflammatory cytokine",
        "allele": "A",
        "frequencies": {"EUR": 0.465, "AFR": 0.290, "EAS": 0.380, "SAS": 0.410, "AMR": 0.380},
    },
    # --- Miscellaneous highly informative ---
    "rs2065160": {
        "gene": "LOXL1",
        "description": "Exfoliation glaucoma",
        "allele": "C",
        "frequencies": {"EUR": 0.285, "AFR": 0.900, "EAS": 0.100, "SAS": 0.180, "AMR": 0.400},
    },
    "rs4331426": {
        "gene": "LRRK2",
        "description": "Leprosy susceptibility (GWAS)",
        "allele": "G",
        "frequencies": {"EUR": 0.070, "AFR": 0.280, "EAS": 0.470, "SAS": 0.380, "AMR": 0.180},
    },
    "rs2250072": {
        "gene": "HBB",
        "description": "Beta-globin region",
        "allele": "T",
        "frequencies": {"EUR": 0.180, "AFR": 0.450, "EAS": 0.050, "SAS": 0.230, "AMR": 0.260},
    },
    "rs310644": {
        "gene": "GJB2",
        "description": "Connexin 26 region",
        "allele": "G",
        "frequencies": {"EUR": 0.850, "AFR": 0.510, "EAS": 0.760, "SAS": 0.700, "AMR": 0.680},
    },
    "rs3811801": {
        "gene": "TRPV6",
        "description": "Calcium absorption channel",
        "allele": "G",
        "frequencies": {"EUR": 0.830, "AFR": 0.030, "EAS": 0.710, "SAS": 0.600, "AMR": 0.520},
    },
    "rs1042602": {
        "gene": "TYR",
        "description": "Tyrosinase S192Y",
        "allele": "A",
        "frequencies": {"EUR": 0.380, "AFR": 0.005, "EAS": 0.001, "SAS": 0.050, "AMR": 0.190},
    },
    "rs2228479": {
        "gene": "MC1R",
        "description": "Melanocortin 1 receptor",
        "allele": "A",
        "frequencies": {"EUR": 0.095, "AFR": 0.010, "EAS": 0.230, "SAS": 0.040, "AMR": 0.070},
    },
    "rs1805007": {
        "gene": "MC1R",
        "description": "Red hair / fair skin",
        "allele": "T",
        "frequencies": {"EUR": 0.098, "AFR": 0.001, "EAS": 0.001, "SAS": 0.003, "AMR": 0.035},
    },
    "rs2402130": {
        "gene": "SLC24A4",
        "description": "Hair/eye color",
        "allele": "A",
        "frequencies": {"EUR": 0.530, "AFR": 0.140, "EAS": 0.210, "SAS": 0.290, "AMR": 0.370},
    },
    "rs1800414": {
        "gene": "OCA2",
        "description": "Skin pigmentation (East Asian)",
        "allele": "C",
        "frequencies": {"EUR": 0.001, "AFR": 0.001, "EAS": 0.630, "SAS": 0.060, "AMR": 0.100},
    },
    "rs6497268": {
        "gene": "MFSD12",
        "description": "Skin pigmentation (African)",
        "allele": "T",
        "frequencies": {"EUR": 0.070, "AFR": 0.530, "EAS": 0.040, "SAS": 0.090, "AMR": 0.180},
    },
    "rs1834640": {
        "gene": "SLC24A5",
        "description": "Skin lightening linked variant",
        "allele": "A",
        "frequencies": {"EUR": 0.980, "AFR": 0.060, "EAS": 0.020, "SAS": 0.580, "AMR": 0.540},
    },
    "rs174537": {
        "gene": "FADS1",
        "description": "Fatty acid desaturase",
        "allele": "T",
        "frequencies": {"EUR": 0.660, "AFR": 0.960, "EAS": 0.550, "SAS": 0.610, "AMR": 0.710},
    },
    "rs2187668": {
        "gene": "HLA-DQA1",
        "description": "Celiac disease risk",
        "allele": "T",
        "frequencies": {"EUR": 0.130, "AFR": 0.025, "EAS": 0.020, "SAS": 0.060, "AMR": 0.070},
    },
    "rs4833103": {
        "gene": "NRG1",
        "description": "Neuregulin 1 (ancestry informative)",
        "allele": "T",
        "frequencies": {"EUR": 0.550, "AFR": 0.120, "EAS": 0.780, "SAS": 0.400, "AMR": 0.350},
    },
    "rs7657799": {
        "gene": "DDB1",
        "description": "Skin pigmentation (DDB1)",
        "allele": "T",
        "frequencies": {"EUR": 0.880, "AFR": 0.220, "EAS": 0.960, "SAS": 0.780, "AMR": 0.600},
    },
    "rs260714": {
        "gene": "LILRA3",
        "description": "Immune receptor deletion",
        "allele": "C",
        "frequencies": {"EUR": 0.680, "AFR": 0.400, "EAS": 0.780, "SAS": 0.580, "AMR": 0.560},
    },
    "rs3796332": {
        "gene": "ATRN",
        "description": "Attractin (pigmentation modulator)",
        "allele": "A",
        "frequencies": {"EUR": 0.170, "AFR": 0.510, "EAS": 0.050, "SAS": 0.200, "AMR": 0.250},
    },
    "rs3916235": {
        "gene": "ASPM",
        "description": "Brain size gene",
        "allele": "A",
        "frequencies": {"EUR": 0.460, "AFR": 0.115, "EAS": 0.380, "SAS": 0.340, "AMR": 0.330},
    },
    "rs4680": {
        "gene": "COMT",
        "description": "Catechol-O-methyltransferase",
        "allele": "A",
        "frequencies": {"EUR": 0.498, "AFR": 0.330, "EAS": 0.280, "SAS": 0.370, "AMR": 0.400},
    },
    "rs762551": {
        "gene": "CYP1A2",
        "description": "Caffeine metabolism",
        "allele": "A",
        "frequencies": {"EUR": 0.680, "AFR": 0.520, "EAS": 0.390, "SAS": 0.580, "AMR": 0.580},
    },
    "rs1800562": {
        "gene": "HFE",
        "description": "Hemochromatosis C282Y",
        "allele": "A",
        "frequencies": {"EUR": 0.058, "AFR": 0.001, "EAS": 0.001, "SAS": 0.003, "AMR": 0.020},
    },
    "rs1799971": {
        "gene": "OPRM1",
        "description": "Opioid receptor",
        "allele": "G",
        "frequencies": {"EUR": 0.153, "AFR": 0.029, "EAS": 0.388, "SAS": 0.220, "AMR": 0.130},
    },
    "rs4532": {
        "gene": "DRD1",
        "description": "Dopamine receptor D1",
        "allele": "G",
        "frequencies": {"EUR": 0.565, "AFR": 0.340, "EAS": 0.730, "SAS": 0.500, "AMR": 0.480},
    },
    "rs806377": {
        "gene": "CNR1",
        "description": "Cannabinoid receptor",
        "allele": "A",
        "frequencies": {"EUR": 0.380, "AFR": 0.180, "EAS": 0.410, "SAS": 0.350, "AMR": 0.310},
    },
}

POPULATIONS = ["EUR", "AFR", "EAS", "SAS", "AMR"]

POPULATION_LABELS = {
    "EUR": "European",
    "AFR": "African",
    "EAS": "East Asian",
    "SAS": "South Asian",
    "AMR": "Admixed American",
}


# =============================================================================
# POPULATION-SPECIFIC WARNINGS
# =============================================================================

POPULATION_NOTES = {
    ("ALDH2", "reduced"): [
        "ALDH2*2 is almost exclusively East Asian (30-40% frequency). "
        "If detected in a non-EAS individual, verify the genotype."
    ],
    ("ALDH2", "non_functional"): [
        "Homozygous ALDH2*2 causes severe alcohol intolerance. "
        "Essentially absent outside East Asia."
    ],
    ("ADH1B", "fast"): [
        "ADH1B*2 (fast metabolizer) is common in East Asians (~70%) "
        "but rare in Europeans (~4%). Frequency varies greatly by population."
    ],
    ("SLC24A5", "derived"): [
        "The SLC24A5 A allele (light skin) is near-fixed in Europeans "
        "but also moderately common in South Asians (~60%)."
    ],
    ("HFE", "homozygous"): [
        "HFE C282Y homozygosity is predominantly European (~1 in 200). "
        "Extremely rare in non-European populations."
    ],
    ("HFE", "carrier"): [
        "HFE C282Y carrier frequency is ~10% in Northern Europeans, "
        "near-zero in African and East Asian populations."
    ],
    ("F5", "heterozygous"): [
        "Factor V Leiden is primarily European (~3-5% carrier frequency). "
        "Very rare in African and East Asian populations."
    ],
    ("CYP3A5", "expresser"): [
        "CYP3A5*1 (expresser) is the global majority allele — ~70% in Africans, "
        "~27% in East Asians, but only ~6% in Europeans."
    ],
    ("CYP2D6", "poor"): [
        "CYP2D6 poor metabolizer frequency varies: ~5-10% European, "
        "~1-2% East Asian, ~1-3% African. Genotype interpretation may "
        "depend on population-specific allele frequencies."
    ],
    ("MTHFR", "reduced"): [
        "MTHFR C677T TT frequency varies: ~10% European, ~25% Mexican/Hispanic, "
        "~1-2% African. Clinical impact may differ by folate intake."
    ],
    ("G6PD", "deficient"): [
        "G6PD deficiency is much more common in African, Mediterranean, "
        "and Southeast Asian populations (malaria protection)."
    ],
    ("HBB", "carrier"): [
        "Sickle cell trait frequency: ~8% African American, ~15-25% West African, "
        "rare in Europeans and East Asians."
    ],
    ("CFTR", "carrier"): [
        "CFTR carrier frequency is ~4% in European-descent populations, "
        "much lower in African and East Asian populations."
    ],
    ("APOE", "e4_homozygous"): [
        "APOE e4 frequency varies: ~14% European, ~26% African, ~9% East Asian. "
        "Alzheimer's risk conferred by e4 may differ by ancestry."
    ],
    ("APOE", "e4_heterozygous"): [
        "APOE e4 frequency varies substantially by population. "
        "Risk interpretation is most validated in European-descent studies."
    ],
    ("BRCA1", "pathogenic"): [
        "Certain BRCA1/2 founder mutations are more common in specific populations "
        "(e.g., Ashkenazi Jewish). Panel testing may be more informative."
    ],
    ("BRCA2", "pathogenic"): [
        "Certain BRCA1/2 founder mutations are more common in specific populations "
        "(e.g., Ashkenazi Jewish). Panel testing may be more informative."
    ],
    ("UGT1A1", "reduced"): [
        "UGT1A1*28 (Gilbert syndrome) frequency varies: ~12% European, "
        "~35% African, ~2% East Asian. Drug dosing recommendations may differ."
    ],
    ("DPYD", "reduced"): [
        "DPYD variant frequencies differ by population. "
        "DPYD*2A is ~1% in Europeans, rarer in other populations."
    ],
    ("TPMT", "intermediate"): [
        "TPMT poor metabolizer variants are more common in Europeans (~5%) "
        "than in East Asians (~2%) or Africans (~3%)."
    ],
}


def get_population_warnings(gene: str, status: str) -> list:
    """Return population-specific warnings for a gene/status combination."""
    return list(POPULATION_NOTES.get((gene, status), []))


# =============================================================================
# ANCESTRY ESTIMATION
# =============================================================================

def _count_allele(genotype: str, informative_allele: str) -> int:
    """Count copies of informative allele in genotype (0, 1, or 2)."""
    return sum(1 for a in genotype if a == informative_allele)


def _marker_informativeness(frequencies: dict) -> float:
    """Calculate informativeness of a marker using Fst-like delta metric.

    Higher values = more discriminating between populations.
    Returns the max absolute frequency difference between any two populations.
    """
    freqs = list(frequencies.values())
    return max(freqs) - min(freqs)


def _softmax(log_likelihoods: dict) -> dict:
    """Softmax normalization of log-likelihoods to proportions."""
    max_ll = max(log_likelihoods.values())
    exps = {pop: math.exp(ll - max_ll) for pop, ll in log_likelihoods.items()}
    total = sum(exps.values())
    return {pop: val / total for pop, val in exps.items()}


def estimate_ancestry(genome_by_rsid: dict) -> dict:
    """Estimate ancestry proportions from genome data.

    Algorithm:
        1. For each AIM present in genome, count informative allele copies (0/1/2)
        2. Per population: sum log-likelihood = n*log(freq) + (2-n)*log(1-freq)
           with freq clamped to [0.001, 0.999]
        3. Softmax normalize to proportions
        4. Confidence based on marker count

    Args:
        genome_by_rsid: Dict mapping rsID -> {genotype, chromosome, position}

    Returns:
        Dict with keys: proportions, markers_found, confidence,
        top_ancestry, details
    """
    log_likelihoods = {pop: 0.0 for pop in POPULATIONS}
    details = []
    markers_found = 0
    total_informativeness = 0.0

    for rsid, marker in ANCESTRY_MARKERS.items():
        if rsid not in genome_by_rsid:
            continue

        genotype = genome_by_rsid[rsid]["genotype"]
        allele = marker["allele"]
        n = _count_allele(genotype, allele)
        markers_found += 1
        info_score = _marker_informativeness(marker["frequencies"])
        total_informativeness += info_score

        marker_detail = {
            "rsid": rsid,
            "gene": marker["gene"],
            "description": marker["description"],
            "genotype": genotype,
            "allele_count": n,
        }

        for pop in POPULATIONS:
            freq = marker["frequencies"][pop]
            # Clamp to avoid log(0)
            freq = max(0.001, min(0.999, freq))
            ll = n * math.log(freq) + (2 - n) * math.log(1 - freq)
            log_likelihoods[pop] += ll

        details.append(marker_detail)

    if markers_found == 0:
        return {
            "proportions": {pop: 1.0 / len(POPULATIONS) for pop in POPULATIONS},
            "markers_found": 0,
            "confidence": "none",
            "top_ancestry": "Unknown",
            "details": [],
        }

    proportions = _softmax(log_likelihoods)

    # Determine confidence using both marker count and total informativeness
    # Average informativeness per marker (scale 0-1; >0.5 is highly informative)
    avg_info = total_informativeness / markers_found if markers_found else 0
    # Weighted effective marker count (markers * avg discriminative power)
    effective_markers = markers_found * avg_info

    if markers_found >= 40 and effective_markers >= 15:
        confidence = "high"
    elif markers_found >= 20 and effective_markers >= 8:
        confidence = "moderate"
    else:
        confidence = "low"

    # Top ancestry
    top_pop = max(proportions, key=proportions.get)
    top_ancestry = POPULATION_LABELS.get(top_pop, top_pop)

    return {
        "proportions": proportions,
        "markers_found": markers_found,
        "confidence": confidence,
        "top_ancestry": top_ancestry,
        "details": details,
    }
