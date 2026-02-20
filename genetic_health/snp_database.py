"""
Comprehensive SNP Database for Total Health Optimization
Covers: Drug metabolism, methylation, fitness, nutrition, sleep, cardiovascular,
cognition, longevity, inflammation, and lifestyle factors.

The optional ``freq`` field on some entries contains **minor (variant) allele
frequencies** per 1000 Genomes superpopulation (EUR, AFR, EAS, SAS, AMR).
Values represent the proportion of chromosomes carrying the variant allele in
that population.  They do NOT sum to 1.0 across populations — each population
is an independent measurement.
"""

COMPREHENSIVE_SNPS = {

    # =========================================================================
    # SECTION 1: DRUG METABOLISM (from original)
    # =========================================================================

    "rs762551": {
        "gene": "CYP1A2", "category": "Drug Metabolism",
        "variants": {
            "AA": {"status": "fast", "desc": "Fast caffeine metabolizer - clears caffeine quickly, lower cardiovascular risk from coffee", "magnitude": 1},
            "AC": {"status": "intermediate", "desc": "Intermediate caffeine metabolizer - moderate clearance, ~5-6hr half-life", "magnitude": 2},
            "CC": {"status": "slow", "desc": "Slow caffeine metabolizer - caffeine lingers 8-12hrs, increased cardiovascular risk with high intake", "magnitude": 3},
        }
    },
    "rs4244285": {
        "gene": "CYP2C19", "category": "Drug Metabolism",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal CYP2C19 - standard drug metabolism", "magnitude": 0},
            "GA": {"status": "intermediate", "desc": "Intermediate CYP2C19 (*2 carrier) - clopidogrel less effective", "magnitude": 3},
            "AA": {"status": "poor", "desc": "Poor CYP2C19 (*2/*2) - clopidogrel ineffective, use alternative antiplatelet", "magnitude": 4},
        },
        "freq": {"EUR": 0.15, "AFR": 0.18, "EAS": 0.30, "SAS": 0.34, "AMR": 0.12},
    },
    "rs12248560": {
        "gene": "CYP2C19", "category": "Drug Metabolism",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal CYP2C19 metabolism", "magnitude": 0},
            "CT": {"status": "rapid", "desc": "Rapid CYP2C19 (*17) - faster metabolism of PPIs, some antidepressants, may need higher doses", "magnitude": 2},
            "TT": {"status": "ultrarapid", "desc": "Ultrarapid CYP2C19 (*17/*17) - significantly faster drug metabolism", "magnitude": 3},
        }
    },
    "rs1799853": {
        "gene": "CYP2C9", "category": "Drug Metabolism",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal CYP2C9 - standard warfarin/NSAID metabolism", "magnitude": 0},
            "CT": {"status": "intermediate", "desc": "Intermediate CYP2C9 (*2) - warfarin dose reduction needed", "magnitude": 3},
            "TT": {"status": "poor", "desc": "Poor CYP2C9 (*2/*2) - significant warfarin sensitivity", "magnitude": 4},
        }
    },
    "rs1057910": {
        "gene": "CYP2C9", "category": "Drug Metabolism",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal CYP2C9 function", "magnitude": 0},
            "AC": {"status": "intermediate", "desc": "Intermediate CYP2C9 (*3) - warfarin dose reduction", "magnitude": 3},
            "CC": {"status": "poor", "desc": "Poor CYP2C9 (*3/*3) - high warfarin sensitivity", "magnitude": 4},
        }
    },
    "rs9923231": {
        "gene": "VKORC1", "category": "Drug Metabolism",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal warfarin sensitivity", "magnitude": 0},
            "GA": {"status": "sensitive", "desc": "Increased warfarin sensitivity - lower doses needed", "magnitude": 3},
            "AG": {"status": "sensitive", "desc": "Increased warfarin sensitivity - lower doses needed", "magnitude": 3},
            "AA": {"status": "highly_sensitive", "desc": "Highly warfarin sensitive - significantly lower doses", "magnitude": 4},
        },
        "freq": {"EUR": 0.39, "AFR": 0.11, "EAS": 0.90, "SAS": 0.34, "AMR": 0.45},
    },
    "rs4149056": {
        "gene": "SLCO1B1", "category": "Drug Metabolism",
        "variants": {
            "TT": {"status": "normal", "desc": "Normal statin transport - standard myopathy risk", "magnitude": 0},
            "TC": {"status": "intermediate", "desc": "Intermediate statin transporter - 4x myopathy risk with simvastatin", "magnitude": 3},
            "CT": {"status": "intermediate", "desc": "Intermediate statin transporter - 4x myopathy risk with simvastatin", "magnitude": 3},
            "CC": {"status": "poor", "desc": "Poor statin transporter - 17x myopathy risk, avoid simvastatin", "magnitude": 4},
        },
        "freq": {"EUR": 0.15, "AFR": 0.02, "EAS": 0.11, "SAS": 0.05, "AMR": 0.08},
    },
    "rs3892097": {
        "gene": "CYP2D6", "category": "Drug Metabolism",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal CYP2D6 - standard codeine/tramadol metabolism", "magnitude": 0},
            "GA": {"status": "intermediate", "desc": "Intermediate CYP2D6 - reduced opioid activation", "magnitude": 3},
            "AG": {"status": "intermediate", "desc": "Intermediate CYP2D6 - reduced opioid activation", "magnitude": 3},
            "AA": {"status": "poor", "desc": "Poor CYP2D6 (*4/*4) - codeine ineffective, tramadol reduced", "magnitude": 4},
        }
    },
    "rs776746": {
        "gene": "CYP3A5", "category": "Drug Metabolism",
        "variants": {
            "TT": {"status": "expressor", "desc": "CYP3A5 expressor - may need higher tacrolimus doses", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "Intermediate CYP3A5 expression", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "Intermediate CYP3A5 expression", "magnitude": 1},
            "CC": {"status": "non_expressor", "desc": "CYP3A5 non-expressor (*3/*3) - standard tacrolimus dosing", "magnitude": 1},
        }
    },
    "rs3918290": {
        "gene": "DPYD", "category": "Drug Metabolism",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal DPYD - standard fluoropyrimidine tolerance", "magnitude": 0},
            "CT": {"status": "intermediate", "desc": "Reduced DPYD - 50% dose reduction for 5-FU/capecitabine", "magnitude": 5},
            "TT": {"status": "deficient", "desc": "DPYD deficient - fluoropyrimidines contraindicated (can be fatal)", "magnitude": 6},
        },
        "freq": {"EUR": 0.01, "AFR": 0.001, "EAS": 0.001, "SAS": 0.001, "AMR": 0.005},
    },
    "rs1800460": {
        "gene": "TPMT", "category": "Drug Metabolism",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal TPMT - standard thiopurine tolerance", "magnitude": 0},
            "CT": {"status": "intermediate", "desc": "Intermediate TPMT - thiopurine dose reduction needed", "magnitude": 4},
            "TC": {"status": "intermediate", "desc": "Intermediate TPMT - thiopurine dose reduction needed", "magnitude": 4},
            "TT": {"status": "poor", "desc": "Poor TPMT - thiopurines can cause severe myelosuppression", "magnitude": 5},
        }
    },
    "rs2395029": {
        "gene": "HLA-B", "category": "Drug Metabolism",
        "variants": {
            "TT": {"status": "normal", "desc": "Low abacavir hypersensitivity risk", "magnitude": 0},
            "TG": {"status": "carrier", "desc": "HLA-B*5701 carrier - abacavir contraindicated", "magnitude": 5},
            "GT": {"status": "carrier", "desc": "HLA-B*5701 carrier - abacavir contraindicated", "magnitude": 5},
            "GG": {"status": "positive", "desc": "HLA-B*5701 positive - abacavir contraindicated", "magnitude": 5},
        }
    },

    # =========================================================================
    # SECTION 2: METHYLATION & DETOXIFICATION
    # =========================================================================

    "rs1801133": {
        "gene": "MTHFR", "category": "Methylation",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal MTHFR C677 - full methylation capacity", "magnitude": 0},
            "AG": {"status": "reduced", "desc": "MTHFR C677T heterozygous - ~35% reduced activity, may benefit from methylfolate", "magnitude": 2},
            "GA": {"status": "reduced", "desc": "MTHFR C677T heterozygous - ~35% reduced activity, may benefit from methylfolate", "magnitude": 2},
            "AA": {"status": "significantly_reduced", "desc": "MTHFR C677T homozygous - ~70% reduced activity, methylfolate recommended over folic acid", "magnitude": 3},
        },
        "freq": {"EUR": 0.36, "AFR": 0.11, "EAS": 0.33, "SAS": 0.15, "AMR": 0.46},
    },
    "rs1801131": {
        "gene": "MTHFR", "category": "Methylation",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal MTHFR A1298 function", "magnitude": 0},
            "AC": {"status": "reduced", "desc": "MTHFR A1298C heterozygous - mild reduction", "magnitude": 1},
            "CA": {"status": "reduced", "desc": "MTHFR A1298C heterozygous - mild reduction", "magnitude": 1},
            "CC": {"status": "reduced", "desc": "MTHFR A1298C homozygous - moderate reduction in BH4 recycling", "magnitude": 2},
            "TT": {"status": "normal", "desc": "Normal MTHFR A1298 function (23andMe orientation)", "magnitude": 0},
            "TG": {"status": "reduced", "desc": "MTHFR A1298C heterozygous (23andMe orientation)", "magnitude": 1},
            "GT": {"status": "reduced", "desc": "MTHFR A1298C heterozygous (23andMe orientation)", "magnitude": 1},
            "GG": {"status": "reduced", "desc": "MTHFR A1298C homozygous (23andMe orientation)", "magnitude": 2},
        }
    },
    "rs1805087": {
        "gene": "MTR", "category": "Methylation",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal methionine synthase function", "magnitude": 0},
            "AG": {"status": "reduced", "desc": "MTR A2756G heterozygous - reduced B12 utilization", "magnitude": 1},
            "GA": {"status": "reduced", "desc": "MTR A2756G heterozygous - reduced B12 utilization", "magnitude": 1},
            "GG": {"status": "significantly_reduced", "desc": "MTR A2756G homozygous - may need higher B12, check homocysteine", "magnitude": 2},
        }
    },
    "rs1801394": {
        "gene": "MTRR", "category": "Methylation",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal MTRR function - efficient B12 recycling", "magnitude": 0},
            "AG": {"status": "reduced", "desc": "MTRR A66G heterozygous - reduced B12 regeneration", "magnitude": 1},
            "GA": {"status": "reduced", "desc": "MTRR A66G heterozygous - reduced B12 regeneration", "magnitude": 1},
            "GG": {"status": "significantly_reduced", "desc": "MTRR A66G homozygous - impaired B12 recycling, consider methylcobalamin", "magnitude": 2},
        }
    },
    "rs234706": {
        "gene": "CBS", "category": "Methylation",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal CBS function - standard homocysteine processing", "magnitude": 0},
            "GA": {"status": "upregulated", "desc": "CBS C699T heterozygous - may have faster transsulfuration", "magnitude": 1},
            "AG": {"status": "upregulated", "desc": "CBS C699T heterozygous - may have faster transsulfuration", "magnitude": 1},
            "AA": {"status": "upregulated", "desc": "CBS C699T homozygous - upregulated CBS, may deplete homocysteine/SAMe faster", "magnitude": 2},
        }
    },
    "rs7946": {
        "gene": "PEMT", "category": "Methylation",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal PEMT - adequate choline synthesis", "magnitude": 0},
            "CT": {"status": "reduced", "desc": "PEMT G5765A heterozygous - may need more dietary choline", "magnitude": 1},
            "TC": {"status": "reduced", "desc": "PEMT G5765A heterozygous - may need more dietary choline", "magnitude": 1},
            "TT": {"status": "significantly_reduced", "desc": "PEMT G5765A homozygous - higher choline requirements, especially for women", "magnitude": 2},
        }
    },

    # NAT2 - Acetylation
    "rs1801280": {
        "gene": "NAT2", "category": "Detoxification",
        "variants": {
            "TT": {"status": "fast", "desc": "Fast NAT2 acetylator", "magnitude": 0},
            "TC": {"status": "intermediate", "desc": "Intermediate NAT2 acetylator", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "Intermediate NAT2 acetylator", "magnitude": 1},
            "CC": {"status": "slow", "desc": "Slow NAT2 acetylator - increased drug/toxin sensitivity", "magnitude": 2},
        }
    },
    "rs1799930": {
        "gene": "NAT2", "category": "Detoxification",
        "variants": {
            "GG": {"status": "fast", "desc": "Fast acetylator at this locus", "magnitude": 0},
            "GA": {"status": "intermediate", "desc": "Intermediate acetylator", "magnitude": 1},
            "AG": {"status": "intermediate", "desc": "Intermediate acetylator", "magnitude": 1},
            "AA": {"status": "slow", "desc": "Slow acetylator at this locus", "magnitude": 2},
        }
    },

    # Glutathione
    "rs1695": {
        "gene": "GSTP1", "category": "Detoxification",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal GSTP1 - good glutathione conjugation", "magnitude": 0},
            "AG": {"status": "reduced", "desc": "GSTP1 Ile105Val heterozygous - reduced detox capacity", "magnitude": 1},
            "GA": {"status": "reduced", "desc": "GSTP1 Ile105Val heterozygous - reduced detox capacity", "magnitude": 1},
            "GG": {"status": "significantly_reduced", "desc": "GSTP1 Val/Val - reduced glutathione conjugation, may benefit from NAC/glutathione support", "magnitude": 2},
        }
    },
    "rs1138272": {
        "gene": "GSTP1", "category": "Detoxification",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal GSTP1 Ala114 function", "magnitude": 0},
            "CT": {"status": "reduced", "desc": "GSTP1 Ala114Val heterozygous", "magnitude": 1},
            "TC": {"status": "reduced", "desc": "GSTP1 Ala114Val heterozygous", "magnitude": 1},
            "TT": {"status": "significantly_reduced", "desc": "GSTP1 Val/Val at 114 - further reduced activity", "magnitude": 2},
        }
    },
    "rs4880": {
        "gene": "SOD2", "category": "Detoxification",
        "variants": {
            "AA": {"status": "high_activity", "desc": "High SOD2 activity (Ala/Ala) - efficient mitochondrial antioxidant", "magnitude": 1},
            "AG": {"status": "intermediate", "desc": "Intermediate SOD2 activity", "magnitude": 0},
            "GA": {"status": "intermediate", "desc": "Intermediate SOD2 activity", "magnitude": 0},
            "GG": {"status": "low_activity", "desc": "Lower SOD2 activity (Val/Val) - may benefit from antioxidant support", "magnitude": 2},
        }
    },

    # =========================================================================
    # SECTION 3: NEUROTRANSMITTERS & COGNITION
    # =========================================================================

    "rs4680": {
        "gene": "COMT", "category": "Neurotransmitters",
        "variants": {
            "GG": {"status": "fast", "desc": "Fast COMT (Val/Val) - clears dopamine quickly, better stress resilience, may need more stimulation", "magnitude": 2},
            "AG": {"status": "intermediate", "desc": "Intermediate COMT (Val/Met) - balanced dopamine clearance", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "Intermediate COMT (Val/Met) - balanced dopamine clearance", "magnitude": 1},
            "AA": {"status": "slow", "desc": "Slow COMT (Met/Met) - higher dopamine, better working memory but more stress-sensitive, stimulants hit harder", "magnitude": 3},
        }
    },
    "rs4633": {
        "gene": "COMT", "category": "Neurotransmitters",
        "variants": {
            "CC": {"status": "fast", "desc": "Fast COMT haplotype marker", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "Intermediate COMT", "magnitude": 0},
            "TC": {"status": "intermediate", "desc": "Intermediate COMT", "magnitude": 0},
            "TT": {"status": "slow", "desc": "Slow COMT haplotype marker", "magnitude": 1},
        }
    },
    "rs6265": {
        "gene": "BDNF", "category": "Neurotransmitters",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal BDNF Val66 - standard neuroplasticity and memory", "magnitude": 0},
            "CT": {"status": "reduced", "desc": "BDNF Val66Met heterozygous - reduced activity-dependent BDNF secretion, exercise especially beneficial", "magnitude": 2},
            "TC": {"status": "reduced", "desc": "BDNF Val66Met heterozygous - reduced activity-dependent BDNF secretion, exercise especially beneficial", "magnitude": 2},
            "TT": {"status": "significantly_reduced", "desc": "BDNF Met/Met - reduced neuroplasticity, higher depression risk, exercise strongly recommended", "magnitude": 3},
        }
    },
    "rs25531": {
        "gene": "SLC6A4", "category": "Neurotransmitters",
        "variants": {
            "AA": {"status": "low_expression", "desc": "Low serotonin transporter expression (La/La) - may be more responsive to SSRIs", "magnitude": 1},
            "AG": {"status": "intermediate", "desc": "Intermediate serotonin transporter", "magnitude": 0},
            "GA": {"status": "intermediate", "desc": "Intermediate serotonin transporter", "magnitude": 0},
            "GG": {"status": "high_expression", "desc": "Higher serotonin transporter expression - may need higher SSRI doses", "magnitude": 1},
        }
    },
    "rs1800497": {
        "gene": "ANKK1/DRD2", "category": "Neurotransmitters",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal D2 receptor density", "magnitude": 0},
            "CT": {"status": "reduced", "desc": "Taq1A heterozygous - reduced D2 receptors, may seek more stimulation/rewards", "magnitude": 2},
            "TC": {"status": "reduced", "desc": "Taq1A heterozygous - reduced D2 receptors, may seek more stimulation/rewards", "magnitude": 2},
            "TT": {"status": "significantly_reduced", "desc": "Taq1A homozygous - ~40% fewer D2 receptors, higher addiction susceptibility", "magnitude": 3},
        }
    },
    "rs1799971": {
        "gene": "OPRM1", "category": "Neurotransmitters",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal mu-opioid receptor function", "magnitude": 0},
            "AG": {"status": "altered", "desc": "OPRM1 A118G heterozygous - altered opioid/alcohol response", "magnitude": 2},
            "GA": {"status": "altered", "desc": "OPRM1 A118G heterozygous - altered opioid/alcohol response", "magnitude": 2},
            "GG": {"status": "significantly_altered", "desc": "OPRM1 G/G - reduced opioid sensitivity, may need higher doses for pain, altered alcohol reward", "magnitude": 3},
        }
    },

    # =========================================================================
    # SECTION 4: CAFFEINE & ADENOSINE
    # =========================================================================

    "rs5751876": {
        "gene": "ADORA2A", "category": "Caffeine Response",
        "variants": {
            "CC": {"status": "lower_sensitivity", "desc": "Lower caffeine sensitivity - less anxiety from caffeine", "magnitude": 1},
            "CT": {"status": "normal", "desc": "Normal caffeine sensitivity", "magnitude": 0},
            "TC": {"status": "normal", "desc": "Normal caffeine sensitivity", "magnitude": 0},
            "TT": {"status": "high_sensitivity", "desc": "High caffeine sensitivity - more prone to caffeine anxiety and sleep disruption", "magnitude": 2},
        }
    },
    "rs2298383": {
        "gene": "ADORA2A", "category": "Caffeine Response",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal anxiety response to caffeine", "magnitude": 0},
            "CT": {"status": "intermediate", "desc": "Intermediate caffeine-anxiety response", "magnitude": 1},
            "TC": {"status": "intermediate", "desc": "Intermediate caffeine-anxiety response", "magnitude": 1},
            "TT": {"status": "anxiety_prone", "desc": "Increased anxiety response to caffeine - consider lower doses or alternatives", "magnitude": 2},
        }
    },
    "rs73598374": {
        "gene": "ADA", "category": "Caffeine Response",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal adenosine deaminase - standard sleep pressure", "magnitude": 0},
            "CT": {"status": "reduced", "desc": "ADA G22A heterozygous - slower adenosine clearance, deeper sleep need", "magnitude": 1},
            "TC": {"status": "reduced", "desc": "ADA G22A heterozygous - slower adenosine clearance, deeper sleep need", "magnitude": 1},
            "TT": {"status": "significantly_reduced", "desc": "ADA G22A homozygous - high sleep pressure, may need more sleep", "magnitude": 2},
        }
    },

    # =========================================================================
    # SECTION 5: SLEEP & CIRCADIAN RHYTHM
    # =========================================================================

    "rs1801260": {
        "gene": "CLOCK", "category": "Sleep/Circadian",
        "variants": {
            "TT": {"status": "normal", "desc": "Normal CLOCK gene - standard circadian timing", "magnitude": 0},
            "TC": {"status": "evening_tendency", "desc": "CLOCK 3111C heterozygous - slight evening preference", "magnitude": 1},
            "CT": {"status": "evening_tendency", "desc": "CLOCK 3111C heterozygous - slight evening preference", "magnitude": 1},
            "CC": {"status": "evening_type", "desc": "CLOCK 3111C homozygous - evening chronotype, may have delayed sleep phase", "magnitude": 2},
        }
    },
    "rs57875989": {
        "gene": "PER2", "category": "Sleep/Circadian",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal PER2 - standard circadian period", "magnitude": 0},
            "CG": {"status": "morning_tendency", "desc": "PER2 variant - tendency toward morning chronotype", "magnitude": 1},
            "GC": {"status": "morning_tendency", "desc": "PER2 variant - tendency toward morning chronotype", "magnitude": 1},
            "GG": {"status": "morning_type", "desc": "PER2 variant homozygous - strong morning chronotype", "magnitude": 2},
        }
    },
    "rs12649507": {
        "gene": "ARNTL", "category": "Sleep/Circadian",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal BMAL1 function", "magnitude": 0},
            "AG": {"status": "altered", "desc": "BMAL1 variant - may affect circadian amplitude", "magnitude": 1},
            "GA": {"status": "altered", "desc": "BMAL1 variant - may affect circadian amplitude", "magnitude": 1},
            "GG": {"status": "significantly_altered", "desc": "BMAL1 variant homozygous - may have weaker circadian rhythm", "magnitude": 2},
        }
    },
    "rs28532698": {
        "gene": "MTNR1B", "category": "Sleep/Circadian",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal melatonin receptor", "magnitude": 0},
            "CT": {"status": "altered", "desc": "MTNR1B variant - altered melatonin signaling, higher T2D risk with late eating", "magnitude": 2},
            "TC": {"status": "altered", "desc": "MTNR1B variant - altered melatonin signaling, higher T2D risk with late eating", "magnitude": 2},
            "TT": {"status": "significantly_altered", "desc": "MTNR1B variant homozygous - avoid late-night eating, higher diabetes risk", "magnitude": 3},
        }
    },

    # =========================================================================
    # SECTION 6: FITNESS & EXERCISE RESPONSE
    # =========================================================================

    "rs1815739": {
        "gene": "ACTN3", "category": "Fitness",
        "variants": {
            "CC": {"status": "power", "desc": "ACTN3 R/R (power type) - fast-twitch muscle fiber advantage, suited for sprinting/power sports", "magnitude": 2},
            "CT": {"status": "mixed", "desc": "ACTN3 R/X (mixed) - balanced muscle fiber composition", "magnitude": 1},
            "TC": {"status": "mixed", "desc": "ACTN3 R/X (mixed) - balanced muscle fiber composition", "magnitude": 1},
            "TT": {"status": "endurance", "desc": "ACTN3 X/X (endurance type) - no alpha-actinin-3, better suited for endurance sports", "magnitude": 2},
        },
        "freq": {"EUR": 0.42, "AFR": 0.10, "EAS": 0.52, "SAS": 0.38, "AMR": 0.35},
    },
    "rs4994": {
        "gene": "ADRB3", "category": "Fitness",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal beta-3 adrenergic receptor - standard fat mobilization", "magnitude": 0},
            "CT": {"status": "reduced", "desc": "ADRB3 Trp64Arg heterozygous - reduced fat mobilization, may resist weight loss", "magnitude": 2},
            "TC": {"status": "reduced", "desc": "ADRB3 Trp64Arg heterozygous - reduced fat mobilization, may resist weight loss", "magnitude": 2},
            "TT": {"status": "significantly_reduced", "desc": "ADRB3 Arg/Arg - lower metabolic rate, weight loss more difficult", "magnitude": 3},
        }
    },
    "rs1042713": {
        "gene": "ADRB2", "category": "Fitness",
        "variants": {
            "GG": {"status": "gly16", "desc": "ADRB2 Gly16 - enhanced lipolysis response to exercise", "magnitude": 1},
            "GA": {"status": "heterozygous", "desc": "ADRB2 Gly16Arg heterozygous - intermediate response", "magnitude": 0},
            "AG": {"status": "heterozygous", "desc": "ADRB2 Gly16Arg heterozygous - intermediate response", "magnitude": 0},
            "AA": {"status": "arg16", "desc": "ADRB2 Arg16 - reduced exercise-induced lipolysis", "magnitude": 1},
        }
    },
    "rs8192678": {
        "gene": "PPARGC1A", "category": "Fitness",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal PGC-1alpha - standard mitochondrial biogenesis", "magnitude": 0},
            "CT": {"status": "enhanced", "desc": "PPARGC1A Gly482Ser heterozygous - may have enhanced endurance adaptation", "magnitude": 1},
            "TC": {"status": "enhanced", "desc": "PPARGC1A Gly482Ser heterozygous - may have enhanced endurance adaptation", "magnitude": 1},
            "TT": {"status": "altered", "desc": "PPARGC1A Ser/Ser - altered mitochondrial response, may need more training volume", "magnitude": 2},
        }
    },
    "rs4253778": {
        "gene": "PPARA", "category": "Fitness",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal PPAR-alpha - standard fat oxidation", "magnitude": 0},
            "GC": {"status": "enhanced", "desc": "PPARA intron 7 C allele - enhanced fat oxidation capacity", "magnitude": 1},
            "CG": {"status": "enhanced", "desc": "PPARA intron 7 C allele - enhanced fat oxidation capacity", "magnitude": 1},
            "CC": {"status": "highly_enhanced", "desc": "PPARA C/C - superior fat oxidation, endurance advantage", "magnitude": 2},
        }
    },
    "rs1799752": {
        "gene": "ACE", "category": "Fitness",
        "variants": {
            "DD": {"status": "power", "desc": "ACE D/D - higher ACE activity, power/strength advantage", "magnitude": 2},
            "DI": {"status": "mixed", "desc": "ACE D/I - balanced ACE activity", "magnitude": 1},
            "ID": {"status": "mixed", "desc": "ACE I/D - balanced ACE activity", "magnitude": 1},
            "II": {"status": "endurance", "desc": "ACE I/I - lower ACE activity, endurance advantage, better altitude adaptation", "magnitude": 2},
            "GG": {"status": "power", "desc": "ACE D/D equivalent - power advantage", "magnitude": 2},
            "GT": {"status": "mixed", "desc": "ACE heterozygous", "magnitude": 1},
            "TG": {"status": "mixed", "desc": "ACE heterozygous", "magnitude": 1},
            "TT": {"status": "endurance", "desc": "ACE I/I equivalent - endurance advantage", "magnitude": 2},
        }
    },
    "rs7181866": {
        "gene": "COL5A1", "category": "Fitness",
        "variants": {
            "CC": {"status": "flexible", "desc": "COL5A1 C/C - more flexible tendons, lower injury risk", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "COL5A1 heterozygous - intermediate tendon properties", "magnitude": 0},
            "TC": {"status": "intermediate", "desc": "COL5A1 heterozygous - intermediate tendon properties", "magnitude": 0},
            "TT": {"status": "stiff", "desc": "COL5A1 T/T - stiffer tendons, higher injury risk, more warm-up needed", "magnitude": 2},
        }
    },
    "rs1800012": {
        "gene": "COL1A1", "category": "Fitness",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal collagen type I - standard bone/tendon strength", "magnitude": 0},
            "GT": {"status": "reduced", "desc": "COL1A1 Sp1 heterozygous - slightly reduced collagen, watch for overuse injuries", "magnitude": 1},
            "TG": {"status": "reduced", "desc": "COL1A1 Sp1 heterozygous - slightly reduced collagen, watch for overuse injuries", "magnitude": 1},
            "TT": {"status": "significantly_reduced", "desc": "COL1A1 Sp1 T/T - reduced collagen, higher fracture/injury risk", "magnitude": 2},
        }
    },

    # =========================================================================
    # SECTION 7: NUTRITION & METABOLISM
    # =========================================================================

    "rs9939609": {
        "gene": "FTO", "category": "Nutrition",
        "variants": {
            "TT": {"status": "normal", "desc": "Normal FTO - standard obesity risk", "magnitude": 0},
            "TA": {"status": "increased", "desc": "FTO risk allele heterozygous - 1.3x obesity risk, may benefit from high protein", "magnitude": 1},
            "AT": {"status": "increased", "desc": "FTO risk allele heterozygous - 1.3x obesity risk, may benefit from high protein", "magnitude": 1},
            "AA": {"status": "elevated", "desc": "FTO A/A - 1.7x obesity risk, responds well to exercise and high-protein diet", "magnitude": 2},
        }
    },
    "rs1801282": {
        "gene": "PPARG", "category": "Nutrition",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal PPAR-gamma - standard insulin sensitivity", "magnitude": 0},
            "CG": {"status": "protective", "desc": "PPARG Pro12Ala heterozygous - improved insulin sensitivity, lower T2D risk", "magnitude": 1},
            "GC": {"status": "protective", "desc": "PPARG Pro12Ala heterozygous - improved insulin sensitivity, lower T2D risk", "magnitude": 1},
            "GG": {"status": "highly_protective", "desc": "PPARG Ala/Ala - enhanced insulin sensitivity", "magnitude": 2},
        }
    },
    "rs7903146": {
        "gene": "TCF7L2", "category": "Nutrition",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal TCF7L2 - standard diabetes risk", "magnitude": 0},
            "CT": {"status": "increased", "desc": "TCF7L2 risk allele heterozygous - 1.4x T2D risk, carb restriction helpful", "magnitude": 2},
            "TC": {"status": "increased", "desc": "TCF7L2 risk allele heterozygous - 1.4x T2D risk, carb restriction helpful", "magnitude": 2},
            "TT": {"status": "elevated", "desc": "TCF7L2 T/T - 2x T2D risk, low-glycemic diet strongly recommended", "magnitude": 3},
        },
        "freq": {"EUR": 0.30, "AFR": 0.28, "EAS": 0.03, "SAS": 0.30, "AMR": 0.25},
    },
    "rs5082": {
        "gene": "APOA2", "category": "Nutrition",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal APOA2 - standard saturated fat response", "magnitude": 0},
            "GA": {"status": "intermediate", "desc": "APOA2 -265T>C heterozygous - intermediate sat fat sensitivity", "magnitude": 1},
            "AG": {"status": "intermediate", "desc": "APOA2 -265T>C heterozygous - intermediate sat fat sensitivity", "magnitude": 1},
            "AA": {"status": "sensitive", "desc": "APOA2 C/C - saturated fat intake strongly linked to obesity, limit sat fat", "magnitude": 2},
        }
    },
    "rs174547": {
        "gene": "FADS1", "category": "Nutrition",
        "variants": {
            "TT": {"status": "high_conversion", "desc": "FADS1 T/T - efficient omega-3/6 conversion, plant sources adequate", "magnitude": 1},
            "TC": {"status": "intermediate", "desc": "FADS1 heterozygous - moderate conversion efficiency", "magnitude": 0},
            "CT": {"status": "intermediate", "desc": "FADS1 heterozygous - moderate conversion efficiency", "magnitude": 0},
            "CC": {"status": "low_conversion", "desc": "FADS1 C/C - poor ALA to EPA/DHA conversion, direct fish oil/algae preferred", "magnitude": 2},
        }
    },
    "rs4988235": {
        "gene": "MCM6/LCT", "category": "Nutrition",
        "variants": {
            "AA": {"status": "lactose_intolerant", "desc": "Lactase non-persistence - lactose intolerance in adulthood", "magnitude": 2},
            "AG": {"status": "tolerant", "desc": "Lactase persistence heterozygous - likely lactose tolerant", "magnitude": 0},
            "GA": {"status": "tolerant", "desc": "Lactase persistence heterozygous - likely lactose tolerant", "magnitude": 0},
            "GG": {"status": "tolerant", "desc": "Lactase persistence - maintains lactase production, lactose tolerant", "magnitude": 0},
        },
        "freq": {"EUR": 0.75, "AFR": 0.10, "EAS": 0.01, "SAS": 0.25, "AMR": 0.40},
    },
    "rs2282679": {
        "gene": "GC", "category": "Nutrition",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal vitamin D binding protein - adequate vitamin D transport", "magnitude": 0},
            "GT": {"status": "reduced", "desc": "GC variant heterozygous - lower vitamin D levels common", "magnitude": 1},
            "TG": {"status": "reduced", "desc": "GC variant heterozygous - lower vitamin D levels common", "magnitude": 1},
            "TT": {"status": "low", "desc": "GC T/T - genetically low vitamin D, supplementation often needed especially in northern latitudes", "magnitude": 2},
        }
    },
    "rs12934922": {
        "gene": "BCMO1", "category": "Nutrition",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal beta-carotene conversion to vitamin A", "magnitude": 0},
            "AT": {"status": "reduced", "desc": "BCMO1 A379V heterozygous - ~30% reduced conversion, consider preformed vitamin A", "magnitude": 1},
            "TA": {"status": "reduced", "desc": "BCMO1 A379V heterozygous - ~30% reduced conversion, consider preformed vitamin A", "magnitude": 1},
            "TT": {"status": "significantly_reduced", "desc": "BCMO1 T/T - ~70% reduced beta-carotene conversion, need preformed vitamin A sources", "magnitude": 2},
        }
    },
    "rs2228570": {
        "gene": "VDR", "category": "Nutrition",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal VDR FokI - standard vitamin D receptor function", "magnitude": 0},
            "CT": {"status": "reduced", "desc": "VDR FokI heterozygous - slightly reduced vitamin D receptor activity", "magnitude": 1},
            "TC": {"status": "reduced", "desc": "VDR FokI heterozygous - slightly reduced vitamin D receptor activity", "magnitude": 1},
            "TT": {"status": "significantly_reduced", "desc": "VDR FokI T/T - reduced vitamin D receptor function, may need higher vitamin D levels", "magnitude": 2},
        },
        "note": "FokI variant affects VDR protein length. T allele produces a longer, less active receptor."
    },
    "rs1544410": {
        "gene": "VDR", "category": "Nutrition",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal VDR BsmI - standard bone mineral density", "magnitude": 0},
            "CT": {"status": "reduced", "desc": "VDR BsmI heterozygous - associated with lower bone mineral density", "magnitude": 1},
            "TC": {"status": "reduced", "desc": "VDR BsmI heterozygous - associated with lower bone mineral density", "magnitude": 1},
            "TT": {"status": "low_bmd", "desc": "VDR BsmI T/T - associated with reduced bone mineral density, ensure adequate vitamin D + calcium", "magnitude": 2},
        },
        "note": "BsmI polymorphism affects VDR mRNA stability and bone metabolism."
    },
    "rs602662": {
        "gene": "FUT2", "category": "Nutrition",
        "variants": {
            "GG": {"status": "secretor", "desc": "Secretor status - normal B12 absorption, standard gut microbiome", "magnitude": 0},
            "GA": {"status": "secretor", "desc": "Secretor heterozygous - normal B12 absorption", "magnitude": 0},
            "AG": {"status": "secretor", "desc": "Secretor heterozygous - normal B12 absorption", "magnitude": 0},
            "AA": {"status": "non_secretor", "desc": "Non-secretor - may have lower B12 levels, different gut microbiome, consider B12 monitoring", "magnitude": 2},
        }
    },

    # =========================================================================
    # SECTION 8: CARDIOVASCULAR
    # =========================================================================

    "rs429358": {
        "gene": "APOE", "category": "Cardiovascular",
        "note": "Combine with rs7412 for APOE type",
        "variants": {
            "TT": {"status": "e2_or_e3", "desc": "APOE not e4 at this position", "magnitude": 0},
            "TC": {"status": "e4_carrier", "desc": "APOE e4 carrier (one copy) - increased CVD and Alzheimer's risk", "magnitude": 3},
            "CT": {"status": "e4_carrier", "desc": "APOE e4 carrier (one copy) - increased CVD and Alzheimer's risk", "magnitude": 3},
            "CC": {"status": "e4_homozygous", "desc": "APOE e4/e4 - significantly elevated Alzheimer's risk (10-15x)", "magnitude": 5},
        },
        "freq": {"EUR": 0.15, "AFR": 0.27, "EAS": 0.09, "SAS": 0.12, "AMR": 0.11},
    },
    "rs7412": {
        "gene": "APOE", "category": "Cardiovascular",
        "note": "Combine with rs429358 for APOE type",
        "variants": {
            "CC": {"status": "e3_or_e4", "desc": "Not APOE e2 at this position", "magnitude": 0},
            "CT": {"status": "e2_carrier", "desc": "APOE e2 carrier - may be protective against Alzheimer's", "magnitude": 1},
            "TC": {"status": "e2_carrier", "desc": "APOE e2 carrier - may be protective against Alzheimer's", "magnitude": 1},
            "TT": {"status": "e2_homozygous", "desc": "APOE e2/e2 - protective vs Alzheimer's but watch Type III hyperlipoproteinemia", "magnitude": 2},
        }
    },
    "rs6025": {
        "gene": "F5", "category": "Cardiovascular",
        "variants": {
            "CC": {"status": "normal", "desc": "No Factor V Leiden - normal clotting", "magnitude": 0},
            "CT": {"status": "carrier", "desc": "Factor V Leiden heterozygous - 5-10x DVT risk, avoid estrogen contraceptives", "magnitude": 4},
            "TC": {"status": "carrier", "desc": "Factor V Leiden heterozygous - 5-10x DVT risk, avoid estrogen contraceptives", "magnitude": 4},
            "TT": {"status": "homozygous", "desc": "Factor V Leiden homozygous - 50-100x DVT risk", "magnitude": 5},
        },
        "freq": {"EUR": 0.05, "AFR": 0.01, "EAS": 0.001, "SAS": 0.01, "AMR": 0.02},
    },
    "rs1799963": {
        "gene": "F2", "category": "Cardiovascular",
        "variants": {
            "GG": {"status": "normal", "desc": "No prothrombin mutation - normal clotting", "magnitude": 0},
            "GA": {"status": "carrier", "desc": "Prothrombin G20210A heterozygous - 3x DVT risk", "magnitude": 3},
            "AG": {"status": "carrier", "desc": "Prothrombin G20210A heterozygous - 3x DVT risk", "magnitude": 3},
            "AA": {"status": "homozygous", "desc": "Prothrombin G20210A homozygous - significantly elevated clot risk", "magnitude": 4},
        }
    },
    "rs5186": {
        "gene": "AGTR1", "category": "Cardiovascular",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal angiotensin II receptor - standard blood pressure response", "magnitude": 0},
            "AC": {"status": "increased", "desc": "AGTR1 A1166C heterozygous - increased hypertension risk", "magnitude": 2},
            "CA": {"status": "increased", "desc": "AGTR1 A1166C heterozygous - increased hypertension risk", "magnitude": 2},
            "CC": {"status": "elevated", "desc": "AGTR1 C/C - elevated hypertension risk, salt-sensitive, responds well to ARBs", "magnitude": 3},
        }
    },
    "rs699": {
        "gene": "AGT", "category": "Cardiovascular",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal angiotensinogen - standard blood pressure", "magnitude": 0},
            "AG": {"status": "increased", "desc": "AGT M235T heterozygous - ~20% higher AGT, slightly elevated BP risk", "magnitude": 1},
            "GA": {"status": "increased", "desc": "AGT M235T heterozygous - ~20% higher AGT, slightly elevated BP risk", "magnitude": 1},
            "GG": {"status": "elevated", "desc": "AGT T/T - ~40% higher AGT levels, elevated hypertension risk, salt restriction helpful", "magnitude": 2},
        }
    },
    "rs4343": {
        "gene": "ACE", "category": "Cardiovascular",
        "variants": {
            "AA": {"status": "low", "desc": "Lower ACE activity - better endurance, lower BP tendency", "magnitude": 1},
            "AG": {"status": "intermediate", "desc": "Intermediate ACE activity", "magnitude": 0},
            "GA": {"status": "intermediate", "desc": "Intermediate ACE activity", "magnitude": 0},
            "GG": {"status": "high", "desc": "Higher ACE activity - power advantage but higher BP risk, responds well to ACE inhibitors", "magnitude": 2},
        }
    },
    "rs5443": {
        "gene": "GNB3", "category": "Cardiovascular",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal G-protein signaling", "magnitude": 0},
            "CT": {"status": "increased", "desc": "GNB3 C825T heterozygous - increased hypertension and obesity risk", "magnitude": 1},
            "TC": {"status": "increased", "desc": "GNB3 C825T heterozygous - increased hypertension and obesity risk", "magnitude": 1},
            "TT": {"status": "elevated", "desc": "GNB3 T/T - elevated hypertension risk, responds well to diuretics", "magnitude": 2},
        }
    },
    "rs1801253": {
        "gene": "ADRB1", "category": "Cardiovascular",
        "variants": {
            "CC": {"status": "arg389", "desc": "ADRB1 Arg389 - enhanced beta-blocker response", "magnitude": 1},
            "CG": {"status": "heterozygous", "desc": "ADRB1 Arg389Gly heterozygous - intermediate beta-blocker response", "magnitude": 0},
            "GC": {"status": "heterozygous", "desc": "ADRB1 Arg389Gly heterozygous - intermediate beta-blocker response", "magnitude": 0},
            "GG": {"status": "gly389", "desc": "ADRB1 Gly389 - reduced beta-blocker efficacy, may need dose adjustment", "magnitude": 2},
        }
    },
    "rs1800629": {
        "gene": "TNF", "category": "Inflammation",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal TNF-alpha levels", "magnitude": 0},
            "GA": {"status": "increased", "desc": "TNF-308 G>A heterozygous - higher TNF-alpha, increased inflammation", "magnitude": 2},
            "AG": {"status": "increased", "desc": "TNF-308 G>A heterozygous - higher TNF-alpha, increased inflammation", "magnitude": 2},
            "AA": {"status": "elevated", "desc": "TNF-308 A/A - significantly elevated TNF-alpha, chronic inflammation risk", "magnitude": 3},
        }
    },
    "rs1800795": {
        "gene": "IL6", "category": "Inflammation",
        "variants": {
            "GG": {"status": "high", "desc": "IL-6 -174 G/G - higher baseline IL-6, more inflammatory response", "magnitude": 2},
            "GC": {"status": "intermediate", "desc": "IL-6 -174 heterozygous - intermediate IL-6 levels", "magnitude": 1},
            "CG": {"status": "intermediate", "desc": "IL-6 -174 heterozygous - intermediate IL-6 levels", "magnitude": 1},
            "CC": {"status": "low", "desc": "IL-6 -174 C/C - lower baseline inflammation", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION 9: IRON & MINERALS
    # =========================================================================

    "rs1800562": {
        "gene": "HFE", "category": "Iron Metabolism",
        "variants": {
            "GG": {"status": "normal", "desc": "No C282Y HFE mutation - normal iron regulation", "magnitude": 0},
            "GA": {"status": "carrier", "desc": "HFE C282Y carrier - monitor iron levels periodically", "magnitude": 2},
            "AG": {"status": "carrier", "desc": "HFE C282Y carrier - monitor iron levels periodically", "magnitude": 2},
            "AA": {"status": "at_risk", "desc": "HFE C282Y homozygous - hereditary hemochromatosis risk, regular iron monitoring essential", "magnitude": 4},
        },
        "freq": {"EUR": 0.08, "AFR": 0.001, "EAS": 0.001, "SAS": 0.001, "AMR": 0.03},
    },
    "rs1799945": {
        "gene": "HFE", "category": "Iron Metabolism",
        "variants": {
            "CC": {"status": "normal", "desc": "No H63D HFE mutation", "magnitude": 0},
            "CG": {"status": "carrier", "desc": "HFE H63D carrier - mild iron accumulation possible", "magnitude": 1},
            "GC": {"status": "carrier", "desc": "HFE H63D carrier - mild iron accumulation possible", "magnitude": 1},
            "GG": {"status": "homozygous", "desc": "HFE H63D homozygous - elevated iron possible, especially with C282Y", "magnitude": 2},
        }
    },

    # =========================================================================
    # SECTION 10: AUTOIMMUNE & DISEASE RISK
    # =========================================================================

    "rs2187668": {
        "gene": "HLA-DQA1", "category": "Autoimmune",
        "variants": {
            "CC": {"status": "low_risk", "desc": "Lower celiac disease risk", "magnitude": 0},
            "CT": {"status": "increased_risk", "desc": "HLA-DQ2.5 carrier - increased celiac disease risk", "magnitude": 2},
            "TC": {"status": "increased_risk", "desc": "HLA-DQ2.5 carrier - increased celiac disease risk", "magnitude": 2},
            "TT": {"status": "high_risk", "desc": "HLA-DQ2.5 homozygous - highest celiac disease risk", "magnitude": 3},
        }
    },
    "rs7574865": {
        "gene": "STAT4", "category": "Autoimmune",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal autoimmune risk at this locus", "magnitude": 0},
            "GT": {"status": "increased", "desc": "STAT4 risk allele - increased RA, lupus risk", "magnitude": 1},
            "TG": {"status": "increased", "desc": "STAT4 risk allele - increased RA, lupus risk", "magnitude": 1},
            "TT": {"status": "elevated", "desc": "STAT4 T/T - elevated autoimmune disease risk", "magnitude": 2},
        }
    },
    "rs6457620": {
        "gene": "HLA-DRB1", "category": "Autoimmune",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal HLA-DRB1 risk - lower RA/T1D susceptibility", "magnitude": 0},
            "AG": {"status": "increased", "desc": "HLA-DRB1 tag heterozygous - increased rheumatoid arthritis and type 1 diabetes risk", "magnitude": 2},
            "GA": {"status": "increased", "desc": "HLA-DRB1 tag heterozygous - increased rheumatoid arthritis and type 1 diabetes risk", "magnitude": 2},
            "GG": {"status": "elevated", "desc": "HLA-DRB1 tag homozygous - significantly elevated RA and T1D risk", "magnitude": 3},
        },
        "note": "Tag SNP for HLA-DRB1 shared epitope alleles associated with RA and T1D."
    },
    "rs3134792": {
        "gene": "HLA-B27", "category": "Autoimmune",
        "variants": {
            "CC": {"status": "normal", "desc": "HLA-B27 proxy negative - lower ankylosing spondylitis risk", "magnitude": 0},
            "CT": {"status": "carrier", "desc": "HLA-B27 proxy heterozygous - increased ankylosing spondylitis risk (~5-6% develop AS)", "magnitude": 2},
            "TC": {"status": "carrier", "desc": "HLA-B27 proxy heterozygous - increased ankylosing spondylitis risk (~5-6% develop AS)", "magnitude": 2},
            "TT": {"status": "positive", "desc": "HLA-B27 proxy positive - elevated risk for ankylosing spondylitis and reactive arthritis", "magnitude": 3},
        },
        "note": "Proxy for HLA-B27. ~8% of European population is HLA-B27+, but only 5-6% of carriers develop AS."
    },
    "rs2476601": {
        "gene": "PTPN22", "category": "Autoimmune",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal autoimmune risk", "magnitude": 0},
            "GA": {"status": "increased", "desc": "PTPN22 R620W heterozygous - increased T1D, RA, thyroid autoimmunity risk", "magnitude": 2},
            "AG": {"status": "increased", "desc": "PTPN22 R620W heterozygous - increased T1D, RA, thyroid autoimmunity risk", "magnitude": 2},
            "AA": {"status": "elevated", "desc": "PTPN22 W/W - significantly elevated autoimmune risk", "magnitude": 3},
        }
    },

    # =========================================================================
    # SECTION 11: SKIN & AGING
    # =========================================================================

    "rs1805007": {
        "gene": "MC1R", "category": "Skin",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal MC1R - standard sun sensitivity", "magnitude": 0},
            "CT": {"status": "carrier", "desc": "MC1R R151C carrier - increased sun sensitivity, freckling, skin cancer risk", "magnitude": 2},
            "TC": {"status": "carrier", "desc": "MC1R R151C carrier - increased sun sensitivity, freckling, skin cancer risk", "magnitude": 2},
            "TT": {"status": "high_risk", "desc": "MC1R R151C homozygous - red hair phenotype, very high sun sensitivity", "magnitude": 3},
        }
    },
    "rs1805008": {
        "gene": "MC1R", "category": "Skin",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal MC1R R160W", "magnitude": 0},
            "CT": {"status": "carrier", "desc": "MC1R R160W carrier - increased sun sensitivity", "magnitude": 2},
            "TC": {"status": "carrier", "desc": "MC1R R160W carrier - increased sun sensitivity", "magnitude": 2},
            "TT": {"status": "high_risk", "desc": "MC1R R160W homozygous - high sun/skin cancer risk", "magnitude": 3},
        }
    },
    "rs12203592": {
        "gene": "IRF4", "category": "Skin",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal pigmentation regulation", "magnitude": 0},
            "CT": {"status": "lighter", "desc": "IRF4 variant - tendency toward lighter skin, increased sun sensitivity", "magnitude": 1},
            "TC": {"status": "lighter", "desc": "IRF4 variant - tendency toward lighter skin, increased sun sensitivity", "magnitude": 1},
            "TT": {"status": "very_light", "desc": "IRF4 T/T - very light skin, high sun sensitivity, extra sun protection needed", "magnitude": 2},
        }
    },
    "rs2228479": {
        "gene": "MC1R", "category": "Skin",
        "variants": {
            "AA": {"status": "normal", "desc": "Normal MC1R Val92Met", "magnitude": 0},
            "AG": {"status": "variant", "desc": "MC1R V92M heterozygous - slightly increased skin aging risk", "magnitude": 1},
            "GA": {"status": "variant", "desc": "MC1R V92M heterozygous - slightly increased skin aging risk", "magnitude": 1},
            "GG": {"status": "accelerated", "desc": "MC1R V92M homozygous - may show accelerated skin aging", "magnitude": 2},
        }
    },

    # =========================================================================
    # SECTION 12: LONGEVITY & AGING
    # =========================================================================

    "rs2802292": {
        "gene": "FOXO3", "category": "Longevity",
        "variants": {
            "TT": {"status": "normal", "desc": "Normal FOXO3 - standard longevity", "magnitude": 0},
            "TG": {"status": "favorable", "desc": "FOXO3 longevity variant heterozygous - associated with increased lifespan", "magnitude": 1},
            "GT": {"status": "favorable", "desc": "FOXO3 longevity variant heterozygous - associated with increased lifespan", "magnitude": 1},
            "GG": {"status": "highly_favorable", "desc": "FOXO3 G/G - strongly associated with longevity in multiple populations", "magnitude": 2},
        }
    },
    "rs1042522": {
        "gene": "TP53", "category": "Longevity",
        "variants": {
            "GG": {"status": "pro72", "desc": "TP53 Pro72 - more efficient apoptosis, may be protective against cancer", "magnitude": 1},
            "GC": {"status": "heterozygous", "desc": "TP53 Pro72Arg heterozygous - balanced", "magnitude": 0},
            "CG": {"status": "heterozygous", "desc": "TP53 Pro72Arg heterozygous - balanced", "magnitude": 0},
            "CC": {"status": "arg72", "desc": "TP53 Arg72 - less efficient apoptosis", "magnitude": 1},
        }
    },
    "rs2542052": {
        "gene": "CETP", "category": "Longevity",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal CETP activity - standard HDL metabolism", "magnitude": 0},
            "CA": {"status": "favorable", "desc": "CETP I405V heterozygous - higher HDL, longevity association", "magnitude": 1},
            "AC": {"status": "favorable", "desc": "CETP I405V heterozygous - higher HDL, longevity association", "magnitude": 1},
            "AA": {"status": "highly_favorable", "desc": "CETP V/V - significantly higher HDL, associated with longevity", "magnitude": 2},
        }
    },

    # =========================================================================
    # SECTION 13: RESPIRATORY & LUNG
    # =========================================================================

    "rs28929474": {
        "gene": "SERPINA1", "category": "Respiratory",
        "variants": {
            "CC": {"status": "normal", "desc": "Normal alpha-1 antitrypsin", "magnitude": 0},
            "CT": {"status": "carrier", "desc": "Alpha-1 antitrypsin Pi*Z carrier - avoid smoking, monitor lung function", "magnitude": 3},
            "TC": {"status": "carrier", "desc": "Alpha-1 antitrypsin Pi*Z carrier - avoid smoking, monitor lung function", "magnitude": 3},
            "TT": {"status": "deficient", "desc": "Alpha-1 antitrypsin deficiency (Pi*ZZ) - high COPD/liver disease risk", "magnitude": 5},
        },
        "freq": {"EUR": 0.02, "AFR": 0.001, "EAS": 0.001, "SAS": 0.001, "AMR": 0.01},
    },

    # =========================================================================
    # SECTION 14: ALCOHOL METABOLISM
    # =========================================================================

    "rs671": {
        "gene": "ALDH2", "category": "Alcohol",
        "variants": {
            "GG": {"status": "normal", "desc": "Normal ALDH2 - efficient alcohol metabolism", "magnitude": 0},
            "GA": {"status": "reduced", "desc": "ALDH2*2 heterozygous - alcohol flush, increased esophageal cancer risk with drinking", "magnitude": 3},
            "AG": {"status": "reduced", "desc": "ALDH2*2 heterozygous - alcohol flush, increased esophageal cancer risk with drinking", "magnitude": 3},
            "AA": {"status": "deficient", "desc": "ALDH2*2 homozygous - severe alcohol intolerance, avoid alcohol", "magnitude": 4},
        },
        "freq": {"EUR": 0.001, "AFR": 0.001, "EAS": 0.30, "SAS": 0.03, "AMR": 0.01},
    },
    "rs1229984": {
        "gene": "ADH1B", "category": "Alcohol",
        "variants": {
            "CC": {"status": "slow", "desc": "Slower alcohol metabolism - alcohol effects last longer", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "Intermediate alcohol metabolism", "magnitude": 0},
            "TC": {"status": "intermediate", "desc": "Intermediate alcohol metabolism", "magnitude": 0},
            "TT": {"status": "fast", "desc": "Fast alcohol metabolism - protective against alcoholism", "magnitude": 1},
        }
    },

    # =========================================================================
    # SECTION 15: BLOOD TYPE
    # =========================================================================

    "rs505922": {
        "gene": "ABO", "category": "Blood Type",
        "variants": {
            "CC": {"status": "A_or_B", "desc": "ABO proxy: likely blood type A or B (C allele)", "magnitude": 1},
            "CT": {"status": "A_or_B_carrier_O", "desc": "ABO proxy: A/B with one O allele", "magnitude": 1},
            "TC": {"status": "A_or_B_carrier_O", "desc": "ABO proxy: A/B with one O allele", "magnitude": 1},
            "TT": {"status": "likely_O", "desc": "ABO proxy: likely blood type O (TT strongly associated with O type)", "magnitude": 1},
        },
        "note": "rs505922 is a reliable proxy for ABO blood type. T allele strongly associated with O type."
    },
    "rs8176746": {
        "gene": "ABO", "category": "Blood Type",
        "variants": {
            "CC": {"status": "not_B", "desc": "ABO B allele absent - blood type A or O", "magnitude": 0},
            "CT": {"status": "B_carrier", "desc": "ABO B allele heterozygous - possible blood type B or AB", "magnitude": 1},
            "TC": {"status": "B_carrier", "desc": "ABO B allele heterozygous - possible blood type B or AB", "magnitude": 1},
            "TT": {"status": "B_likely", "desc": "ABO B allele homozygous - likely blood type B", "magnitude": 1},
        },
        "note": "T allele at rs8176746 defines the B antigen."
    },
    "rs590787": {
        "gene": "RHD", "category": "Blood Type",
        "variants": {
            "CC": {"status": "Rh_positive", "desc": "Rh factor positive (most common)", "magnitude": 0},
            "CT": {"status": "Rh_positive", "desc": "Rh factor positive (carrier of Rh- allele)", "magnitude": 0},
            "TC": {"status": "Rh_positive", "desc": "Rh factor positive (carrier of Rh- allele)", "magnitude": 0},
            "TT": {"status": "Rh_negative", "desc": "Rh factor likely negative - relevant for pregnancy and transfusions", "magnitude": 2},
        },
        "note": "Proxy SNP for RhD status. Rh-negative is ~15% in Europeans, <1% in East Asians."
    },

    # =========================================================================
    # SECTION 16: BITTER TASTE
    # =========================================================================

    "rs713598": {
        "gene": "TAS2R38", "category": "Taste",
        "variants": {
            "GG": {"status": "taster", "desc": "TAS2R38 A49P: strong bitter taster (PAV haplotype) - perceives PTC/PROP as very bitter", "magnitude": 1},
            "GC": {"status": "intermediate", "desc": "TAS2R38 A49P heterozygous: moderate bitter taste perception", "magnitude": 0},
            "CG": {"status": "intermediate", "desc": "TAS2R38 A49P heterozygous: moderate bitter taste perception", "magnitude": 0},
            "CC": {"status": "non_taster", "desc": "TAS2R38 A49P: non-taster (AVI haplotype) - reduced bitter taste, may eat more cruciferous vegetables", "magnitude": 1},
        },
        "note": "Most important variant for bitter taste. Affects preference for cruciferous vegetables, coffee, beer."
    },
    "rs1726866": {
        "gene": "TAS2R38", "category": "Taste",
        "variants": {
            "GG": {"status": "taster", "desc": "TAS2R38 V262A: taster allele", "magnitude": 0},
            "GA": {"status": "intermediate", "desc": "TAS2R38 V262A heterozygous", "magnitude": 0},
            "AG": {"status": "intermediate", "desc": "TAS2R38 V262A heterozygous", "magnitude": 0},
            "AA": {"status": "non_taster", "desc": "TAS2R38 V262A: non-taster allele", "magnitude": 0},
        }
    },
    "rs10246939": {
        "gene": "TAS2R38", "category": "Taste",
        "variants": {
            "CC": {"status": "taster", "desc": "TAS2R38 I296V: taster allele (I296)", "magnitude": 0},
            "CT": {"status": "intermediate", "desc": "TAS2R38 I296V heterozygous", "magnitude": 0},
            "TC": {"status": "intermediate", "desc": "TAS2R38 I296V heterozygous", "magnitude": 0},
            "TT": {"status": "non_taster", "desc": "TAS2R38 I296V: non-taster allele (V296)", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: LONGEVITY & AGING
    # =========================================================================

    "rs2536": {
        "gene": "MTOR", "category": "Longevity",
        "variants": {
            "CC": {"status": "protective", "desc": "mTOR variant associated with longevity signaling", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "mTOR heterozygous — partial longevity effect", "magnitude": 0},
            "TC": {"status": "intermediate", "desc": "mTOR heterozygous — partial longevity effect", "magnitude": 0},
            "TT": {"status": "reference", "desc": "mTOR reference — standard aging pathway activity", "magnitude": 0},
        }
    },
    "rs7069102": {
        "gene": "SIRT1", "category": "Longevity",
        "variants": {
            "GG": {"status": "protective", "desc": "SIRT1 GG — enhanced NAD+ deacetylase, caloric restriction pathway", "magnitude": 2},
            "GC": {"status": "intermediate", "desc": "SIRT1 heterozygous — partial sirtuin benefit", "magnitude": 1},
            "CG": {"status": "intermediate", "desc": "SIRT1 heterozygous — partial sirtuin benefit", "magnitude": 1},
            "CC": {"status": "reference", "desc": "SIRT1 reference — standard sirtuin activity", "magnitude": 0},
        }
    },
    "rs9536314": {
        "gene": "KLOTHO", "category": "Longevity",
        "variants": {
            "TT": {"status": "protective", "desc": "Klotho TT — enhanced anti-aging hormone, kidney & brain protection", "magnitude": 2},
            "TG": {"status": "intermediate", "desc": "Klotho heterozygous — partial anti-aging benefit", "magnitude": 1},
            "GT": {"status": "intermediate", "desc": "Klotho heterozygous — partial anti-aging benefit", "magnitude": 1},
            "GG": {"status": "reference", "desc": "Klotho reference — standard aging hormone levels", "magnitude": 0},
        }
    },
    "rs2229765": {
        "gene": "IGF1R", "category": "Longevity",
        "variants": {
            "AA": {"status": "protective", "desc": "IGF1R AA — reduced IGF-1 signaling associated with longevity", "magnitude": 2},
            "AG": {"status": "intermediate", "desc": "IGF1R heterozygous — partial growth pathway reduction", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "IGF1R heterozygous — partial growth pathway reduction", "magnitude": 1},
            "GG": {"status": "reference", "desc": "IGF1R reference — standard growth hormone signaling", "magnitude": 0},
        }
    },
    "rs10936599": {
        "gene": "TERC", "category": "Longevity",
        "variants": {
            "CC": {"status": "protective", "desc": "TERC CC — longer telomere maintenance, cellular longevity", "magnitude": 2},
            "CT": {"status": "intermediate", "desc": "TERC heterozygous — moderate telomere benefit", "magnitude": 1},
            "TC": {"status": "intermediate", "desc": "TERC heterozygous — moderate telomere benefit", "magnitude": 1},
            "TT": {"status": "reduced", "desc": "TERC TT — shorter telomeres, accelerated cellular aging", "magnitude": 2},
        }
    },
    "rs2736100": {
        "gene": "TERT", "category": "Longevity",
        "variants": {
            "CC": {"status": "protective", "desc": "TERT CC — enhanced telomerase activity, telomere maintenance", "magnitude": 2},
            "CA": {"status": "intermediate", "desc": "TERT heterozygous — moderate telomerase activity", "magnitude": 1},
            "AC": {"status": "intermediate", "desc": "TERT heterozygous — moderate telomerase activity", "magnitude": 1},
            "AA": {"status": "reduced", "desc": "TERT AA — lower telomerase activity", "magnitude": 1},
        }
    },

    # =========================================================================
    # SECTION: SLEEP & CIRCADIAN
    # =========================================================================

    "rs228697": {
        "gene": "PER3", "category": "Sleep/Circadian",
        "variants": {
            "CC": {"status": "evening", "desc": "PER3 CC — short-sleep tolerance, evening chronotype tendency", "magnitude": 2},
            "CG": {"status": "intermediate", "desc": "PER3 heterozygous — moderate chronotype effect", "magnitude": 1},
            "GC": {"status": "intermediate", "desc": "PER3 heterozygous — moderate chronotype effect", "magnitude": 1},
            "GG": {"status": "morning", "desc": "PER3 GG — morning chronotype, needs full sleep duration", "magnitude": 1},
        }
    },
    "rs2305160": {
        "gene": "NPAS2", "category": "Sleep/Circadian",
        "variants": {
            "GG": {"status": "evening", "desc": "NPAS2 GG — evening chronotype preference", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "NPAS2 heterozygous — moderate chronotype effect", "magnitude": 0},
            "AG": {"status": "intermediate", "desc": "NPAS2 heterozygous — moderate chronotype effect", "magnitude": 0},
            "AA": {"status": "morning", "desc": "NPAS2 AA — morning chronotype preference", "magnitude": 1},
        }
    },

    # =========================================================================
    # SECTION: MENTAL HEALTH
    # =========================================================================

    "rs6295": {
        "gene": "HTR1A", "category": "Mental Health",
        "variants": {
            "GG": {"status": "risk", "desc": "HTR1A GG — reduced serotonin autoreceptor feedback, depression/SSRI resistance risk", "magnitude": 3},
            "GC": {"status": "intermediate", "desc": "HTR1A heterozygous — moderate serotonin receptor effect", "magnitude": 1},
            "CG": {"status": "intermediate", "desc": "HTR1A heterozygous — moderate serotonin receptor effect", "magnitude": 1},
            "CC": {"status": "normal", "desc": "HTR1A CC — normal serotonin autoreceptor function", "magnitude": 0},
        }
    },
    "rs1800532": {
        "gene": "TPH1", "category": "Mental Health",
        "variants": {
            "AA": {"status": "reduced", "desc": "TPH1 AA — reduced serotonin synthesis, depression susceptibility", "magnitude": 2},
            "AC": {"status": "intermediate", "desc": "TPH1 heterozygous — moderate serotonin synthesis reduction", "magnitude": 1},
            "CA": {"status": "intermediate", "desc": "TPH1 heterozygous — moderate serotonin synthesis reduction", "magnitude": 1},
            "CC": {"status": "normal", "desc": "TPH1 CC — normal serotonin synthesis rate", "magnitude": 0},
        }
    },
    "rs110402": {
        "gene": "CRHR1", "category": "Mental Health",
        "variants": {
            "AA": {"status": "resilient", "desc": "CRHR1 AA — may be protective after trauma exposure", "magnitude": 2},
            "AG": {"status": "intermediate", "desc": "CRHR1 heterozygous — moderate stress response modulation", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "CRHR1 heterozygous — moderate stress response modulation", "magnitude": 1},
            "GG": {"status": "standard", "desc": "CRHR1 GG — standard HPA axis stress response", "magnitude": 0},
        }
    },
    "rs9296158": {
        "gene": "FKBP5", "category": "Mental Health",
        "variants": {
            "AA": {"status": "risk", "desc": "FKBP5 AA — impaired cortisol recovery, trauma sensitivity", "magnitude": 3},
            "AG": {"status": "intermediate", "desc": "FKBP5 heterozygous — moderate cortisol recovery effect", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "FKBP5 heterozygous — moderate cortisol recovery effect", "magnitude": 1},
            "GG": {"status": "normal", "desc": "FKBP5 GG — normal cortisol feedback", "magnitude": 0},
        }
    },
    "rs279858": {
        "gene": "GABRA2", "category": "Mental Health",
        "variants": {
            "GG": {"status": "risk", "desc": "GABRA2 GG — elevated alcohol dependence susceptibility", "magnitude": 2},
            "GA": {"status": "intermediate", "desc": "GABRA2 heterozygous — moderate risk", "magnitude": 1},
            "AG": {"status": "intermediate", "desc": "GABRA2 heterozygous — moderate risk", "magnitude": 1},
            "AA": {"status": "normal", "desc": "GABRA2 AA — typical GABA-A receptor function", "magnitude": 0},
        }
    },
    "rs16969968": {
        "gene": "CHRNA5", "category": "Mental Health",
        "variants": {
            "AA": {"status": "risk", "desc": "CHRNA5 AA — high nicotine dependence risk, increased cigarettes/day", "magnitude": 3},
            "AG": {"status": "intermediate", "desc": "CHRNA5 heterozygous — moderate nicotine dependence risk", "magnitude": 2},
            "GA": {"status": "intermediate", "desc": "CHRNA5 heterozygous — moderate nicotine dependence risk", "magnitude": 2},
            "GG": {"status": "normal", "desc": "CHRNA5 GG — typical nicotine receptor sensitivity", "magnitude": 0},
        }
    },
    "rs1611115": {
        "gene": "DBH", "category": "Mental Health",
        "variants": {
            "TT": {"status": "reduced", "desc": "DBH TT — reduced dopamine-to-norepinephrine conversion, ADHD trait association", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "DBH heterozygous — moderate dopamine/norepinephrine balance", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "DBH heterozygous — moderate dopamine/norepinephrine balance", "magnitude": 1},
            "CC": {"status": "normal", "desc": "DBH CC — normal dopamine beta-hydroxylase activity", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: SKIN HEALTH & AGING
    # =========================================================================

    "rs1800255": {
        "gene": "COL3A1", "category": "Skin Health",
        "variants": {
            "GG": {"status": "reduced", "desc": "COL3A1 GG — reduced type III collagen, lower skin elasticity", "magnitude": 2},
            "GA": {"status": "intermediate", "desc": "COL3A1 heterozygous — moderate collagen effect", "magnitude": 1},
            "AG": {"status": "intermediate", "desc": "COL3A1 heterozygous — moderate collagen effect", "magnitude": 1},
            "AA": {"status": "normal", "desc": "COL3A1 AA — normal type III collagen production", "magnitude": 0},
        }
    },
    "rs11204681": {
        "gene": "FLG", "category": "Skin Health",
        "variants": {
            "AA": {"status": "risk", "desc": "FLG AA — filaggrin loss, impaired skin barrier, eczema/atopic dermatitis risk", "magnitude": 3},
            "AT": {"status": "carrier", "desc": "FLG heterozygous — carrier for impaired skin barrier", "magnitude": 2},
            "TA": {"status": "carrier", "desc": "FLG heterozygous — carrier for impaired skin barrier", "magnitude": 2},
            "TT": {"status": "normal", "desc": "FLG TT — normal filaggrin, intact skin barrier", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: ATHLETIC PERFORMANCE
    # =========================================================================

    "rs6746030": {
        "gene": "SCN9A", "category": "Fitness",
        "variants": {
            "AA": {"status": "low_pain", "desc": "SCN9A AA — reduced pain sensitivity (sodium channel variant)", "magnitude": 2},
            "AG": {"status": "intermediate", "desc": "SCN9A heterozygous — moderate pain sensitivity", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "SCN9A heterozygous — moderate pain sensitivity", "magnitude": 1},
            "GG": {"status": "normal", "desc": "SCN9A GG — normal pain sensitivity", "magnitude": 0},
        }
    },
    "rs679620": {
        "gene": "MMP3", "category": "Fitness",
        "variants": {
            "AA": {"status": "risk", "desc": "MMP3 AA — increased collagen remodeling, higher tendon/ligament injury risk", "magnitude": 2},
            "AG": {"status": "intermediate", "desc": "MMP3 heterozygous — moderate injury recovery", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "MMP3 heterozygous — moderate injury recovery", "magnitude": 1},
            "GG": {"status": "normal", "desc": "MMP3 GG — normal connective tissue turnover", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: HAIR LOSS / ANDROGENETIC ALOPECIA
    # =========================================================================

    "rs1160312": {
        "gene": "AR/EDA2R", "category": "Hair",
        "variants": {
            "AA": {"status": "high_risk", "desc": "AR/EDA2R AA — strongly associated with male pattern baldness (Xq12 locus)", "magnitude": 3},
            "AG": {"status": "moderate_risk", "desc": "AR/EDA2R heterozygous — moderate hair loss risk", "magnitude": 2},
            "GA": {"status": "moderate_risk", "desc": "AR/EDA2R heterozygous — moderate hair loss risk", "magnitude": 2},
            "GG": {"status": "low_risk", "desc": "AR/EDA2R GG — lower hair loss risk at this locus", "magnitude": 0},
        }
    },
    "rs6152": {
        "gene": "AR", "category": "Hair",
        "variants": {
            "AA": {"status": "risk", "desc": "Androgen receptor AA — associated with androgenetic alopecia sensitivity", "magnitude": 2},
            "AG": {"status": "intermediate", "desc": "AR heterozygous — moderate androgen sensitivity", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "AR heterozygous — moderate androgen sensitivity", "magnitude": 1},
            "GG": {"status": "normal", "desc": "AR GG — typical androgen receptor sensitivity", "magnitude": 0},
        }
    },
    "rs2180439": {
        "gene": "20p11", "category": "Hair",
        "variants": {
            "TT": {"status": "risk", "desc": "20p11 TT — second strongest hair loss locus after AR", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "20p11 heterozygous — moderate hair loss association", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "20p11 heterozygous — moderate hair loss association", "magnitude": 1},
            "CC": {"status": "normal", "desc": "20p11 CC — lower risk at this locus", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: THYROID
    # =========================================================================

    "rs1443434": {
        "gene": "TSHR", "category": "Thyroid",
        "variants": {
            "GG": {"status": "risk", "desc": "TSHR GG — TSH receptor variant, Graves' disease / hyperthyroidism susceptibility", "magnitude": 2},
            "GA": {"status": "intermediate", "desc": "TSHR heterozygous — moderate thyroid autoimmune risk", "magnitude": 1},
            "AG": {"status": "intermediate", "desc": "TSHR heterozygous — moderate thyroid autoimmune risk", "magnitude": 1},
            "AA": {"status": "normal", "desc": "TSHR AA — typical TSH receptor function", "magnitude": 0},
        }
    },
    "rs965513": {
        "gene": "FOXE1/9q22", "category": "Thyroid",
        "variants": {
            "AA": {"status": "risk", "desc": "FOXE1 AA — thyroid cancer and autoimmune thyroid disease susceptibility", "magnitude": 3},
            "AG": {"status": "intermediate", "desc": "FOXE1 heterozygous — moderate thyroid disease risk", "magnitude": 2},
            "GA": {"status": "intermediate", "desc": "FOXE1 heterozygous — moderate thyroid disease risk", "magnitude": 2},
            "GG": {"status": "normal", "desc": "FOXE1 GG — typical thyroid cancer risk", "magnitude": 0},
        }
    },
    # =========================================================================
    # SECTION: GOUT / URIC ACID
    # =========================================================================

    "rs2231142": {
        "gene": "ABCG2", "category": "Gout",
        "variants": {
            "AA": {"status": "high_risk", "desc": "ABCG2 Q141K AA — major urate transporter defect, strong gout risk", "magnitude": 4},
            "AC": {"status": "risk", "desc": "ABCG2 Q141K heterozygous — elevated serum urate and gout risk", "magnitude": 3},
            "CA": {"status": "risk", "desc": "ABCG2 Q141K heterozygous — elevated serum urate and gout risk", "magnitude": 3},
            "CC": {"status": "normal", "desc": "ABCG2 CC — normal urate excretion", "magnitude": 0},
        }
    },
    "rs1014290": {
        "gene": "SLC2A9", "category": "Gout",
        "variants": {
            "TT": {"status": "risk", "desc": "SLC2A9 TT — reduced urate reabsorption regulator, gout susceptibility", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "SLC2A9 heterozygous — moderate urate effect", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "SLC2A9 heterozygous — moderate urate effect", "magnitude": 1},
            "CC": {"status": "normal", "desc": "SLC2A9 CC — typical urate handling", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: EYE HEALTH
    # =========================================================================

    "rs10483727": {
        "gene": "LOXL1", "category": "Eye Health",
        "variants": {
            "CC": {"status": "risk", "desc": "LOXL1 CC — exfoliation glaucoma susceptibility", "magnitude": 3},
            "CT": {"status": "intermediate", "desc": "LOXL1 heterozygous — moderate glaucoma risk", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "LOXL1 heterozygous — moderate glaucoma risk", "magnitude": 2},
            "TT": {"status": "normal", "desc": "LOXL1 TT — lower exfoliation glaucoma risk", "magnitude": 0},
        }
    },
    "rs10811661": {
        "gene": "CDKN2B-AS1", "category": "Eye Health",
        "variants": {
            "TT": {"status": "risk", "desc": "CDKN2B-AS1 TT — primary open-angle glaucoma susceptibility (also T2D risk locus)", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "CDKN2B-AS1 heterozygous — moderate glaucoma risk", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "CDKN2B-AS1 heterozygous — moderate glaucoma risk", "magnitude": 1},
            "CC": {"status": "normal", "desc": "CDKN2B-AS1 CC — typical glaucoma risk", "magnitude": 0},
        }
    },
    "rs12913832": {
        "gene": "HERC2/OCA2", "category": "Eye Health",
        "variants": {
            "AA": {"status": "blue_eyes", "desc": "HERC2 AA — blue eye color, higher UV sensitivity, slightly higher AMD risk", "magnitude": 1},
            "AG": {"status": "mixed", "desc": "HERC2 heterozygous — green/hazel eyes possible", "magnitude": 0},
            "GA": {"status": "mixed", "desc": "HERC2 heterozygous — green/hazel eyes possible", "magnitude": 0},
            "GG": {"status": "brown_eyes", "desc": "HERC2 GG — brown eyes, more UV-protective melanin", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: DENTAL HEALTH
    # =========================================================================

    "rs2274327": {
        "gene": "CA6", "category": "Dental",
        "variants": {
            "TT": {"status": "risk", "desc": "CA6 TT — reduced salivary carbonic anhydrase, higher caries susceptibility", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "CA6 heterozygous — moderate caries risk", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "CA6 heterozygous — moderate caries risk", "magnitude": 1},
            "CC": {"status": "normal", "desc": "CA6 CC — normal salivary buffering capacity", "magnitude": 0},
        }
    },
    "rs17878486": {
        "gene": "AMELX", "category": "Dental",
        "variants": {
            "TT": {"status": "risk", "desc": "AMELX TT — amelogenin variant, enamel formation defect risk", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "AMELX heterozygous — moderate enamel risk", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "AMELX heterozygous — moderate enamel risk", "magnitude": 1},
            "CC": {"status": "normal", "desc": "AMELX CC — normal enamel formation", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: KIDNEY HEALTH
    # =========================================================================

    "rs4293393": {
        "gene": "UMOD", "category": "Kidney",
        "variants": {
            "TT": {"status": "risk", "desc": "UMOD TT — uromodulin variant, chronic kidney disease susceptibility", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "UMOD heterozygous — moderate CKD risk", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "UMOD heterozygous — moderate CKD risk", "magnitude": 1},
            "CC": {"status": "normal", "desc": "UMOD CC — typical uromodulin levels", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: BONE DENSITY
    # =========================================================================

    "rs3736228": {
        "gene": "LRP5", "category": "Bone Health",
        "variants": {
            "TT": {"status": "risk", "desc": "LRP5 TT — reduced Wnt signaling, lower bone mineral density, fracture risk", "magnitude": 3},
            "TC": {"status": "intermediate", "desc": "LRP5 heterozygous — moderate bone density effect", "magnitude": 2},
            "CT": {"status": "intermediate", "desc": "LRP5 heterozygous — moderate bone density effect", "magnitude": 2},
            "CC": {"status": "normal", "desc": "LRP5 CC — normal Wnt signaling and bone density", "magnitude": 0},
        }
    },
    "rs4355801": {
        "gene": "TNFRSF11B", "category": "Bone Health",
        "variants": {
            "AA": {"status": "risk", "desc": "OPG AA — reduced osteoprotegerin, increased bone resorption", "magnitude": 2},
            "AG": {"status": "intermediate", "desc": "OPG heterozygous — moderate bone density effect", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "OPG heterozygous — moderate bone density effect", "magnitude": 1},
            "GG": {"status": "normal", "desc": "OPG GG — normal bone remodeling", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: EXPANDED CARDIOVASCULAR
    # =========================================================================

    "rs1799983": {
        "gene": "NOS3", "category": "Cardiovascular",
        "variants": {
            "TT": {"status": "risk", "desc": "eNOS TT (Glu298Asp) — reduced nitric oxide, endothelial dysfunction, hypertension risk", "magnitude": 3},
            "TG": {"status": "intermediate", "desc": "eNOS heterozygous — moderate endothelial function", "magnitude": 2},
            "GT": {"status": "intermediate", "desc": "eNOS heterozygous — moderate endothelial function", "magnitude": 2},
            "GG": {"status": "normal", "desc": "eNOS GG — normal nitric oxide production", "magnitude": 0},
        }
    },
    "rs1800775": {
        "gene": "CETP", "category": "Cardiovascular",
        "variants": {
            "AA": {"status": "elevated_hdl", "desc": "CETP AA — reduced CETP activity, higher HDL cholesterol (cardioprotective)", "magnitude": 1},
            "AC": {"status": "intermediate", "desc": "CETP heterozygous — moderate HDL effect", "magnitude": 0},
            "CA": {"status": "intermediate", "desc": "CETP heterozygous — moderate HDL effect", "magnitude": 0},
            "CC": {"status": "normal", "desc": "CETP CC — normal HDL cholesterol transfer", "magnitude": 0},
        }
    },
    "rs10455872": {
        "gene": "LPA", "category": "Cardiovascular",
        "variants": {
            "GG": {"status": "high_risk", "desc": "LPA GG — elevated Lp(a), strong independent cardiovascular risk factor", "magnitude": 4},
            "GA": {"status": "risk", "desc": "LPA heterozygous — elevated Lp(a), increased CVD risk", "magnitude": 3},
            "AG": {"status": "risk", "desc": "LPA heterozygous — elevated Lp(a), increased CVD risk", "magnitude": 3},
            "AA": {"status": "normal", "desc": "LPA AA — normal Lp(a) levels", "magnitude": 0},
        },
        "note": "Lp(a) is an independent causal risk factor for heart disease. Test serum Lp(a) if carrier."
    },
    "rs4420638": {
        "gene": "APOC1/APOE", "category": "Cardiovascular",
        "variants": {
            "GG": {"status": "risk", "desc": "APOC1 GG — elevated LDL and triglycerides, linked to APOE e4 haplotype", "magnitude": 2},
            "GA": {"status": "intermediate", "desc": "APOC1 heterozygous — moderate lipid effect", "magnitude": 1},
            "AG": {"status": "intermediate", "desc": "APOC1 heterozygous — moderate lipid effect", "magnitude": 1},
            "AA": {"status": "normal", "desc": "APOC1 AA — typical lipid levels", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: EXPANDED INFLAMMATION & IMMUNE
    # =========================================================================

    "rs1205": {
        "gene": "CRP", "category": "Inflammation",
        "variants": {
            "CC": {"status": "elevated", "desc": "CRP CC — higher baseline C-reactive protein, systemic inflammation marker", "magnitude": 2},
            "CT": {"status": "intermediate", "desc": "CRP heterozygous — moderate CRP tendency", "magnitude": 1},
            "TC": {"status": "intermediate", "desc": "CRP heterozygous — moderate CRP tendency", "magnitude": 1},
            "TT": {"status": "normal", "desc": "CRP TT — lower baseline CRP levels", "magnitude": 0},
        }
    },
    "rs3087243": {
        "gene": "CTLA4", "category": "Autoimmune",
        "variants": {
            "GG": {"status": "risk", "desc": "CTLA4 GG — reduced immune checkpoint, autoimmune disease susceptibility (T1D, thyroid, RA)", "magnitude": 2},
            "GA": {"status": "intermediate", "desc": "CTLA4 heterozygous — moderate autoimmune risk", "magnitude": 1},
            "AG": {"status": "intermediate", "desc": "CTLA4 heterozygous — moderate autoimmune risk", "magnitude": 1},
            "AA": {"status": "normal", "desc": "CTLA4 AA — normal immune checkpoint function", "magnitude": 0},
        }
    },
    "rs2104286": {
        "gene": "IL2RA", "category": "Autoimmune",
        "variants": {
            "AA": {"status": "risk", "desc": "IL2RA AA — reduced IL-2 signaling, T1D and MS susceptibility", "magnitude": 2},
            "AG": {"status": "intermediate", "desc": "IL2RA heterozygous — moderate autoimmune risk", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "IL2RA heterozygous — moderate autoimmune risk", "magnitude": 1},
            "GG": {"status": "normal", "desc": "IL2RA GG — normal regulatory T cell function", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: EXPANDED NUTRITION / METABOLISM
    # =========================================================================

    "rs12325817": {
        "gene": "PEMT", "category": "Nutrition",
        "variants": {
            "CC": {"status": "reduced", "desc": "PEMT CC — reduced phosphatidylcholine synthesis, higher dietary choline requirement", "magnitude": 2},
            "CG": {"status": "intermediate", "desc": "PEMT heterozygous — moderate choline synthesis reduction", "magnitude": 1},
            "GC": {"status": "intermediate", "desc": "PEMT heterozygous — moderate choline synthesis reduction", "magnitude": 1},
            "GG": {"status": "normal", "desc": "PEMT GG — normal endogenous choline production", "magnitude": 0},
        }
    },
    "rs174546": {
        "gene": "FADS1", "category": "Nutrition",
        "variants": {
            "TT": {"status": "low_conversion", "desc": "FADS1 TT — reduced ALA→EPA/DHA conversion, higher need for preformed omega-3", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "FADS1 heterozygous — moderate omega-3 conversion", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "FADS1 heterozygous — moderate omega-3 conversion", "magnitude": 1},
            "CC": {"status": "normal", "desc": "FADS1 CC — efficient omega-3 conversion from plant sources", "magnitude": 0},
        }
    },

    "rs7041": {
        "gene": "GC", "category": "Nutrition",
        "variants": {
            "TT": {"status": "reduced", "desc": "GC TT — Gc1f haplotype, lower vitamin D binding capacity in some populations", "magnitude": 1},
            "TG": {"status": "intermediate", "desc": "GC heterozygous — moderate vitamin D binding", "magnitude": 0},
            "GT": {"status": "intermediate", "desc": "GC heterozygous — moderate vitamin D binding", "magnitude": 0},
            "GG": {"status": "normal", "desc": "GC GG — Gc1s haplotype, typical vitamin D binding", "magnitude": 0},
        }
    },
    "rs492602": {
        "gene": "FUT2", "category": "Nutrition",
        "variants": {
            "GG": {"status": "non_secretor", "desc": "FUT2 GG — non-secretor, ~35% lower B12 absorption from food, altered gut microbiome", "magnitude": 2},
            "GA": {"status": "secretor", "desc": "FUT2 heterozygous — secretor status, normal B12 absorption", "magnitude": 0},
            "AG": {"status": "secretor", "desc": "FUT2 heterozygous — secretor status, normal B12 absorption", "magnitude": 0},
            "AA": {"status": "secretor", "desc": "FUT2 AA — secretor status, normal B12 absorption and gut flora", "magnitude": 0},
        },
        "note": "Non-secretors have different gut microbiome composition and are resistant to norovirus."
    },

    # =========================================================================
    # SECTION: EXPANDED DRUG METABOLISM
    # =========================================================================

    "rs1045642": {
        "gene": "ABCB1", "category": "Drug Metabolism",
        "variants": {
            "TT": {"status": "reduced_transport", "desc": "ABCB1/MDR1 TT (C3435T) — reduced P-glycoprotein efflux, higher drug levels for digoxin, cyclosporine", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "ABCB1 heterozygous — moderate P-glycoprotein activity", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "ABCB1 heterozygous — moderate P-glycoprotein activity", "magnitude": 1},
            "CC": {"status": "normal", "desc": "ABCB1 CC — normal P-glycoprotein drug efflux pump", "magnitude": 0},
        }
    },
    "rs662": {
        "gene": "PON1", "category": "Drug Metabolism",
        "variants": {
            "TT": {"status": "reduced", "desc": "PON1 Q192R TT — reduced paraoxonase, lower organophosphate detoxification, may affect clopidogrel", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "PON1 heterozygous — moderate detoxification capacity", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "PON1 heterozygous — moderate detoxification capacity", "magnitude": 1},
            "CC": {"status": "normal", "desc": "PON1 CC — efficient paraoxonase activity", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: EXPANDED DETOXIFICATION
    # =========================================================================

    "rs1056806": {
        "gene": "CYP1B1", "category": "Detoxification",
        "variants": {
            "CC": {"status": "fast", "desc": "CYP1B1 Leu432Val CC — faster estrogen hydroxylation, potentially more reactive metabolites", "magnitude": 2},
            "CG": {"status": "intermediate", "desc": "CYP1B1 heterozygous — moderate estrogen metabolism", "magnitude": 1},
            "GC": {"status": "intermediate", "desc": "CYP1B1 heterozygous — moderate estrogen metabolism", "magnitude": 1},
            "GG": {"status": "normal", "desc": "CYP1B1 GG — normal estrogen metabolism", "magnitude": 0},
        }
    },
    "rs2606345": {
        "gene": "CYP1A1", "category": "Detoxification",
        "variants": {
            "CC": {"status": "induced", "desc": "CYP1A1 CC — highly inducible by PAHs (grilled/smoked foods, tobacco), more reactive metabolites", "magnitude": 2},
            "CA": {"status": "intermediate", "desc": "CYP1A1 heterozygous — moderate inducibility", "magnitude": 1},
            "AC": {"status": "intermediate", "desc": "CYP1A1 heterozygous — moderate inducibility", "magnitude": 1},
            "AA": {"status": "normal", "desc": "CYP1A1 AA — normal PAH metabolism inducibility", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: EXPANDED FITNESS / ATHLETIC
    # =========================================================================

    "rs1800169": {
        "gene": "CILP", "category": "Fitness",
        "variants": {
            "TT": {"status": "risk", "desc": "CILP TT — cartilage intermediate layer protein variant, disc degeneration risk", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "CILP heterozygous — moderate disc/joint risk", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "CILP heterozygous — moderate disc/joint risk", "magnitude": 1},
            "CC": {"status": "normal", "desc": "CILP CC — normal cartilage integrity", "magnitude": 0},
        }
    },

    "rs1110400": {
        "gene": "MC1R", "category": "Skin",
        "variants": {
            "CC": {"status": "risk", "desc": "MC1R D84E CC — red/fair pigmentation, increased UV damage susceptibility", "magnitude": 2},
            "CT": {"status": "intermediate", "desc": "MC1R D84E carrier — mild fair pigmentation tendency", "magnitude": 1},
            "TC": {"status": "intermediate", "desc": "MC1R D84E carrier — mild fair pigmentation tendency", "magnitude": 1},
            "TT": {"status": "normal", "desc": "MC1R TT — typical melanin at this position", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: EXPANDED HAIR
    # =========================================================================

    "rs4846480": {
        "gene": "7p21.1", "category": "Hair",
        "variants": {
            "TT": {"status": "risk", "desc": "7p21.1 TT — additional male pattern baldness locus, independent of androgen receptor", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "7p21.1 heterozygous — moderate hair loss association", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "7p21.1 heterozygous — moderate hair loss association", "magnitude": 1},
            "CC": {"status": "normal", "desc": "7p21.1 CC — lower risk at this locus", "magnitude": 0},
        }
    },
    "rs12565727": {
        "gene": "TARDBP", "category": "Hair",
        "variants": {
            "CC": {"status": "risk", "desc": "TARDBP CC — hair graying onset, associated with premature graying", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "TARDBP heterozygous — moderate graying timeline", "magnitude": 0},
            "TC": {"status": "intermediate", "desc": "TARDBP heterozygous — moderate graying timeline", "magnitude": 0},
            "TT": {"status": "normal", "desc": "TARDBP TT — typical graying timeline", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: EXPANDED COGNITION / MEMORY
    # =========================================================================

    "rs53576": {
        "gene": "OXTR", "category": "Cognition",
        "variants": {
            "GG": {"status": "empathic", "desc": "Oxytocin receptor GG — higher empathy, social cognition, and stress resilience", "magnitude": 1},
            "GA": {"status": "intermediate", "desc": "OXTR heterozygous — moderate social cognition effect", "magnitude": 0},
            "AG": {"status": "intermediate", "desc": "OXTR heterozygous — moderate social cognition effect", "magnitude": 0},
            "AA": {"status": "reduced", "desc": "OXTR AA — lower oxytocin receptor density, reduced empathic accuracy under stress", "magnitude": 1},
        }
    },

    "rs13146355": {
        "gene": "UMOD", "category": "Kidney",
        "variants": {
            "CC": {"status": "risk", "desc": "UMOD CC — elevated uromodulin, hypertension-mediated kidney damage risk", "magnitude": 2},
            "CT": {"status": "intermediate", "desc": "UMOD heterozygous — moderate kidney risk", "magnitude": 1},
            "TC": {"status": "intermediate", "desc": "UMOD heterozygous — moderate kidney risk", "magnitude": 1},
            "TT": {"status": "normal", "desc": "UMOD TT — normal uromodulin levels", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: EXPANDED THYROID
    # =========================================================================

    "rs3104413": {
        "gene": "HLA-DR", "category": "Thyroid",
        "variants": {
            "CC": {"status": "risk", "desc": "HLA-DR CC — Graves' disease and Hashimoto's thyroiditis susceptibility", "magnitude": 2},
            "CT": {"status": "intermediate", "desc": "HLA-DR heterozygous — moderate thyroid autoimmune risk", "magnitude": 1},
            "TC": {"status": "intermediate", "desc": "HLA-DR heterozygous — moderate thyroid autoimmune risk", "magnitude": 1},
            "TT": {"status": "normal", "desc": "HLA-DR TT — typical thyroid autoimmune risk", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: BLOOD CLOTTING / THROMBOSIS
    # =========================================================================

    "rs8176719": {
        "gene": "ABO", "category": "Blood Clotting",
        "variants": {
            "TT": {"status": "type_O", "desc": "ABO TT — blood type O, 25% lower VTE risk (protective against clots)", "magnitude": 1},
            "TC": {"status": "non_O", "desc": "ABO heterozygous — non-O blood type, higher VTE risk than type O", "magnitude": 1},
            "CT": {"status": "non_O", "desc": "ABO heterozygous — non-O blood type, higher VTE risk than type O", "magnitude": 1},
            "CC": {"status": "non_O", "desc": "ABO CC — non-O blood type, ~25% higher VTE risk", "magnitude": 1},
        },
        "note": "Non-O blood types have 2-4x higher risk of venous thromboembolism."
    },

    # =========================================================================
    # SECTION: MIGRAINE
    # =========================================================================

    "rs10166942": {
        "gene": "TRPM8", "category": "Migraine",
        "variants": {
            "TT": {"status": "protective", "desc": "TRPM8 TT — cold/menthol receptor variant, reduced migraine susceptibility", "magnitude": 1},
            "TC": {"status": "intermediate", "desc": "TRPM8 heterozygous — moderate migraine susceptibility", "magnitude": 0},
            "CT": {"status": "intermediate", "desc": "TRPM8 heterozygous — moderate migraine susceptibility", "magnitude": 0},
            "CC": {"status": "risk", "desc": "TRPM8 CC — typical migraine susceptibility at this locus", "magnitude": 1},
        }
    },
    "rs1835740": {
        "gene": "MTDH", "category": "Migraine",
        "variants": {
            "AA": {"status": "risk", "desc": "MTDH AA — glutamate regulation variant, first GWAS-confirmed migraine locus", "magnitude": 2},
            "AC": {"status": "intermediate", "desc": "MTDH heterozygous — moderate migraine association", "magnitude": 1},
            "CA": {"status": "intermediate", "desc": "MTDH heterozygous — moderate migraine association", "magnitude": 1},
            "CC": {"status": "normal", "desc": "MTDH CC — typical migraine risk at this locus", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: LIVER HEALTH
    # =========================================================================

    "rs738409": {
        "gene": "PNPLA3", "category": "Liver",
        "variants": {
            "GG": {"status": "high_risk", "desc": "PNPLA3 I148M GG — strong non-alcoholic fatty liver disease (NAFLD) and cirrhosis risk", "magnitude": 4},
            "GC": {"status": "risk", "desc": "PNPLA3 I148M heterozygous — elevated NAFLD risk, liver fat accumulation", "magnitude": 3},
            "CG": {"status": "risk", "desc": "PNPLA3 I148M heterozygous — elevated NAFLD risk, liver fat accumulation", "magnitude": 3},
            "CC": {"status": "normal", "desc": "PNPLA3 CC — normal hepatic lipid metabolism", "magnitude": 0},
        },
        "note": "Most common genetic risk factor for fatty liver. Weight loss is the primary intervention."
    },
    "rs58542926": {
        "gene": "TM6SF2", "category": "Liver",
        "variants": {
            "TT": {"status": "high_risk", "desc": "TM6SF2 E167K TT — increased hepatic fat retention, NAFLD/NASH risk", "magnitude": 3},
            "TC": {"status": "risk", "desc": "TM6SF2 heterozygous — moderate liver fat accumulation risk", "magnitude": 2},
            "CT": {"status": "risk", "desc": "TM6SF2 heterozygous — moderate liver fat accumulation risk", "magnitude": 2},
            "CC": {"status": "normal", "desc": "TM6SF2 CC — normal hepatic lipid secretion", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: HEARING
    # =========================================================================

    "rs7598759": {
        "gene": "GRHL2", "category": "Hearing",
        "variants": {
            "TT": {"status": "risk", "desc": "GRHL2 TT — age-related hearing loss susceptibility (presbycusis)", "magnitude": 2},
            "TC": {"status": "intermediate", "desc": "GRHL2 heterozygous — moderate hearing loss risk", "magnitude": 1},
            "CT": {"status": "intermediate", "desc": "GRHL2 heterozygous — moderate hearing loss risk", "magnitude": 1},
            "CC": {"status": "normal", "desc": "GRHL2 CC — typical age-related hearing trajectory", "magnitude": 0},
        }
    },

    # =========================================================================
    # SECTION: FERTILITY / REPRODUCTIVE
    # =========================================================================

    "rs10835638": {
        "gene": "FSHB", "category": "Fertility",
        "variants": {
            "GG": {"status": "reduced", "desc": "FSHB GG — reduced FSH levels, associated with lower sperm count in males", "magnitude": 2},
            "GT": {"status": "intermediate", "desc": "FSHB heterozygous — moderate FSH effect", "magnitude": 1},
            "TG": {"status": "intermediate", "desc": "FSHB heterozygous — moderate FSH effect", "magnitude": 1},
            "TT": {"status": "normal", "desc": "FSHB TT — normal follicle-stimulating hormone levels", "magnitude": 0},
        }
    },
    "rs2736108": {
        "gene": "TERT", "category": "Fertility",
        "variants": {
            "CC": {"status": "risk", "desc": "TERT CC — reduced telomerase in reproductive cells, associated with earlier reproductive aging", "magnitude": 2},
            "CT": {"status": "intermediate", "desc": "TERT heterozygous — moderate reproductive aging effect", "magnitude": 1},
            "TC": {"status": "intermediate", "desc": "TERT heterozygous — moderate reproductive aging effect", "magnitude": 1},
            "TT": {"status": "normal", "desc": "TERT TT — typical reproductive aging trajectory", "magnitude": 0},
        }
    },
}


def _validate_snp_database():
    """Validate all SNP entries at import time."""
    for rsid, info in COMPREHENSIVE_SNPS.items():
        assert "gene" in info, f"{rsid} missing 'gene'"
        assert "category" in info, f"{rsid} missing 'category'"
        assert "variants" in info, f"{rsid} missing 'variants'"
        for gt, variant in info["variants"].items():
            mag = variant.get("magnitude", 0)
            assert 0 <= mag <= 6, f"{rsid}/{gt} has invalid magnitude {mag}"
            assert "status" in variant, f"{rsid}/{gt} missing 'status'"
            assert "desc" in variant, f"{rsid}/{gt} missing 'desc'"

_validate_snp_database()
