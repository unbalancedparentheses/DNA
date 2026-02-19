"""Research-backed genomic insights.

Provides curated, referenced interpretations for single-gene findings and
multi-gene narrative patterns.  Each entry links to published research
(PMIDs) and includes practical, personalized context.
"""


# =========================================================================
# SINGLE-GENE INSIGHTS
# =========================================================================
# Keyed by (gene, status).  Each value is a list of insight dicts.

SINGLE_GENE_INSIGHTS = {
    # --- Pharmacogenomics ---
    ("CYP2C9", "intermediate"): [
        {
            "title": "CYP2C9 intermediate metabolism and warfarin sensitivity",
            "finding": (
                "CYP2C9 *1/*2 carriers require ~20% lower warfarin doses on average. "
                "The *2 allele reduces enzyme activity to ~12% of normal for that copy, "
                "slowing metabolism of S-warfarin, NSAIDs like ibuprofen, and phenytoin."
            ),
            "reference": "Johnson et al. 2017, Clin Pharmacol Ther, PMID: 28198005",
            "practical": (
                "If prescribed warfarin, your genotype predicts lower dose requirements. "
                "Standard starting doses may cause over-anticoagulation. NSAIDs should "
                "be used at the lowest effective dose."
            ),
        },
    ],
    ("CYP2C9", "poor"): [
        {
            "title": "CYP2C9 poor metabolism — significant drug impact",
            "finding": (
                "CYP2C9 poor metabolizers have substantially reduced clearance of "
                "warfarin, NSAIDs, and several other drugs. Warfarin dose reductions "
                "of 30-50% are typical."
            ),
            "reference": "Johnson et al. 2017, Clin Pharmacol Ther, PMID: 28198005",
            "practical": (
                "This is a clinically significant finding. Carry a pharmacogenomic "
                "card listing your CYP2C9 status for any prescribing encounter."
            ),
        },
    ],
    ("DPYD", "intermediate"): [
        {
            "title": "DPYD and fluoropyrimidine chemotherapy safety",
            "finding": (
                "DPYD encodes dihydropyrimidine dehydrogenase, which metabolizes >80% "
                "of 5-FU. Intermediate metabolizers have increased risk of severe "
                "toxicity (grade 3-4) from standard-dose fluoropyrimidines."
            ),
            "reference": "Meulendijks et al. 2015, J Clin Oncol, PMID: 29152729",
            "practical": (
                "If you ever need 5-FU or capecitabine chemotherapy, your oncologist "
                "should start at 50% dose and titrate up based on tolerance. Pre-treatment "
                "DPYD testing is now recommended by CPIC guidelines."
            ),
        },
    ],
    ("DPYD", "poor"): [
        {
            "title": "DPYD deficiency — critical pharmacogenomic finding",
            "finding": (
                "DPYD poor metabolizers are at high risk of life-threatening toxicity "
                "from fluoropyrimidine chemotherapy. Standard doses can cause fatal "
                "myelosuppression and mucositis."
            ),
            "reference": "Meulendijks et al. 2015, J Clin Oncol, PMID: 29152729",
            "practical": (
                "Fluoropyrimidines (5-FU, capecitabine) should be avoided or used at "
                "drastically reduced doses with extreme caution. Ensure this is in "
                "your medical record."
            ),
        },
    ],
    ("CYP1A2", "slow"): [
        {
            "title": "Slow caffeine metabolism and cardiovascular risk",
            "finding": (
                "Slow CYP1A2 metabolizers who drink 2+ cups of coffee/day have a "
                "significantly increased risk of myocardial infarction compared to "
                "fast metabolizers. Caffeine half-life is ~8-10 hours vs ~4 hours."
            ),
            "reference": "Cornelis et al. 2006, JAMA, PMID: 16522833",
            "practical": (
                "Limiting coffee to 1 cup/day before 10am is prudent. Afternoon "
                "caffeine will still be circulating at bedtime, disrupting sleep "
                "architecture even if you fall asleep fine."
            ),
        },
    ],
    ("CYP1A2", "intermediate"): [
        {
            "title": "Moderate caffeine metabolism",
            "finding": (
                "Intermediate CYP1A2 metabolizers clear caffeine at a moderate rate "
                "(half-life ~5-6 hours). The cardiovascular risk from coffee is "
                "intermediate between slow and fast metabolizers."
            ),
            "reference": "Cornelis et al. 2006, JAMA, PMID: 16522833",
            "practical": (
                "Moderate caffeine intake (1-2 cups before early afternoon) is "
                "reasonable. Pay attention to sleep quality as a personal barometer."
            ),
        },
    ],

    # --- Longevity / Protection ---
    ("APOE", "e2_carrier"): [
        {
            "title": "APOE e2 and cardiovascular protection",
            "finding": (
                "APOE e2 carriers have ~20% lower LDL cholesterol on average and "
                "reduced risk of coronary artery disease. The e2 allele is associated "
                "with enhanced clearance of atherogenic lipoproteins."
            ),
            "reference": "Bennet et al. 2007, J Biol Chem, PMID: 17142461",
            "practical": (
                "Your e2 allele is associated with favorable lipid profiles. This "
                "is a protective finding, though standard cardiovascular prevention "
                "still applies."
            ),
        },
        {
            "title": "APOE e2 and reduced Alzheimer's risk",
            "finding": (
                "APOE e2 carriers have roughly half the Alzheimer's disease risk "
                "compared to the common e3/e3 genotype. The e2 allele is associated "
                "with less amyloid-beta aggregation."
            ),
            "reference": "Corder et al. 1994, Nat Genet, PMID: 7894489",
            "practical": (
                "This is protective but not absolute. Maintain cardiovascular health, "
                "regular exercise, quality sleep, and cognitive engagement — all of "
                "which further reduce dementia risk."
            ),
        },
    ],
    ("CETP", "favorable"): [
        {
            "title": "CETP I405V and longevity association",
            "finding": (
                "The CETP I405V variant (VV genotype) is associated with larger HDL "
                "particle size, higher HDL cholesterol, and was significantly enriched "
                "in Ashkenazi Jewish centenarians."
            ),
            "reference": "Barzilai et al. 2003, JAMA, PMID: 18235225",
            "practical": (
                "This is a favorable variant associated with cardiovascular protection "
                "and longevity. Support with a heart-healthy lifestyle."
            ),
        },
    ],
    ("CETP", "highly_favorable"): [
        {
            "title": "CETP homozygous longevity variant",
            "finding": (
                "Homozygous CETP I405V carriers have the strongest association with "
                "elevated HDL particle size and longevity. This genotype was found at "
                "2-3x higher frequency in centenarians."
            ),
            "reference": "Barzilai et al. 2003, JAMA, PMID: 18235225",
            "practical": (
                "This is one of the strongest single-gene longevity signals. Your "
                "lipid profile likely trends favorable."
            ),
        },
    ],

    # --- Blood Pressure ---
    ("AGTR1", "increased"): [
        {
            "title": "AGTR1 and angiotensin receptor blocker response",
            "finding": (
                "The AGTR1 A1166C variant is associated with enhanced angiotensin II "
                "receptor activity and hypertension. Carriers show particularly good "
                "blood pressure response to ARBs (losartan, valsartan, irbesartan)."
            ),
            "reference": "Bonnardeaux et al. 1994, Hypertension, PMID: 16129819",
            "practical": (
                "If you need blood pressure medication, ARBs target the exact "
                "receptor your variant overactivates — a genetically informed choice."
            ),
        },
    ],
    ("ACE", "high"): [
        {
            "title": "ACE I/D polymorphism and exercise blood pressure response",
            "finding": (
                "The ACE DD genotype is associated with higher ACE activity, greater "
                "exercise-induced blood pressure rise, but also better strength/power "
                "performance. Regular aerobic exercise reduces the BP impact."
            ),
            "reference": "Myerson et al. 1999, J Appl Physiol, PMID: 10973849",
            "practical": (
                "Consistent aerobic exercise is especially important for your "
                "genotype. If medication is needed, ACE inhibitors target your "
                "specific pathway."
            ),
        },
    ],

    # --- Nutrition ---
    ("GC", "low"): [
        {
            "title": "GC variants and vitamin D status at different latitudes",
            "finding": (
                "GC encodes vitamin D binding protein. Your variant reduces circulating "
                "25-OH vitamin D levels. At latitudes above ~35N (e.g. most of the US, "
                "Europe), winter UVB is insufficient, compounding genetic predisposition."
            ),
            "reference": "Wang et al. 2010, Lancet, PMID: 20541252",
            "practical": (
                "Supplement vitamin D3 (2,000-5,000 IU/day) year-round, especially "
                "in winter. Target 40-60 ng/mL. Take with a fat-containing meal."
            ),
        },
    ],
    ("PEMT", "reduced"): [
        {
            "title": "PEMT and dietary choline requirements",
            "finding": (
                "PEMT creates phosphatidylcholine using SAMe as a methyl donor. "
                "Reduced PEMT function increases dietary choline needs. 75% of the "
                "population does not meet adequate intake for choline (550mg men, "
                "425mg women)."
            ),
            "reference": "da Costa et al. 2006, Am J Clin Nutr, PMID: 17510091",
            "practical": (
                "Prioritize choline-rich foods: eggs (147mg/egg), beef liver "
                "(356mg/3oz), salmon, chicken. Choline is critical for liver "
                "health, brain function, and cell membranes."
            ),
        },
    ],
    ("PEMT", "significantly_reduced"): [
        {
            "title": "PEMT homozygous — high choline dependency",
            "finding": (
                "Homozygous PEMT variants substantially reduce endogenous "
                "phosphatidylcholine synthesis. This genotype developed organ "
                "dysfunction on low-choline diets in clinical studies."
            ),
            "reference": "da Costa et al. 2006, Am J Clin Nutr, PMID: 17510091",
            "practical": (
                "Dietary choline is essential for you. Consider supplementation "
                "(CDP-choline 250-500mg) if dietary intake is insufficient. "
                "Particularly important during pregnancy."
            ),
        },
    ],
    ("BCMO1", "reduced"): [
        {
            "title": "BCMO1 and beta-carotene conversion efficiency",
            "finding": (
                "BCMO1 converts plant beta-carotene to retinol (vitamin A). Your "
                "variant reduces conversion efficiency by ~30-70%. Relying solely on "
                "plant carotenoids for vitamin A may lead to subclinical deficiency."
            ),
            "reference": "Leung et al. 2009, FASEB J, PMID: 19103647",
            "practical": (
                "Include preformed vitamin A sources: eggs, liver, dairy, fatty fish. "
                "Don't rely on carrots and sweet potatoes alone for vitamin A."
            ),
        },
    ],

    # --- Fitness ---
    ("ACTN3", "endurance"): [
        {
            "title": "ACTN3 R577X and muscle fiber composition",
            "finding": (
                "ACTN3 XX genotype means absence of alpha-actinin-3 in fast-twitch "
                "muscle fibers. This shifts muscle composition toward slow-twitch "
                "(endurance). XX is found in ~18% of Europeans but is rare in elite "
                "sprinters."
            ),
            "reference": "Yang et al. 2003, Am J Hum Genet, PMID: 18043716",
            "practical": (
                "Your genetics favor endurance over sprint/power activities. You may "
                "find you excel at distance running, cycling, or swimming. Strength "
                "training still builds muscle, just with a different fiber emphasis."
            ),
        },
    ],
    ("ACTN3", "power"): [
        {
            "title": "ACTN3 RR and sprint/power performance",
            "finding": (
                "ACTN3 RR genotype preserves full alpha-actinin-3 in fast-twitch "
                "fibers. This genotype is over-represented in elite sprinters and "
                "power athletes across multiple populations."
            ),
            "reference": "Yang et al. 2003, Am J Hum Genet, PMID: 18043716",
            "practical": (
                "Your fast-twitch fiber composition favors explosive activities. "
                "Include plyometrics, sprints, and heavy compound lifts."
            ),
        },
    ],
    ("ACTN3", "mixed"): [
        {
            "title": "ACTN3 heterozygous — versatile muscle profile",
            "finding": (
                "ACTN3 RX gives a balanced fast/slow-twitch fiber mix. Heterozygotes "
                "can develop either power or endurance based on training focus."
            ),
            "reference": "Yang et al. 2003, Am J Hum Genet, PMID: 18043716",
            "practical": (
                "You have the most versatile muscle genetics — adapt training to "
                "your goals. Both strength and endurance respond well."
            ),
        },
    ],

    # --- Iron ---
    ("HFE", "carrier"): [
        {
            "title": "HFE H63D and iron regulation",
            "finding": (
                "HFE H63D is a common variant (~15% carrier frequency in Europeans) "
                "that mildly increases iron absorption. Alone it rarely causes clinical "
                "hemochromatosis, but iron can accumulate over decades."
            ),
            "reference": "Feder et al. 1996, Nat Genet; EASL Guidelines, PMID: 20471131",
            "practical": (
                "Monitor ferritin every 1-2 years. Avoid iron supplements unless "
                "deficiency is confirmed. Regular blood donation helps regulate iron."
            ),
        },
    ],

    # --- Taste / Behavior ---
    ("TAS2R38", "non_taster"): [
        {
            "title": "TAS2R38 non-taster and food preferences",
            "finding": (
                "TAS2R38 non-tasters cannot detect bitter compounds (PTC/PROP) in "
                "cruciferous vegetables. This may reduce aversion to bitter greens "
                "but is also associated with higher preference for sweet/salty foods."
            ),
            "reference": "Kim et al. 2003, Science, PMID: 16051168",
            "practical": (
                "You likely find broccoli, kale, and Brussels sprouts palatable. "
                "Be mindful of the tendency toward sweet/salty preferences."
            ),
        },
    ],

    # --- Skin ---
    ("MC1R", "accelerated"): [
        {
            "title": "MC1R variants and perceived age",
            "finding": (
                "MC1R variants are associated with looking ~2 years older than "
                "chronological age, independent of sun exposure and skin color. "
                "The effect is mediated through melanin pathway impacts on skin "
                "structure."
            ),
            "reference": "Liu et al. 2016, Curr Biol, PMID: 27189542",
            "practical": (
                "Consistent sun protection (SPF 30+) and topical retinoids can "
                "partially offset the photoaging acceleration."
            ),
        },
    ],

    # --- Alcohol ---
    ("ADH1B", "slow"): [
        {
            "title": "ADH1B and alcohol metabolism patterns",
            "finding": (
                "ADH1B slow metabolizers convert alcohol to acetaldehyde more slowly. "
                "Alcohol effects last longer per drink. This is distinct from the "
                "ALDH2 'Asian flush' variant."
            ),
            "reference": "Edenberg 2007, Alcohol Res Health, PMID: 17718403",
            "practical": (
                "Pace drinks more slowly. Effects persist longer, so factor timing "
                "into driving and next-day planning."
            ),
        },
    ],
    ("ADH1B", "fast"): [
        {
            "title": "ADH1B fast metabolism and alcohol response",
            "finding": (
                "ADH1B fast metabolizers rapidly convert alcohol to acetaldehyde. "
                "This can cause flushing and nausea with alcohol, which is actually "
                "protective against alcohol use disorder."
            ),
            "reference": "Edenberg 2007, Alcohol Res Health, PMID: 17718403",
            "practical": (
                "Your aversive response to alcohol is genetically protective. "
                "Don't push through discomfort — it's a signal."
            ),
        },
    ],

    # --- Methylation ---
    ("MTHFR", "reduced"): [
        {
            "title": "MTHFR C677T heterozygous and folate metabolism",
            "finding": (
                "One copy of C677T reduces MTHFR enzyme activity by ~35%. This is "
                "very common (~40% of many populations carry at least one copy) and "
                "has modest effects on homocysteine that are easily managed."
            ),
            "reference": "Frosst et al. 1995, Nat Genet; Klerk et al. 2002, PMID: 12145525",
            "practical": (
                "Use methylfolate (5-MTHF) instead of folic acid when supplementing. "
                "Eat folate-rich foods: leafy greens, legumes, citrus."
            ),
        },
    ],
    ("MTHFR", "significantly_reduced"): [
        {
            "title": "MTHFR C677T homozygous — 70% reduced activity",
            "finding": (
                "Homozygous C677T reduces MTHFR activity by ~70%, significantly "
                "impairing conversion of folic acid to the active 5-MTHF form. "
                "Homocysteine levels are typically elevated without intervention."
            ),
            "reference": "Frosst et al. 1995, Nat Genet; Klerk et al. 2002, PMID: 12145525",
            "practical": (
                "Methylfolate supplementation (400-800mcg) is recommended. Monitor "
                "homocysteine levels. Avoid synthetic folic acid in fortified foods."
            ),
        },
    ],

    # --- Autoimmune ---
    ("HLA-DQA1", "increased_risk"): [
        {
            "title": "HLA-DQ2.5 and celiac disease susceptibility",
            "finding": (
                "HLA-DQ2.5 is present in ~95% of celiac patients, but only ~3% of "
                "carriers develop the disease. It is necessary but not sufficient."
            ),
            "reference": "Sollid et al. 1989, J Exp Med; Lundin et al. 2015",
            "practical": (
                "No need for preventive gluten-free diet. If GI symptoms, fatigue, "
                "or unexplained iron deficiency develop, request tTG-IgA testing "
                "while still eating gluten."
            ),
        },
    ],

    # --- Neurotransmitters ---
    ("COMT", "slow"): [
        {
            "title": "COMT Val158Met and stress-cognition tradeoff",
            "finding": (
                "Slow COMT (Met/Met) provides higher baseline dopamine in prefrontal "
                "cortex — better working memory and focus in calm conditions, but "
                "vulnerability to stress-induced cognitive decline (the 'worrier' "
                "vs 'warrior' model)."
            ),
            "reference": "Egan et al. 2001, PNAS, PMID: 11381111",
            "practical": (
                "Leverage calm-state cognitive advantages. Build robust stress "
                "management (meditation, exercise). Stimulants have amplified effects."
            ),
        },
    ],
    ("COMT", "fast"): [
        {
            "title": "COMT Val/Val — the 'warrior' genotype",
            "finding": (
                "Fast COMT rapidly clears dopamine and catecholamines. This provides "
                "stress resilience and steady performance under pressure, but lower "
                "baseline prefrontal dopamine in calm conditions."
            ),
            "reference": "Egan et al. 2001, PNAS, PMID: 11381111",
            "practical": (
                "You perform well under pressure. May benefit from mild stimulants "
                "(caffeine) for focus during routine tasks."
            ),
        },
    ],
}


# =========================================================================
# MULTI-GENE NARRATIVE PATTERNS
# =========================================================================

MULTI_GENE_NARRATIVES = [
    {
        "id": "bp_cluster",
        "title": "Your Blood Pressure Genetic Profile",
        "required_genes": {"AGTR1": ["elevated", "increased"],
                          "ACE": ["high"]},
        "optional_genes": {"AGT": ["increased"], "GNB3": ["increased"],
                          "ADRB1": ["arg389"]},
        "min_matches": 2,
        "narrative": (
            "You carry variants in {matched_count} genes of the renin-angiotensin-"
            "aldosterone system (RAAS), the body's primary blood pressure control "
            "pathway. Each variant independently nudges BP upward; together they "
            "compound the effect. Your AGTR1 variant overactivates the angiotensin "
            "receptor, while ACE DD increases angiotensin II production."
        ),
        "references": [
            "Bonnardeaux et al. 1994, Hypertension, PMID: 8021009",
            "Rigat et al. 1990, J Clin Invest, PMID: 10973849",
        ],
        "practical": (
            "ARBs (losartan, valsartan) target the exact receptor your AGTR1 variant "
            "overactivates. ACE inhibitors address the upstream enzyme. Sodium "
            "restriction (<2g/day) and regular aerobic exercise have outsized benefit "
            "for RAAS-pathway hypertension."
        ),
    },
    {
        "id": "caffeine_profile",
        "title": "Your Caffeine Response Profile",
        "required_genes": {"CYP1A2": ["slow", "intermediate"]},
        "optional_genes": {"COMT": ["slow"], "ADORA2A": ["anxiety_prone"]},
        "min_matches": 1,
        "narrative": (
            "Your caffeine response is shaped by {matched_count} genetic variants. "
            "CYP1A2 determines how fast you clear caffeine from your bloodstream. "
            "Slow/intermediate metabolism means caffeine lingers longer, amplifying "
            "both alertness and side effects (jitters, insomnia, anxiety)."
        ),
        "references": [
            "Cornelis et al. 2006, JAMA, PMID: 16522833",
        ],
        "practical": (
            "Limit to 1-2 cups before noon. Green tea provides gentler caffeine "
            "with L-theanine to smooth the curve. If you also carry COMT slow, "
            "the catecholamine amplification makes caffeine particularly stimulating."
        ),
    },
    {
        "id": "methylation_choline",
        "title": "Methylation & Choline Needs",
        "required_genes": {"PEMT": ["reduced", "significantly_reduced"]},
        "optional_genes": {"MTHFR": ["reduced", "significantly_reduced"],
                          "MTRR": ["reduced", "significantly_reduced"],
                          "COMT": ["slow"]},
        "min_matches": 1,
        "narrative": (
            "Your methylation cycle has {matched_count} variant(s) affecting methyl "
            "donor availability and choline synthesis. PEMT creates phosphatidylcholine "
            "using SAMe — reduced function means higher dietary choline dependency."
        ),
        "references": [
            "da Costa et al. 2006, Am J Clin Nutr, PMID: 17510091",
            "Frosst et al. 1995, Nat Genet, PMID: 12145525",
        ],
        "practical": (
            "Prioritize choline (eggs, liver, fish) and methylfolate. If MTHFR is "
            "also reduced, SAMe production drops further, increasing PEMT's reliance "
            "on dietary choline. Slow COMT means start methyl donors at low doses."
        ),
    },
    {
        "id": "iron_profile",
        "title": "Your Iron Metabolism Profile",
        "required_genes": {"HFE": ["carrier", "at_risk", "high_risk"]},
        "optional_genes": {},
        "min_matches": 1,
        "narrative": (
            "You carry {matched_count} HFE variant(s) that increase intestinal iron "
            "absorption. Over decades, excess iron accumulates in the liver, heart, "
            "and pancreas. Historically, this may have been advantageous in iron-poor "
            "diets, but modern iron-fortified foods and red meat availability reverse "
            "the benefit."
        ),
        "references": [
            "EASL Clinical Practice Guidelines, PMID: 20471131",
        ],
        "practical": (
            "Monitor ferritin annually. Avoid iron supplements unless proven "
            "deficient. Regular blood donation is therapeutic and community-beneficial."
        ),
    },
    {
        "id": "cardiovascular_lipids",
        "title": "Your Cardiovascular Lipid Profile",
        "required_genes": {"CETP": ["favorable", "highly_favorable"]},
        "optional_genes": {"APOE": ["e2_carrier", "protective"],
                          "APOA2": ["sensitive"]},
        "min_matches": 1,
        "narrative": (
            "Your lipid genetics include {matched_count} variant(s) influencing "
            "cholesterol transport. The CETP longevity variant promotes larger, "
            "more protective HDL particles."
        ),
        "references": [
            "Barzilai et al. 2003, JAMA, PMID: 18235225",
        ],
        "practical": (
            "Your genetic lipid profile trends favorable. Support with "
            "Mediterranean diet, regular exercise, and omega-3 fatty acids."
        ),
    },
    {
        "id": "vitamin_d_profile",
        "title": "Your Vitamin D Metabolism",
        "required_genes": {"GC": ["low"]},
        "optional_genes": {"BCMO1": ["reduced"]},
        "min_matches": 1,
        "narrative": (
            "Your GC variant reduces vitamin D binding protein levels, lowering "
            "circulating 25-OH vitamin D. If you also have BCMO1 reduced, "
            "beta-carotene conversion to vitamin A is impaired — a double hit "
            "on fat-soluble vitamin availability."
        ),
        "references": [
            "Wang et al. 2010, Lancet, PMID: 20541252",
        ],
        "practical": (
            "Supplement D3 year-round (2,000-5,000 IU). Test 25-OH-D levels. "
            "If BCMO1 is also reduced, include preformed vitamin A sources."
        ),
    },
    {
        "id": "skin_aging",
        "title": "Your Skin & Sun Sensitivity Profile",
        "required_genes": {"MC1R": ["accelerated", "fair_skin", "red_hair",
                                    "high_risk"]},
        "optional_genes": {},
        "min_matches": 1,
        "narrative": (
            "MC1R variants reduce eumelanin (dark pigment) production, increasing "
            "UV sensitivity and photoaging. These variants are associated with "
            "looking ~2 years older and elevated melanoma risk."
        ),
        "references": [
            "Liu et al. 2016, Curr Biol, PMID: 27189542",
        ],
        "practical": (
            "Daily SPF 30+ is essential. Annual dermatology screening. Topical "
            "retinoids can partially offset accelerated photoaging."
        ),
    },
    {
        "id": "inflammation_profile",
        "title": "Your Inflammatory Response Profile",
        "required_genes": {"IL6": ["high"]},
        "optional_genes": {"SOD2": ["reduced", "high_activity"]},
        "min_matches": 1,
        "narrative": (
            "IL-6 is a key inflammatory cytokine. Your variant produces higher "
            "baseline IL-6 levels, driving a more pronounced inflammatory response "
            "to triggers like poor sleep, stress, and processed foods."
        ),
        "references": [
            "Fishman et al. 1998, J Clin Invest",
        ],
        "practical": (
            "Anti-inflammatory diet (omega-3s, colorful vegetables, turmeric), "
            "quality sleep, and regular exercise powerfully counteract this "
            "predisposition."
        ),
    },
]


# =========================================================================
# INSIGHT GENERATOR
# =========================================================================

def generate_insights(lifestyle_findings, apoe=None, star_alleles=None,
                      ancestry_results=None, epistasis_results=None,
                      disease_findings=None):
    """Generate research-backed genomic insights from analysis results.

    Returns dict with:
        single_gene: matched single-gene insight entries
        narratives: matched multi-gene narrative stories
        genome_highlights: top 5 most interesting/unique findings
        protective_findings: good-news items with research context
    """
    # Build gene -> finding lookup
    gene_status = {}
    gene_mag = {}
    for f in (lifestyle_findings or []):
        gene = f.get("gene", "")
        status = f.get("status", "")
        mag = f.get("magnitude", 0)
        if gene and status:
            if gene not in gene_status or mag > gene_mag.get(gene, 0):
                gene_status[gene] = status
                gene_mag[gene] = mag

    # Add APOE status
    if apoe and apoe.get("apoe_type", "Unknown") != "Unknown":
        apoe_type = apoe["apoe_type"]
        if "e2" in apoe_type:
            gene_status["APOE"] = "e2_carrier"
            gene_mag["APOE"] = 2

    # Add star allele statuses
    if star_alleles:
        star_status_map = {
            "poor": "poor", "intermediate": "intermediate",
            "rapid": "rapid", "ultrarapid": "ultrarapid",
        }
        for gene, r in star_alleles.items():
            phenotype = r.get("phenotype", "")
            if phenotype in star_status_map:
                gene_status[gene] = star_status_map[phenotype]
                gene_mag[gene] = 3

    # --- Single-gene insights ---
    single_gene = []
    for (gene, status), entries in SINGLE_GENE_INSIGHTS.items():
        if gene_status.get(gene) == status:
            for entry in entries:
                single_gene.append({
                    **entry,
                    "gene": gene,
                    "status": status,
                    "magnitude": gene_mag.get(gene, 0),
                })

    single_gene.sort(key=lambda x: -x["magnitude"])

    # --- Multi-gene narratives ---
    narratives = []
    for pattern in MULTI_GENE_NARRATIVES:
        matched_genes = []

        # Check required genes
        for gene, valid_statuses in pattern["required_genes"].items():
            if gene_status.get(gene) in valid_statuses:
                matched_genes.append(gene)

        # Check optional genes
        for gene, valid_statuses in pattern["optional_genes"].items():
            if gene_status.get(gene) in valid_statuses:
                matched_genes.append(gene)

        if len(matched_genes) >= pattern["min_matches"]:
            # Check at least one required gene is present
            has_required = any(
                gene_status.get(g) in statuses
                for g, statuses in pattern["required_genes"].items()
            )
            if not has_required:
                continue

            try:
                narrative_text = pattern["narrative"].format(
                    matched_count=len(matched_genes))
            except KeyError:
                narrative_text = pattern["narrative"]
            practical_text = pattern["practical"]

            narratives.append({
                "id": pattern["id"],
                "title": pattern["title"],
                "matched_genes": matched_genes,
                "narrative": narrative_text,
                "references": pattern["references"],
                "practical": practical_text,
            })

    # --- Genome highlights (top 5 most interesting) ---
    highlights = []

    # APOE e2 is protective — interesting
    if apoe and "e2" in apoe.get("apoe_type", ""):
        highlights.append({
            "title": f"APOE {apoe['apoe_type']} — Alzheimer's protection",
            "detail": "APOE e2 carriers have roughly half the Alzheimer's risk.",
            "type": "protective",
        })

    # CETP longevity
    if gene_status.get("CETP") in ("favorable", "highly_favorable"):
        highlights.append({
            "title": "CETP longevity variant",
            "detail": "Associated with exceptional longevity and favorable HDL.",
            "type": "protective",
        })

    # Pathogenic findings
    if disease_findings:
        for f in disease_findings.get("pathogenic", [])[:2]:
            gene = f.get("gene", "Unknown")
            traits = f.get("traits", "")
            first_trait = traits.split(";")[0].strip() if traits else "Pathogenic"
            highlights.append({
                "title": f"{gene} — {first_trait}",
                "detail": f"Pathogenic variant detected ({f.get('gold_stars', 0)}/4 stars).",
                "type": "clinical",
            })

    # Actionable star alleles
    if star_alleles:
        for gene, r in star_alleles.items():
            if r.get("phenotype") in ("poor", "intermediate") and len(highlights) < 5:
                phenotype = r["phenotype"].replace("_", " ").title()
                highlights.append({
                    "title": f"{gene} {r['diplotype']} — {phenotype} Metabolizer",
                    "detail": r.get("clinical_note", ""),
                    "type": "pharmacogenomic",
                })

    # High-magnitude lifestyle findings
    high_mag = sorted(
        [f for f in (lifestyle_findings or []) if f.get("magnitude", 0) >= 3],
        key=lambda x: -x["magnitude"]
    )
    for f in high_mag:
        if len(highlights) >= 5:
            break
        gene = f.get("gene", "")
        # Skip if already covered
        if any(gene in h.get("title", "") for h in highlights):
            continue
        highlights.append({
            "title": f"{gene} — {f.get('status', '').replace('_', ' ').title()}",
            "detail": f.get("description", ""),
            "type": "lifestyle",
        })

    highlights = highlights[:5]

    # --- Protective findings ---
    protective = []
    protective_statuses = {"protective", "longevity", "optimal", "favorable",
                          "highly_favorable", "fast", "e2_carrier"}
    for gene, status in gene_status.items():
        if status in protective_statuses:
            entries = SINGLE_GENE_INSIGHTS.get((gene, status), [])
            if entries:
                protective.append({
                    "gene": gene,
                    "status": status,
                    "title": entries[0]["title"],
                    "finding": entries[0]["finding"],
                    "reference": entries[0]["reference"],
                })

    return {
        "single_gene": single_gene,
        "narratives": narratives,
        "genome_highlights": highlights,
        "protective_findings": protective,
    }
