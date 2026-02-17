"""Clinical context database and pathway groupings.

Provides deeper interpretation for each gene/status combination,
extracted from the original generate_exhaustive_report.py.
"""


def get_clinical_context(gene, status):
    """Get clinical context for a gene/status combination."""
    return CLINICAL_CONTEXT.get((gene, status), None)


def get_related_pathways(gene):
    """Find all pathways a gene belongs to."""
    return [pathway for pathway, genes in PATHWAYS.items() if gene in genes]


CLINICAL_CONTEXT = {
    # METHYLATION
    ("MTHFR", "significantly_reduced"): {
        "mechanism": "The MTHFR enzyme converts folic acid to methylfolate (5-MTHF), the active form used by the body. With 70% reduced activity, you produce significantly less methylfolate, affecting over 200 methylation-dependent reactions.",
        "implications": [
            "Elevated homocysteine levels (cardiovascular risk marker)",
            "Reduced production of SAMe (S-adenosylmethionine), the universal methyl donor",
            "Potential impacts on neurotransmitter synthesis, DNA repair, and detoxification",
            "Folic acid from fortified foods may accumulate as unmetabolized folic acid (UMFA)"
        ],
        "actions": [
            "Use methylfolate (5-MTHF) instead of folic acid - typical dose 400-800mcg",
            "Consider methylcobalamin (methyl-B12) rather than cyanocobalamin",
            "Support with riboflavin (B2) which is a cofactor for MTHFR",
            "Monitor homocysteine levels periodically",
            "Avoid high-dose folic acid supplements and limit heavily fortified foods"
        ],
        "interactions": ["Synergizes with MTRR status - both impaired compounds effect", "COMT status affects methylation demand"]
    },
    ("MTRR", "significantly_reduced"): {
        "mechanism": "MTRR regenerates methylcobalamin (methyl-B12), the active form of B12. Impaired MTRR means B12 recycling is less efficient, potentially leading to functional B12 deficiency even with normal blood levels.",
        "implications": [
            "May need higher B12 intake to maintain adequate methylcobalamin levels",
            "Can compound MTHFR issues (both affect methylation)",
            "Potential neurological and cognitive effects if B12 recycling insufficient"
        ],
        "actions": [
            "Use methylcobalamin or adenosylcobalamin forms of B12",
            "Consider sublingual or injectable B12 for better absorption",
            "Typical supportive dose: 1000-5000mcg methylcobalamin",
            "Check serum B12 and methylmalonic acid (MMA) for functional status"
        ],
        "interactions": ["Compounds with MTHFR C677T", "Affects homocysteine pathway"]
    },
    ("PEMT", "reduced"): {
        "mechanism": "PEMT creates phosphatidylcholine from phosphatidylethanolamine using SAMe. Reduced function means you rely more on dietary choline intake.",
        "implications": [
            "Higher dietary choline requirement",
            "Particularly important during pregnancy (fetal brain development)",
            "May affect liver health and lipid metabolism"
        ],
        "actions": [
            "Increase dietary choline: eggs (especially yolks), liver, beef, fish",
            "Consider choline supplementation: 250-500mg phosphatidylcholine or CDP-choline",
            "Ensure adequate methyl donors (folate, B12) as PEMT uses SAMe"
        ],
        "interactions": ["Increased demand when MTHFR is impaired (SAMe dependent)"]
    },
    # NEUROTRANSMITTERS
    ("COMT", "slow"): {
        "mechanism": "COMT breaks down catecholamines (dopamine, norepinephrine, epinephrine). Slow COMT means these neurotransmitters remain active longer, leading to higher baseline levels.",
        "implications": [
            "Higher baseline dopamine - better working memory, focus in calm conditions",
            "More sensitive to stress - catecholamines accumulate faster under pressure",
            "Stimulants (caffeine, medications) have stronger, longer-lasting effects",
            "May be more prone to anxiety and rumination",
            "Can have advantages in pain tolerance"
        ],
        "actions": [
            "Prioritize stress management: meditation, breathwork, regular exercise",
            "Use lower doses of stimulants; caffeine half-life effectively longer for you",
            "Magnesium (glycinate, 300-400mg) supports COMT function",
            "Consider adaptogens (ashwagandha, rhodiola) over stimulants for energy",
            "L-theanine may help modulate catecholamine effects",
            "Avoid combining multiple stimulants"
        ],
        "interactions": ["Synergizes with caffeine metabolism genes", "Affects response to ADHD medications"]
    },
    ("OPRM1", "altered"): {
        "mechanism": "OPRM1 encodes the mu-opioid receptor. The A118G variant alters receptor binding, affecting response to endorphins, opioid medications, and reward from alcohol.",
        "implications": [
            "May require different opioid dosing for pain management",
            "Altered reward response to alcohol (some studies show increased craving)",
            "Potentially different response to naltrexone treatment"
        ],
        "actions": [
            "Inform anesthesiologists before surgery",
            "May need adjusted opioid dosing for adequate pain control",
            "Be aware of potentially altered alcohol reward response"
        ],
        "interactions": ["Relevant for pain management planning"]
    },
    # CAFFEINE
    ("CYP1A2", "intermediate"): {
        "mechanism": "CYP1A2 is the primary enzyme metabolizing caffeine. Intermediate metabolizers clear caffeine at a moderate rate, with ~5-6 hour half-life.",
        "implications": [
            "Moderate caffeine clearance - effects last several hours",
            "Afternoon caffeine may affect sleep more than fast metabolizers",
            "Cardiovascular effects of caffeine are intermediate"
        ],
        "actions": [
            "Limit caffeine to morning/early afternoon (before 2pm ideally)",
            "Moderate intake (~200-300mg/day) typically well-tolerated",
            "Wait 90+ minutes after waking for first caffeine (cortisol awakening response)"
        ],
        "interactions": ["ADORA2A affects anxiety response to caffeine"]
    },
    ("ADORA2A", "anxiety_prone"): {
        "mechanism": "ADORA2A encodes adenosine receptors. This variant is associated with increased anxiety response to caffeine.",
        "implications": [
            "More likely to experience jitteriness, anxiety, or panic from caffeine",
            "May be more sensitive to sleep disruption from caffeine"
        ],
        "actions": [
            "Consider lower caffeine doses or slower-release forms",
            "Pair caffeine with L-theanine (green tea naturally has this)",
            "Matcha or green tea may be better tolerated than coffee",
            "Consider caffeine alternatives: yerba mate, guayusa"
        ],
        "interactions": ["Compounds with slow COMT for stress sensitivity"]
    },
    ("ADORA2A", "lower_sensitivity"): {
        "mechanism": "This variant is associated with lower caffeine sensitivity regarding anxiety.",
        "implications": [
            "Less likely to experience anxiety from caffeine",
            "May tolerate higher doses without jitteriness"
        ],
        "actions": [
            "Can likely enjoy caffeine without anxiety issues",
            "Still respect timing for sleep quality"
        ],
        "interactions": []
    },
    # CARDIOVASCULAR
    ("AGTR1", "increased"): {
        "mechanism": "AGTR1 encodes the angiotensin II type 1 receptor. This variant is associated with increased receptor activity and hypertension risk.",
        "implications": [
            "Higher risk of developing hypertension",
            "May have enhanced response to angiotensin receptor blockers (ARBs)"
        ],
        "actions": [
            "Regular blood pressure monitoring",
            "Sodium restriction particularly beneficial",
            "ARBs (losartan, valsartan) may be especially effective if BP medication needed"
        ],
        "interactions": ["Compounds with ACE and AGT variants"]
    },
    ("ACE", "high"): {
        "mechanism": "Higher ACE activity means more conversion of angiotensin I to angiotensin II, a potent vasoconstrictor.",
        "implications": [
            "Increased hypertension risk",
            "May confer advantage in power/sprint athletics",
            "Good response expected to ACE inhibitors"
        ],
        "actions": [
            "Blood pressure monitoring essential",
            "ACE inhibitors (lisinopril, enalapril) likely very effective if needed",
            "Potassium-rich diet supports blood pressure"
        ],
        "interactions": ["Synergizes with AGTR1 and AGT variants for BP risk"]
    },
    ("AGT", "increased"): {
        "mechanism": "AGT M235T is associated with higher angiotensinogen levels, feeding into the renin-angiotensin system.",
        "implications": [
            "Slightly elevated blood pressure risk",
            "Contributes to overall cardiovascular profile"
        ],
        "actions": [
            "Part of overall BP monitoring strategy",
            "Responds to general cardiovascular lifestyle measures"
        ],
        "interactions": ["Additive with ACE and AGTR1 variants"]
    },
    ("GNB3", "increased"): {
        "mechanism": "GNB3 C825T affects G-protein signaling, associated with hypertension and obesity risk.",
        "implications": [
            "Increased hypertension risk",
            "May be more prone to weight gain",
            "Can affect response to certain medications"
        ],
        "actions": [
            "Weight management particularly important",
            "Regular blood pressure monitoring",
            "May respond differently to beta-blockers"
        ],
        "interactions": ["Compounds with other BP variants"]
    },
    # NUTRITION
    ("APOA2", "sensitive"): {
        "mechanism": "APOA2 affects how saturated fat influences body weight. The CC genotype shows strong correlation between saturated fat intake and obesity.",
        "implications": [
            "Saturated fat intake more strongly linked to weight gain for you",
            "Limiting saturated fat may be more impactful than for others"
        ],
        "actions": [
            "Limit saturated fat to <7% of calories",
            "Replace with unsaturated fats (olive oil, nuts, avocado)",
            "Minimize: butter, fatty red meat, full-fat dairy, coconut oil",
            "Prioritize: olive oil, fatty fish, nuts, avocados"
        ],
        "interactions": []
    },
    ("GC", "low"): {
        "mechanism": "GC encodes vitamin D binding protein. This variant results in lower circulating 25-OH vitamin D levels.",
        "implications": [
            "Genetically predisposed to lower vitamin D status",
            "Higher supplementation needs, especially at northern latitudes",
            "May need to supplement year-round"
        ],
        "actions": [
            "Supplement vitamin D3: 2,000-5,000 IU/day depending on season and sun exposure",
            "Take with fat-containing meal for absorption",
            "Test 25-OH vitamin D after 2-3 months",
            "Target blood level: 40-60 ng/mL (100-150 nmol/L)",
            "Consider vitamin K2 (MK-7) alongside D3"
        ],
        "interactions": []
    },
    ("BCMO1", "reduced"): {
        "mechanism": "BCMO1 converts beta-carotene to vitamin A. Reduced activity means less efficient conversion from plant sources.",
        "implications": [
            "Plant carotenoids (carrots, sweet potatoes) less efficiently converted to vitamin A",
            "May benefit from preformed vitamin A sources"
        ],
        "actions": [
            "Include preformed vitamin A: liver, eggs, dairy, fatty fish",
            "Don't rely solely on beta-carotene for vitamin A needs",
            "Still consume carotenoid-rich foods for other benefits (antioxidant)"
        ],
        "interactions": []
    },
    # INFLAMMATION
    ("IL6", "high"): {
        "mechanism": "IL-6 -174 G/G is associated with higher baseline IL-6 production, an inflammatory cytokine.",
        "implications": [
            "Higher baseline inflammation",
            "More pronounced inflammatory response to triggers",
            "May affect recovery, aging, chronic disease risk"
        ],
        "actions": [
            "Anti-inflammatory diet: omega-3s, colorful vegetables, low processed foods",
            "Omega-3 fatty acids (EPA/DHA): 2-3g/day",
            "Regular exercise (but don't overtrain)",
            "Adequate sleep - sleep deprivation spikes IL-6",
            "Consider curcumin, SPMs (specialized pro-resolving mediators)"
        ],
        "interactions": ["Affects recovery and chronic disease risk"]
    },
    # AUTOIMMUNE
    ("HLA-DQA1", "increased_risk"): {
        "mechanism": "HLA-DQ2.5 is strongly associated with celiac disease susceptibility. Not deterministic but increases risk.",
        "implications": [
            "Increased risk of celiac disease (but not guaranteed)",
            "~3% of HLA-DQ2.5 carriers develop celiac disease",
            "Should be aware of celiac symptoms"
        ],
        "actions": [
            "Know celiac symptoms: GI issues, fatigue, anemia, nutrient deficiencies",
            "If symptoms arise, get celiac antibody testing (tTG-IgA) while still eating gluten",
            "No need for preventive gluten-free diet unless symptomatic",
            "Inform healthcare providers of this risk"
        ],
        "interactions": []
    },
    # SKIN
    ("MC1R", "accelerated"): {
        "mechanism": "MC1R variants affect melanin production and skin aging. V92M is associated with accelerated skin aging.",
        "implications": [
            "May show earlier signs of skin aging (wrinkles, photoaging)",
            "Importance of sun protection heightened"
        ],
        "actions": [
            "Daily broad-spectrum SPF 30+ sunscreen",
            "Topical retinoids (tretinoin, retinol) for anti-aging",
            "Antioxidant serums (vitamin C, E)",
            "Avoid excessive sun exposure and tanning"
        ],
        "interactions": []
    },
    # SLEEP
    ("ARNTL", "significantly_altered"): {
        "mechanism": "BMAL1 (ARNTL) is a core circadian clock gene. This variant may weaken circadian rhythm strength.",
        "implications": [
            "May have less robust circadian rhythm",
            "Potentially more susceptible to jet lag, shift work effects",
            "Sleep timing may be less consistent"
        ],
        "actions": [
            "Strong light exposure in morning (10,000 lux or sunlight)",
            "Consistent sleep/wake times, even weekends",
            "Blue light reduction in evening",
            "Consider melatonin 0.5-1mg 30-60 min before bed if needed",
            "Keep bedroom cool, dark, quiet"
        ],
        "interactions": []
    },
    # DETOXIFICATION
    ("NAT2", "intermediate"): {
        "mechanism": "NAT2 acetylates various drugs and toxins. Intermediate acetylators have moderate metabolism of NAT2 substrates.",
        "implications": [
            "Moderate metabolism of drugs like isoniazid, hydralazine, some sulfonamides",
            "Between slow and fast acetylator phenotypes"
        ],
        "actions": [
            "Generally standard drug dosing appropriate",
            "Inform physicians if taking NAT2-metabolized drugs"
        ],
        "interactions": []
    },
    ("SOD2", "high_activity"): {
        "mechanism": "SOD2 (MnSOD) is a key mitochondrial antioxidant enzyme. High activity (Ala/Ala) provides efficient superoxide neutralization.",
        "implications": [
            "Efficient mitochondrial antioxidant defense",
            "Good protection against mitochondrial oxidative stress"
        ],
        "actions": [
            "This is a favorable variant - no specific intervention needed",
            "Continue supporting mitochondrial health: exercise, CoQ10, adequate sleep"
        ],
        "interactions": []
    },
    # FITNESS
    ("ACTN3", "mixed"): {
        "mechanism": "ACTN3 R577X affects fast-twitch muscle fiber composition. R/X (heterozygous) provides balanced fiber distribution.",
        "implications": [
            "Versatile muscle fiber composition",
            "Can develop both power and endurance capacity",
            "Neither extreme sprinter nor ultra-endurance genotype"
        ],
        "actions": [
            "Train for both power and endurance based on goals",
            "Respond well to varied training programs",
            "Can optimize for either direction with proper training"
        ],
        "interactions": []
    },
    ("ADRB2", "gly16"): {
        "mechanism": "ADRB2 Gly16 is associated with enhanced lipolysis (fat burning) response to exercise and catecholamines.",
        "implications": [
            "Better fat mobilization during exercise",
            "Enhanced response to beta-agonists (bronchodilators)"
        ],
        "actions": [
            "May see good fat loss response to exercise",
            "Cardio and HIIT can be particularly effective for body composition"
        ],
        "interactions": []
    },
    # IRON
    ("HFE", "carrier"): {
        "mechanism": "HFE H63D is a minor hemochromatosis variant. Carriers have slightly increased iron absorption.",
        "implications": [
            "Mild increase in iron absorption",
            "Usually not clinically significant alone",
            "Should be aware of iron accumulation over time"
        ],
        "actions": [
            "Periodic ferritin checks (every 1-2 years)",
            "Avoid unnecessary iron supplements",
            "Blood donation if ferritin runs high - helps regulate iron"
        ],
        "interactions": ["Compound with C282Y for hemochromatosis risk"]
    },
    # LONGEVITY
    ("TP53", "arg72"): {
        "mechanism": "TP53 R72P affects p53 apoptotic efficiency. Arg72 has less efficient apoptosis induction.",
        "implications": [
            "Slightly less efficient programmed cell death",
            "Complex effects on cancer and aging"
        ],
        "actions": [
            "Standard cancer screening appropriate for age",
            "Anti-aging lifestyle: exercise, sleep, stress management, healthy diet"
        ],
        "interactions": []
    },
    ("CETP", "favorable"): {
        "mechanism": "CETP I405V affects cholesterol transfer between lipoproteins. This variant is associated with higher HDL and longevity.",
        "implications": [
            "Favorable lipid profile tendency",
            "Associated with exceptional longevity in studies"
        ],
        "actions": [
            "This is a favorable variant",
            "Support with heart-healthy lifestyle"
        ],
        "interactions": []
    },
    # ALCOHOL
    ("ADH1B", "slow"): {
        "mechanism": "ADH1B affects the first step of alcohol metabolism (alcohol to acetaldehyde). Slow metabolizers have longer alcohol effects.",
        "implications": [
            "Alcohol effects last longer",
            "May feel effects at lower doses",
            "Slower acetaldehyde production (not the \"flush\" gene)"
        ],
        "actions": [
            "Moderate alcohol consumption appropriate",
            "Effects may persist longer - factor into timing",
            "Space drinks and stay hydrated"
        ],
        "interactions": ["ALDH2 affects second step (acetaldehyde clearance)"]
    },
    # DRUG METABOLISM
    ("CYP2C19", "rapid"): {
        "mechanism": "CYP2C19*17 causes ultra-rapid metabolism of many drugs including PPIs, some antidepressants, and clopidogrel.",
        "implications": [
            "Faster breakdown of PPIs (may need higher doses or alternatives)",
            "Faster clopidogrel activation (actually better for this prodrug)",
            "Some antidepressants metabolized faster"
        ],
        "actions": [
            "PPIs (omeprazole): may need higher doses or alternatives",
            "Clopidogrel: good metabolizer, likely effective",
            "Inform prescribers about CYP2C19 status"
        ],
        "interactions": ["Pharmacogenomic testing reference"]
    },
    ("CYP3A5", "non_expressor"): {
        "mechanism": "CYP3A5*3/*3 means you don't express CYP3A5. CYP3A4 handles the work. This is the most common genotype in many populations.",
        "implications": [
            "Standard dosing for CYP3A4/5 substrates",
            "Tacrolimus dosing follows standard protocols"
        ],
        "actions": [
            "Standard drug dosing typically appropriate",
            "CYP3A4 remains the primary metabolizer for many drugs"
        ],
        "interactions": []
    },
}

# Gene pathway groupings for cross-referencing
PATHWAYS = {
    "Methylation Cycle": ["MTHFR", "MTRR", "MTR", "CBS", "BHMT", "PEMT"],
    "Catecholamine Metabolism": ["COMT", "MAO-A", "MAO-B", "DBH"],
    "Caffeine Response": ["CYP1A2", "ADORA2A", "ADA"],
    "Blood Pressure": ["ACE", "AGT", "AGTR1", "GNB3", "ADRB1"],
    "Inflammation": ["IL6", "TNF", "CRP", "IL10"],
    "Drug Metabolism - Phase I": ["CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4", "CYP3A5"],
    "Drug Metabolism - Phase II": ["NAT2", "GSTP1", "UGT1A1"],
    "Lipid Metabolism": ["APOE", "APOA2", "CETP", "PPARG", "PPARA"],
    "Vitamin Metabolism": ["GC", "BCMO1", "FUT2", "MTHFR"],
    "Circadian Rhythm": ["ARNTL", "PER2", "CLOCK"],
    "Muscle & Exercise": ["ACTN3", "PPARGC1A", "ADRB2", "ACE"],
}
