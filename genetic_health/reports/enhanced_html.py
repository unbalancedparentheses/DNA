"""
Interactive All-in-One Genetic Health Report Generator — Redesigned

12-section layout: actionable information first, clinical details last.
ELI5 plain-language explanations throughout. Paper links on every finding.

Output: reports/GENETIC_HEALTH_REPORT.html
  - Zero external dependencies (all CSS/JS/SVG inline)
  - Dark/light mode via prefers-color-scheme
  - Collapsible sections (vanilla JS)
  - Print-optimized doctor card
  - Database links for every rsID
  - Curated paper references
"""

import html as html_mod
import json
import math
import sys
from datetime import datetime
from collections import defaultdict
from pathlib import Path

from ..config import REPORTS_DIR


# =============================================================================
# COLOR PALETTE — single source of truth for all inline SVG/badge colors
# =============================================================================

C = {
    # Semantic
    "red": "#bf4040",       # warm red — danger, high risk, pathogenic
    "orange": "#c47a2b",    # burnt orange — likely pathogenic, elevated
    "amber": "#c49a20",     # golden — moderate, warning
    "green": "#3a8a5c",     # sage — good, low risk, protective
    "blue": "#4a7da5",      # steel — info, average, drug response
    "purple": "#7d69ac",    # lavender — secondary accent
    "slate": "#8993a0",     # cool gray — muted, informational
    # Chart palette (12 distinguishable, muted tones)
    "chart": [
        "#bf4040", "#c49a20", "#3a8a5c", "#4a7da5", "#7d69ac",
        "#b5577d", "#2d8b7a", "#c47a2b", "#5b68a8", "#6b9c3d",
        "#3d8f9f", "#a85252",
    ],
    # Ancestry
    "eur": "#4a7da5", "afr": "#c49a20", "eas": "#bf4040",
    "sas": "#7d69ac", "amr": "#3a8a5c",
    # Phenotype metabolizers
    "poor": "#bf4040", "intermediate": "#c49a20", "normal": "#3a8a5c",
    "rapid": "#4a7da5", "ultrarapid": "#7d69ac", "unknown_pheno": "#8993a0",
}


def _esc(text):
    """Escape text for safe HTML interpolation."""
    return html_mod.escape(str(text)) if text else ""


def _clean_condition(raw):
    """Clean ClinVar pipe-separated condition text into readable form.

    e.g. 'PI S|Alpha-1-antitrypsin deficiency|not provided' -> 'Alpha-1-antitrypsin deficiency'
    e.g. 'HYPERTENSION, DIASTOLIC, RESISTANCE TO|KCNMB1-related disorder' -> 'Hypertension, diastolic resistance'
    """
    if not raw:
        return "Unknown"
    # Split on pipe and semicolons, filter out junk
    parts = []
    for segment in raw.replace("|", ";").split(";"):
        segment = segment.strip()
        if not segment:
            continue
        low = segment.lower()
        # Skip generic/unhelpful segments
        if low in ("not provided", "not specified", "unknown", "see cases"):
            continue
        if low.endswith("-related disorder") and len(parts) > 0:
            continue  # Skip generic "GENE-related disorder" if we have a real name
        parts.append(segment)
    if not parts:
        return raw.split("|")[0].split(";")[0].strip()
    # Pick the most readable one (prefer longer, lowercase-ish names over ALL CAPS)
    best = parts[0]
    for p in parts:
        if not p.isupper() and len(p) > len(best) // 2:
            best = p
            break
    # Clean ALL CAPS into title case
    if best.isupper():
        best = best.replace(",", ", ").title()
        # Clean up common patterns
        best = best.replace(" To", "").replace("Resistance", "resistance").replace("Susceptibility", "susceptibility")
    return best


def _dedup_phrases(text):
    """Deduplicate repeated semicolon-separated phrases while preserving order.

    e.g. 'ClinVar risk factors related to X; ClinVar risk factors related to X' -> 'ClinVar risk factors related to X'
    """
    if not text or ";" not in text:
        return text
    seen = set()
    unique = []
    for phrase in text.split("; "):
        phrase = phrase.strip()
        if phrase and phrase not in seen:
            seen.add(phrase)
            unique.append(phrase)
    return "; ".join(unique)


def _clean_why(text):
    """Clean pipe-separated ClinVar text from 'why' fields and deduplicate."""
    if not text:
        return text
    if "|" in text:
        text = _clean_condition(text)
    return _dedup_phrases(text)


# =============================================================================
# PAPER REFERENCES — curated, hardcoded
# =============================================================================

PAPER_REFS = {
    "rs1799853": [
        {"title": "CPIC Guideline for Pharmacogenetics-Guided Warfarin Dosing",
         "pmid": "28198005", "year": 2017},
    ],
    "rs1057910": [
        {"title": "CPIC Guideline for Pharmacogenetics-Guided Warfarin Dosing",
         "pmid": "28198005", "year": 2017},
    ],
    "rs5186": [
        {"title": "Association of AGTR1 A1166C polymorphism with hypertension",
         "pmid": "8021009", "year": 1994},
    ],
    "rs28929474": [
        {"title": "ATS/ERS Statement on Alpha-1 Antitrypsin Deficiency",
         "pmid": "14680078", "year": 2003},
    ],
    "rs762551": [
        {"title": "Coffee, CYP1A2 genotype, and risk of myocardial infarction (JAMA)",
         "pmid": "16522833", "year": 2006},
    ],
    "rs7946": [
        {"title": "Choline deficiency: a common variant in PEMT (AJCN)",
         "pmid": "16400055", "year": 2006},
    ],
    "rs1800562": [
        {"title": "EASL Clinical Practice Guidelines for HFE Hemochromatosis",
         "pmid": "20471131", "year": 2010},
    ],
    "rs1799945": [
        {"title": "EASL Clinical Practice Guidelines for HFE Hemochromatosis",
         "pmid": "20471131", "year": 2010},
    ],
    "rs1815739": [
        {"title": "ACTN3 genotype is associated with human elite athletic performance (AJHG)",
         "pmid": "12900424", "year": 2003},
    ],
    "rs429358": [
        {"title": "APOE e4 and risk of Alzheimer's disease (Science)",
         "pmid": "8346443", "year": 1993},
    ],
    "rs2298383": [
        {"title": "Caffeine-induced anxiety and ADORA2A polymorphism",
         "pmid": "12624725", "year": 2003},
    ],
    "rs2228479": [
        {"title": "MC1R variants and skin cancer risk (Nature Genetics)",
         "pmid": "7581446", "year": 1995},
    ],
    "rs1805007": [
        {"title": "MC1R variants and skin cancer risk (Nature Genetics)",
         "pmid": "7581446", "year": 1995},
    ],
    "rs1801133": [
        {"title": "MTHFR C677T polymorphism and cardiovascular disease risk",
         "pmid": "12145525", "year": 2002},
    ],
    "rs4680": [
        {"title": "COMT Val158Met and cognitive function",
         "pmid": "11381111", "year": 2001},
    ],
    "rs9923231": [
        {"title": "VKORC1 pharmacogenetics and warfarin dose requirements",
         "pmid": "15930419", "year": 2005},
    ],
    "rs6025": [
        {"title": "Factor V Leiden and risk of venous thromboembolism",
         "pmid": "7989540", "year": 1994},
    ],
    "rs2802292": [
        {"title": "FOXO3A genotype and human longevity",
         "pmid": "18765803", "year": 2008},
    ],
    "rs699": [
        {"title": "AGT M235T and essential hypertension",
         "pmid": "1528933", "year": 1992},
    ],
}


# =============================================================================
# DATABASE LINK BUILDER
# =============================================================================

def db_links_html(rsid):
    """Generate HTML links to external databases for a given rsID."""
    if not rsid or not rsid.startswith("rs"):
        return ""
    links = [
        f'<a href="https://www.ncbi.nlm.nih.gov/snp/{rsid}" target="_blank" rel="noopener">dbSNP</a>',
        f'<a href="https://www.ncbi.nlm.nih.gov/clinvar/?term={rsid}" target="_blank" rel="noopener">ClinVar</a>',
        f'<a href="https://www.snpedia.com/index.php/{rsid}" target="_blank" rel="noopener">SNPedia</a>',
        f'<a href="https://www.pharmgkb.org/search?query={rsid}" target="_blank" rel="noopener">PharmGKB</a>',
    ]
    return '<span class="db-links">' + " &middot; ".join(links) + "</span>"


def paper_refs_html(rsid):
    """Generate HTML for paper references for a given rsID."""
    refs = PAPER_REFS.get(rsid, [])
    if not refs:
        return ""
    items = []
    for r in refs:
        items.append(
            f'<a href="https://pubmed.ncbi.nlm.nih.gov/{r["pmid"]}/" '
            f'target="_blank" rel="noopener">{r["title"]} ({r["year"]})</a>'
        )
    return '<div class="paper-refs">References: ' + " | ".join(items) + "</div>"


# =============================================================================
# DATA LOADING
# =============================================================================

def load_json(path):
    if not path.exists():
        return {}
    with open(path, "r") as f:
        return json.load(f)


# =============================================================================
# SVG CHART GENERATORS
# =============================================================================

def svg_impact_bar(findings):
    """Horizontal bar chart: impact distribution."""
    high = sum(1 for f in findings if f.get("magnitude", 0) >= 3)
    mod = sum(1 for f in findings if f.get("magnitude", 0) == 2)
    low = sum(1 for f in findings if f.get("magnitude", 0) == 1)
    info = sum(1 for f in findings if f.get("magnitude", 0) == 0)

    max_val = max(high, mod, low, info, 1)
    bar_w = 280

    def bar(y, val, color, label):
        w = int(val / max_val * bar_w) if max_val else 0
        return (
            f'<rect x="110" y="{y}" width="{w}" height="28" rx="4" fill="{color}" opacity="0.85"/>'
            f'<text x="105" y="{y+19}" text-anchor="end" '
            f'fill="currentColor" font-size="13">{label}</text>'
            f'<text x="{115+w}" y="{y+19}" fill="currentColor" '
            f'font-size="13" font-weight="bold">{val}</text>'
        )

    bars = [
        bar(10, high, C["red"], "High"),
        bar(48, mod, C["amber"], "Moderate"),
        bar(86, low, C["green"], "Low"),
        bar(124, info, C["slate"], "Info"),
    ]
    return (
        '<svg viewBox="0 0 440 165" class="chart" role="img" '
        'aria-label="Impact distribution bar chart">'
        f'<title>Impact Distribution</title>{"".join(bars)}</svg>'
    )


def svg_category_donut(findings):
    """Donut chart: findings by category."""
    counts = defaultdict(int)
    for f in findings:
        counts[f.get("category", "Other")] += 1

    if not counts:
        return ""

    colors = C["chart"]
    total = sum(counts.values())
    cx, cy, r = 120, 120, 90
    inner_r = 55
    angle = -90
    paths = []
    legend_items = []

    for idx, (cat, cnt) in enumerate(sorted(counts.items(), key=lambda x: -x[1])):
        color = colors[idx % len(colors)]
        sweep = cnt / total * 360
        start_rad = math.radians(angle)
        end_rad = math.radians(angle + sweep)

        x1 = cx + r * math.cos(start_rad)
        y1 = cy + r * math.sin(start_rad)
        x2 = cx + r * math.cos(end_rad)
        y2 = cy + r * math.sin(end_rad)
        ix1 = cx + inner_r * math.cos(start_rad)
        iy1 = cy + inner_r * math.sin(start_rad)
        ix2 = cx + inner_r * math.cos(end_rad)
        iy2 = cy + inner_r * math.sin(end_rad)

        large = 1 if sweep > 180 else 0
        d = (
            f"M {ix1:.1f} {iy1:.1f} "
            f"L {x1:.1f} {y1:.1f} "
            f"A {r} {r} 0 {large} 1 {x2:.1f} {y2:.1f} "
            f"L {ix2:.1f} {iy2:.1f} "
            f"A {inner_r} {inner_r} 0 {large} 0 {ix1:.1f} {iy1:.1f} Z"
        )
        paths.append(f'<path d="{d}" fill="{color}" opacity="0.85"/>')
        angle += sweep

        ly = 15 + idx * 22
        legend_items.append(
            f'<rect x="260" y="{ly}" width="14" height="14" rx="3" fill="{color}"/>'
            f'<text x="280" y="{ly+12}" fill="currentColor" font-size="12">'
            f'{cat} ({cnt})</text>'
        )

    center_text = (
        f'<text x="{cx}" y="{cy - 5}" text-anchor="middle" '
        f'fill="currentColor" font-size="24" font-weight="bold">{total}</text>'
        f'<text x="{cx}" y="{cy + 15}" text-anchor="middle" '
        f'fill="currentColor" font-size="11">findings</text>'
    )

    height = max(250, 15 + len(counts) * 22 + 10)
    return (
        f'<svg viewBox="0 0 440 {height}" class="chart" role="img" '
        f'aria-label="Category donut chart">'
        f'<title>Findings by Category</title>'
        f'{"".join(paths)}{center_text}{"".join(legend_items)}</svg>'
    )


def svg_metabolism_gauge(label, level, color):
    """Simple gauge for drug metabolism speed. level: 0-2 (slow/intermediate/fast)."""
    cx, cy, r = 70, 65, 50
    start_math = 150
    end_math = 30
    range_math = start_math - end_math

    ptr_math = start_math - (level / 2) * range_math
    ptr_rad = math.radians(ptr_math)
    ptr_r = r - 10
    px = cx + ptr_r * math.cos(ptr_rad)
    py = cy - ptr_r * math.sin(ptr_rad)

    s_rad = math.radians(start_math)
    e_rad = math.radians(end_math)
    sx = cx + r * math.cos(s_rad)
    sy = cy - r * math.sin(s_rad)
    ex = cx + r * math.cos(e_rad)
    ey = cy - r * math.sin(e_rad)

    speed_labels = {0: "Slow", 1: "Intermediate", 2: "Fast"}
    speed_text = speed_labels.get(level, "Unknown")

    return (
        f'<svg viewBox="0 0 140 100" class="gauge" role="img" '
        f'aria-label="{label} metabolism gauge">'
        f'<title>{label}: {speed_text}</title>'
        f'<path d="M {sx:.1f} {sy:.1f} A {r} {r} 0 0 1 {ex:.1f} {ey:.1f}" '
        f'fill="none" stroke="var(--border)" stroke-width="8" stroke-linecap="round"/>'
        f'<circle cx="{px:.1f}" cy="{py:.1f}" r="6" fill="{color}"/>'
        f'<text x="15" y="78" fill="currentColor" font-size="8">Slow</text>'
        f'<text x="108" y="78" fill="currentColor" font-size="8">Fast</text>'
        f'<text x="{cx}" y="95" text-anchor="middle" fill="currentColor" '
        f'font-size="11" font-weight="bold">{label}: {speed_text}</text>'
        f'</svg>'
    )


def svg_ancestry_donut(ancestry_results):
    """Donut chart showing ancestry proportions."""
    if not ancestry_results or not ancestry_results.get("proportions"):
        return ""

    proportions = ancestry_results["proportions"]
    colors = {"EUR": C["eur"], "AFR": C["afr"], "EAS": C["eas"],
              "SAS": C["sas"], "AMR": C["amr"]}
    labels = {"EUR": "European", "AFR": "African", "EAS": "East Asian",
              "SAS": "South Asian", "AMR": "Admixed American"}

    cx, cy, r = 110, 110, 85
    inner_r = 50
    angle = -90
    paths = []
    legend_items = []

    sorted_pops = sorted(proportions.items(), key=lambda x: -x[1])

    for idx, (pop, prop) in enumerate(sorted_pops):
        if prop < 0.005:
            continue
        color = colors.get(pop, C["slate"])
        sweep = prop * 360
        start_rad = math.radians(angle)
        end_rad = math.radians(angle + sweep)

        x1 = cx + r * math.cos(start_rad)
        y1 = cy + r * math.sin(start_rad)
        x2 = cx + r * math.cos(end_rad)
        y2 = cy + r * math.sin(end_rad)
        ix1 = cx + inner_r * math.cos(start_rad)
        iy1 = cy + inner_r * math.sin(start_rad)
        ix2 = cx + inner_r * math.cos(end_rad)
        iy2 = cy + inner_r * math.sin(end_rad)

        large = 1 if sweep > 180 else 0
        d = (
            f"M {ix1:.1f} {iy1:.1f} "
            f"L {x1:.1f} {y1:.1f} "
            f"A {r} {r} 0 {large} 1 {x2:.1f} {y2:.1f} "
            f"L {ix2:.1f} {iy2:.1f} "
            f"A {inner_r} {inner_r} 0 {large} 0 {ix1:.1f} {iy1:.1f} Z"
        )
        paths.append(f'<path d="{d}" fill="{color}" opacity="0.85"/>')
        angle += sweep

        ly = 10 + idx * 24
        label = labels.get(pop, pop)
        legend_items.append(
            f'<rect x="240" y="{ly}" width="14" height="14" rx="3" fill="{color}"/>'
            f'<text x="260" y="{ly+12}" fill="currentColor" font-size="12">'
            f'{label} ({prop:.1%})</text>'
        )

    top = ancestry_results.get("top_ancestry", "")
    center = (
        f'<text x="{cx}" y="{cy - 5}" text-anchor="middle" '
        f'fill="currentColor" font-size="14" font-weight="bold">{top}</text>'
        f'<text x="{cx}" y="{cy + 15}" text-anchor="middle" '
        f'fill="currentColor" font-size="11">{ancestry_results.get("confidence","")}</text>'
    )

    height = max(230, 10 + len(sorted_pops) * 24 + 10)
    return (
        f'<svg viewBox="0 0 420 {height}" class="chart" role="img" '
        f'aria-label="Ancestry proportion donut chart">'
        f'<title>Ancestry Proportions</title>'
        f'{"".join(paths)}{center}{"".join(legend_items)}</svg>'
    )


def svg_prs_gauge(label, percentile, category):
    """Semi-circular gauge for PRS percentile with 4 color zones."""
    cx, cy, r = 80, 70, 55

    zones = [
        (0.0, 0.2, C["green"]),
        (0.2, 0.8, C["blue"]),
        (0.8, 0.95, C["amber"]),
        (0.95, 1.0, C["red"]),
    ]

    zone_paths = []
    for start_frac, end_frac, color in zones:
        a1 = math.radians(180 - start_frac * 180)
        a2 = math.radians(180 - end_frac * 180)
        x1 = cx + r * math.cos(a1)
        y1 = cy - r * math.sin(a1)
        x2 = cx + r * math.cos(a2)
        y2 = cy - r * math.sin(a2)
        large = 1 if abs(end_frac - start_frac) > 0.5 else 0
        zone_paths.append(
            f'<path d="M {x1:.1f} {y1:.1f} A {r} {r} 0 {large} 1 {x2:.1f} {y2:.1f}" '
            f'fill="none" stroke="{color}" stroke-width="10" stroke-linecap="butt"/>'
        )

    frac = max(0.0, min(1.0, percentile / 100.0))
    ptr_angle = math.radians(180 - frac * 180)
    ptr_r = r - 18
    px = cx + ptr_r * math.cos(ptr_angle)
    py = cy - ptr_r * math.sin(ptr_angle)

    cat_colors = {"low": C["green"], "average": C["blue"],
                  "elevated": C["amber"], "high": C["red"]}
    ptr_color = cat_colors.get(category, C["slate"])

    return (
        f'<svg viewBox="0 0 160 105" class="prs-gauge" role="img" '
        f'aria-label="{label} PRS gauge">'
        f'<title>{label}: {percentile:.0f}th percentile ({category})</title>'
        f'{"".join(zone_paths)}'
        f'<circle cx="{px:.1f}" cy="{py:.1f}" r="6" fill="{ptr_color}"/>'
        f'<text x="{cx}" y="{cy + 5}" text-anchor="middle" fill="currentColor" '
        f'font-size="16" font-weight="bold">{percentile:.0f}%</text>'
        f'<text x="{cx}" y="{cy + 20}" text-anchor="middle" fill="currentColor" '
        f'font-size="10">{category.title()}</text>'
        f'<text x="{cx}" y="100" text-anchor="middle" fill="currentColor" '
        f'font-size="11" font-weight="bold">{label}</text>'
        f'</svg>'
    )


# =============================================================================
# ELI5 HELPERS — plain language for every gene/condition
# =============================================================================

ELI5_GENES = {
    "MTHFR": "This gene helps your body use B-vitamins (like folate). A variant here means you may need a special form of folate called methylfolate.",
    "COMT": "This gene controls how fast your brain clears stress chemicals (like adrenaline). A slow version means caffeine and stress hit you harder, but you may also think more carefully.",
    "CYP2C9": "This gene controls how fast your liver breaks down certain medicines (like warfarin and ibuprofen). You may need lower doses.",
    "CYP2C19": "This gene affects how you process common drugs like antidepressants and stomach acid pills. Your doctor should know your type.",
    "CYP2D6": "This gene affects how you process painkillers (codeine), antidepressants, and other drugs. Copy number variants can't be detected from saliva tests.",
    "CYP1A2": "This gene controls how fast you break down caffeine. Slow metabolizers should limit coffee to mornings.",
    "APOE": "This gene affects cholesterol transport in your brain and blood. The e4 version raises Alzheimer's risk, but exercise and heart health can help a lot.",
    "AGTR1": "This gene is part of your blood pressure control system. A variant here means salt and stress may raise your blood pressure more than average.",
    "AGT": "Another blood pressure gene. Together with AGTR1, it tells us your blood pressure system needs extra attention.",
    "SERPINA1": "This gene makes a protein that protects your lungs. A variant means your lungs have less protection — avoid smoke and air pollution.",
    "HFE": "This gene controls iron absorption. A variant means your body may absorb too much iron — get your iron levels checked.",
    "FOXO3": "A longevity gene! Certain variants are found more often in people who live past 100.",
    "ACTN3": "The 'speed gene'. One version is found in most Olympic sprinters; the other favors endurance sports.",
    "F5": "Factor V Leiden — a blood clotting gene. A variant raises your risk of blood clots, especially during surgery or long flights.",
    "VKORC1": "This gene affects warfarin (blood thinner) dosing. If you ever need warfarin, your doctor should test this gene.",
    "DPYD": "This gene processes certain chemotherapy drugs. A deficiency can be life-threatening — always test before fluoropyrimidine chemo.",
    "TPMT": "This gene processes immune-suppressing drugs (like azathioprine). Deficiency causes dangerous side effects.",
    "UGT1A1": "This gene processes bilirubin. The *28 variant causes Gilbert syndrome (harmless yellowing) but affects some drug doses.",
    "BRCA1": "A major cancer-protection gene. Harmful variants significantly raise breast and ovarian cancer risk. Screening is critical.",
    "BRCA2": "Like BRCA1, this gene repairs DNA damage. Harmful variants raise breast, ovarian, and other cancer risks.",
    "MC1R": "The 'red hair gene'. Variants affect hair color, skin type, and sun sensitivity.",
    "HERC2": "The main gene controlling eye color. Variants here are the biggest predictor of blue vs. brown eyes.",
    "ABCC11": "Controls earwax type (wet vs. dry) and body odor. The dry type is common in East Asian populations.",
}

ELI5_CONDITIONS = {
    "blood_pressure": "Your genes suggest your blood pressure system runs a bit high. Think of it like a garden hose with the pressure turned up slightly. Eating less salt, exercising, and monitoring can help a lot.",
    "methylation": "Methylation is like your body's maintenance crew — it repairs DNA, makes brain chemicals, and detoxifies. Your genes suggest the crew could use some help (methylfolate supplements).",
    "iron_metabolism": "Your body's iron thermostat may be set a bit high. Too much iron can damage organs over time, so regular blood tests are important.",
    "caffeine": "Your body is slower at breaking down caffeine. That means coffee keeps you wired longer and may affect sleep. Stick to mornings only.",
    "clotting": "Your blood may clot a bit more easily than average. Stay hydrated, move on long flights, and tell surgeons about this.",
    "drug_metabolism": "Your liver processes some medicines differently. This doesn't mean anything is wrong — but your doctor should know so they can pick the right dose.",
}


def _eli5_for_gene(gene):
    """Get a plain-language explanation for a gene."""
    return ELI5_GENES.get(gene, "")


def _eli5_for_condition(condition_id):
    """Get a plain-language explanation for a condition."""
    return ELI5_CONDITIONS.get(condition_id, "")


# =============================================================================
# SECTION 1: YOUR KEY FINDINGS
# =============================================================================

def build_key_findings(findings, recommendations_data, apoe_data, acmg_data,
                       prs_results, star_alleles_data, longevity_data,
                       traits_data=None, blood_type_data=None, sleep_data=None,
                       mental_health_data=None, ancestry_data=None,
                       nutrigenomics_data=None, disease_findings_data=None):
    """Build comprehensive plain-language summary of ALL findings.

    This is the executive summary — everything important, no jargon.
    Each subsequent section provides the supporting detail.
    """
    parts = []
    parts.append(
        '<p class="eli5">Your body has a recipe book called DNA. We read yours and '
        'here is everything important we found, in plain language. '
        'Each topic below links to a detailed section later in the report.</p>'
    )

    # ---- URGENT: things that need action ----
    urgent_bullets = []

    acmg_findings = (acmg_data or {}).get("acmg_findings", [])
    for af in acmg_findings:
        gene = af.get("gene", "Unknown")
        # Clean ClinVar pipe-separated text: take first meaningful segment
        raw = (af.get("traits") or "Unknown")
        condition = _clean_condition(raw)
        eli5 = _eli5_for_gene(gene)
        text = f'<strong>{gene}</strong>: You carry a variant linked to <em>{condition}</em> that doctors consider medically actionable. You should see a genetic counselor.'
        if eli5:
            text += f'<br><span class="eli5-inline">{eli5}</span>'
        urgent_bullets.append(text)

    priorities = (recommendations_data or {}).get("priorities", [])
    for p in priorities:
        if p["priority"] == "high":
            eli5 = _eli5_for_condition(p.get("id", ""))
            # Clean pipe-separated ClinVar text from recommendation reasons
            why = _clean_why(p["why"])
            text = f'<strong>{p["title"]}</strong>: {why}'
            if eli5:
                text += f'<br><span class="eli5-inline">{eli5}</span>'
            urgent_bullets.append(text)

    if urgent_bullets:
        parts.append('<h3>Act On These <a href="#action-plan">&rarr; details</a></h3>')
        parts.append('<ul class="key-findings">')
        for b in urgent_bullets:
            parts.append(f'<li class="kf-red">{b}</li>')
        parts.append('</ul>')

    # ---- HEALTH RISKS: APOE, PRS, disease findings ----
    risk_bullets = []

    if apoe_data and apoe_data.get("apoe_type", "Unknown") != "Unknown":
        risk = apoe_data["risk_level"]
        color = "red" if risk in ("high", "elevated") else "yellow" if risk == "moderate" else "green"
        eli5 = _eli5_for_gene("APOE")
        text = f'<strong>APOE {apoe_data["apoe_type"]}</strong>: {risk.title()} Alzheimer\'s risk (odds ratio: {apoe_data.get("alzheimer_or", "N/A")}x).'
        if eli5:
            text += f'<br><span class="eli5-inline">{eli5}</span>'
        risk_bullets.append((color, text))

    if prs_results:
        elevated_prs = [(cid, r) for cid, r in prs_results.items()
                        if r["risk_category"] in ("elevated", "high")]
        average_prs = [(cid, r) for cid, r in prs_results.items()
                       if r["risk_category"] == "average"]
        low_prs = [(cid, r) for cid, r in prs_results.items()
                   if r["risk_category"] == "low"]
        for cid, r in elevated_prs:
            risk_bullets.append(("yellow",
                f'<strong>{r["name"]}</strong>: {r["percentile"]:.0f}th percentile genetic risk — higher than {r["percentile"]:.0f}% of people.'))
        if average_prs:
            names = ", ".join(r["name"] for _, r in average_prs)
            risk_bullets.append(("green", f'Average genetic risk for: {names}.'))
        if low_prs:
            names = ", ".join(r["name"] for _, r in low_prs)
            risk_bullets.append(("green", f'<em>Lower</em> than average risk for: {names}.'))

    # Pathogenic ClinVar findings count
    if disease_findings_data:
        path_count = len(disease_findings_data.get("pathogenic", []))
        lp_count = len(disease_findings_data.get("likely_pathogenic", []))
        prot_count = len(disease_findings_data.get("protective", []))
        if path_count or lp_count:
            risk_bullets.append(("yellow",
                f'{path_count} pathogenic and {lp_count} likely pathogenic variants found in ClinVar scan.'))
        if prot_count:
            risk_bullets.append(("green",
                f'{prot_count} protective variants detected — these <em>lower</em> your risk for certain diseases.'))

    if risk_bullets:
        parts.append('<h3>Health Risks <a href="#disease-risk">&rarr; details</a></h3>')
        parts.append('<ul class="key-findings">')
        for color, text in risk_bullets:
            parts.append(f'<li class="kf-{color}">{text}</li>')
        parts.append('</ul>')

    # ---- DRUGS & MEDICATIONS ----
    drug_bullets = []
    if star_alleles_data:
        non_normal = [(gene, r) for gene, r in star_alleles_data.items()
                      if r["phenotype"] not in ("normal", "Unknown")]
        normal = [(gene, r) for gene, r in star_alleles_data.items()
                  if r["phenotype"] == "normal"]
        for gene, r in non_normal:
            eli5 = _eli5_for_gene(gene)
            text = f'<strong>{gene} ({r["diplotype"]})</strong>: {r["phenotype"].replace("_"," ").title()} Metabolizer — some drugs need dose adjustment.'
            if eli5:
                text += f'<br><span class="eli5-inline">{eli5}</span>'
            drug_bullets.append(("yellow", text))
        if normal:
            names = ", ".join(g for g, _ in normal)
            drug_bullets.append(("green", f'Normal metabolism for: {names} — standard drug doses should work.'))

    if drug_bullets:
        parts.append('<h3>Drugs &amp; Medications <a href="#drug-guide">&rarr; details</a></h3>')
        parts.append('<ul class="key-findings">')
        for color, text in drug_bullets:
            parts.append(f'<li class="kf-{color}">{text}</li>')
        parts.append('</ul>')

    # ---- NUTRITION & SUPPLEMENTS ----
    # Skip items already shown in "Act On These" (high priority)
    urgent_titles = {p["title"] for p in priorities if p["priority"] == "high"}
    nutrition_bullets = []
    for p in priorities:
        if p["priority"] in ("high", "moderate") and p["title"] not in urgent_titles:
            nutrition_groups = {"methylation", "iron", "caffeine", "nutrition", "vitamin"}
            if any(g in p.get("id", "").lower() for g in nutrition_groups):
                eli5 = _eli5_for_condition(p.get("id", ""))
                why = _clean_why(p["why"])
                text = f'<strong>{p["title"]}</strong>: {why}'
                if eli5:
                    text += f'<br><span class="eli5-inline">{eli5}</span>'
                color = "yellow"
                nutrition_bullets.append((color, text))

    if nutrigenomics_data:
        supps = nutrigenomics_data.get("supplement_priorities", [])
        if supps:
            top_supps = ", ".join(f'{s["nutrient"]} ({s["form"]})' for s in supps[:4])
            nutrition_bullets.append(("yellow",
                f'<strong>Top supplement priorities</strong>: {top_supps}.'))

    if nutrition_bullets:
        parts.append('<h3>Nutrition &amp; Supplements <a href="#nutrigenomics">&rarr; details</a></h3>')
        parts.append('<ul class="key-findings">')
        for color, text in nutrition_bullets:
            parts.append(f'<li class="kf-{color}">{text}</li>')
        parts.append('</ul>')

    # ---- MENTAL HEALTH ----
    if mental_health_data and mental_health_data.get("domains"):
        mh_bullets = []
        domains = mental_health_data["domains"]
        elevated = [d for d, info in domains.items() if info["risk_level"] == "elevated"]
        low = [d for d, info in domains.items() if info["risk_level"] == "low"]
        if elevated:
            names = ", ".join(d.replace("_", " ").title() for d in elevated)
            mh_bullets.append(("yellow",
                f'Elevated genetic susceptibility for: {names}. Lifestyle and support make a huge difference.'))
        if low:
            names = ", ".join(d.replace("_", " ").title() for d in low)
            mh_bullets.append(("green", f'Low genetic susceptibility for: {names}.'))
        if mental_health_data.get("summary"):
            mh_bullets.append(("yellow" if elevated else "green",
                f'{mental_health_data["summary"]}'))

        parts.append('<h3>Mental Health <a href="#mental-health">&rarr; details</a></h3>')
        parts.append(
            '<p class="eli5-inline">Your genes influence how your brain handles '
            'stress, mood, and anxiety — but environment, relationships, and exercise '
            'matter just as much. Genes are the cards; lifestyle is how you play them.</p>'
        )
        parts.append('<ul class="key-findings">')
        for color, text in mh_bullets:
            parts.append(f'<li class="kf-{color}">{text}</li>')
        parts.append('</ul>')

    # ---- YOUR BODY ----
    body_bullets = []
    if traits_data:
        trait_summaries = []
        for tid, t in traits_data.items():
            trait_summaries.append(f'{tid.replace("_"," ").title()}: {t["prediction"]}')
        if trait_summaries:
            body_bullets.append(("green", "<strong>Traits</strong>: " + " &middot; ".join(trait_summaries)))

    if blood_type_data and blood_type_data.get("blood_type") != "Unknown":
        body_bullets.append(("green",
            f'<strong>Blood type</strong>: {blood_type_data["blood_type"]} ({blood_type_data["confidence"]} confidence)'))

    if sleep_data:
        body_bullets.append(("green",
            f'<strong>Chronotype</strong>: {sleep_data.get("chronotype", "Unknown")} '
            f'(score {sleep_data.get("chronotype_score", "?")}). '
            f'Optimal sleep: {sleep_data.get("optimal_sleep_window", "")}.'))

    if longevity_data:
        score = longevity_data.get("longevity_score", 50)
        color = "green" if score >= 60 else "yellow" if score >= 40 else "red"
        body_bullets.append((color,
            f'<strong>Longevity Score</strong>: {score}/100. {longevity_data.get("summary", "")}'))

    # Athletic profile
    for f in findings:
        if f.get("gene") == "ACTN3":
            eli5 = _eli5_for_gene("ACTN3")
            status = f.get("status", "").replace("_", " ").title()
            body_bullets.append(("green",
                f'<strong>Athletic profile (ACTN3)</strong>: {status}.'
                f'{" " + eli5 if eli5 else ""}'))
            break

    if ancestry_data and ancestry_data.get("top_ancestry"):
        body_bullets.append(("green",
            f'<strong>Ancestry</strong>: Primarily {ancestry_data["top_ancestry"]} '
            f'({ancestry_data.get("confidence", "")} confidence, '
            f'{ancestry_data.get("markers_found", 0)} markers).'))

    if body_bullets:
        parts.append('<h3>Your Body <a href="#body-profile">&rarr; details</a></h3>')
        parts.append('<ul class="key-findings">')
        for color, text in body_bullets:
            parts.append(f'<li class="kf-{color}">{text}</li>')
        parts.append('</ul>')

    # ---- GOOD NEWS ----
    good_news = (recommendations_data or {}).get("good_news", [])
    if good_news:
        parts.append('<h3>Good News</h3>')
        parts.append('<div class="good-news-grid">')
        seen_genes = set()
        for g in good_news:
            # Deduplicate by gene
            if g["gene"] in seen_genes:
                continue
            seen_genes.add(g["gene"])
            eli5 = _eli5_for_gene(g["gene"])
            # Clean ClinVar pipe text in descriptions
            desc = _clean_condition(g["description"]) if "|" in g.get("description", "") else g["description"]
            if eli5:
                desc += f' <span class="eli5-inline">({eli5})</span>'
            parts.append(
                f'<div class="good-news-card"><strong>{g["gene"]}</strong>: {desc}</div>'
            )
        parts.append('</div>')

    parts.append(
        '<div class="doctor-callout">Show this report to your doctor! '
        'They can help you use this information wisely. '
        'Each section below has the detailed data behind these findings.</div>'
    )
    return "\n".join(parts)


# =============================================================================
# SECTION 2: YOUR ACTION PLAN
# =============================================================================

def build_action_plan(recommendations_data, insights_data):
    """Build consolidated action plan: supplements, diet, lifestyle, tests, doctors."""
    if not recommendations_data:
        return "<p>No personalized recommendations available.</p>"

    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Actions grouped by urgency, with genetic rationale for each.</p>'
    )

    priorities = recommendations_data.get("priorities", [])
    clinical_insights = recommendations_data.get("clinical_insights", [])
    insights_by_gene = {ci["gene"]: ci for ci in clinical_insights}

    if not priorities:
        parts.append("<p>No convergent risk patterns detected.</p>")
    else:
        priority_colors = {"high": C["red"], "moderate": C["amber"], "low": C["green"]}

        for p in priorities:
            color = priority_colors.get(p["priority"], "var(--border)")
            eli5 = _eli5_for_condition(p.get("id", ""))

            parts.append(
                f'<details open class="rec-card" style="border-left:4px solid {color};">'
                f'<summary>'
                f'<span class="mag-badge" style="background:{color};color:#fff">'
                f'{p["priority"].upper()}</span> '
                f'<strong>{p["title"]}</strong></summary>'
            )

            if eli5:
                parts.append(f'<p class="eli5">{eli5}</p>')

            parts.append(f'<p><strong>Why:</strong> {_clean_why(p["why"])}</p>')

            parts.append("<p><strong>Actions:</strong></p><ol>")
            for action in p["actions"]:
                parts.append(f"<li>{action}</li>")
            parts.append("</ol>")

            if p.get("clinical_actions"):
                parts.append("<p><strong>Clinical actions:</strong></p><ul>")
                for ca in p["clinical_actions"]:
                    parts.append(f"<li>{ca}</li>")
                parts.append("</ul>")

            if p.get("doctor_note"):
                parts.append(
                    f'<div class="doctor-callout" style="text-align:left;font-weight:normal">'
                    f'<strong>Doctor note:</strong> {p["doctor_note"]}</div>'
                )

            if p.get("monitoring"):
                parts.append(
                    '<table><tr><th>Test</th><th>Frequency</th><th>Reason</th></tr>'
                )
                for m in p["monitoring"]:
                    parts.append(
                        f'<tr><td>{m["test"]}</td><td>{m["frequency"]}</td>'
                        f'<td>{m["reason"]}</td></tr>'
                    )
                parts.append("</table>")

            # Match clinical insights
            matched = [ci for gene, ci in insights_by_gene.items()
                       if gene.lower() in p["title"].lower()
                       or gene.lower() in p.get("why", "").lower()]
            if matched:
                parts.append(
                    '<details style="margin-top:.5em">'
                    '<summary style="font-size:.85em;color:var(--accent2)">'
                    'Clinical Context</summary>'
                )
                for ci in matched:
                    parts.append(f'<p><strong>{ci["gene"]}:</strong> {ci["mechanism"]}</p>')
                    if ci.get("actions"):
                        parts.append("<ul>")
                        for a in ci["actions"]:
                            parts.append(f"<li>{a}</li>")
                        parts.append("</ul>")
                parts.append("</details>")

            parts.append("</details>")

    # Consolidated monitoring schedule
    schedule = recommendations_data.get("monitoring_schedule", [])
    if schedule:
        parts.append("<h3>Monitoring Schedule</h3>")
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Print this and bring it to your next doctor visit.</p>'
        )
        freq_colors = {
            "weekly": C["red"], "monthly": C["orange"], "quarterly": C["amber"],
            "semi-annually": C["blue"], "annually": C["blue"], "baseline": C["green"],
        }
        parts.append('<table><tr><th>Test</th><th>Frequency</th><th>Reason</th></tr>')
        for m in schedule:
            freq = m.get("frequency", "").lower()
            color = next((c for key, c in freq_colors.items() if key in freq), "var(--accent2)")
            parts.append(
                f'<tr><td><strong>{m["test"]}</strong></td>'
                f'<td><span class="mag-badge" style="background:{color};color:#fff">'
                f'{m["frequency"]}</span></td>'
                f'<td>{m["reason"]}</td></tr>'
            )
        parts.append("</table>")

    # Specialist referrals
    referrals = recommendations_data.get("specialist_referrals", [])
    if referrals:
        parts.append("<h3>Doctors to See</h3>")
        urgency_colors = {"soon": C["red"], "routine": C["amber"]}
        for ref in referrals:
            color = urgency_colors.get(ref.get("urgency", ""), "var(--border)")
            parts.append(
                f'<div class="finding-card" style="border-left-color:{color}">'
                f'<strong>{ref["specialist"]}</strong>: {_clean_why(ref["reason"])}'
                f' <span class="badge" style="background:{color};color:#fff">'
                f'{ref.get("urgency", "routine")}</span></div>'
            )

    # Good news
    good_news = recommendations_data.get("good_news", [])
    if good_news:
        parts.append("<h3>Good News</h3>")
        parts.append('<div class="good-news-grid">')
        seen_genes = set()
        for g in good_news:
            if g["gene"] in seen_genes:
                continue
            seen_genes.add(g["gene"])
            desc = _clean_condition(g["description"]) if "|" in g.get("description", "") else g["description"]
            parts.append(
                f'<div class="good-news-card">'
                f'<strong>{g["gene"]}</strong>: {desc}</div>'
            )
        parts.append("</div>")

    return "\n".join(parts)


# =============================================================================
# SECTION 3: DRUG & MEDICATION GUIDE
# =============================================================================

def build_drug_guide(star_alleles_data, findings, pharmgkb_findings, polypharmacy_data):
    """Unified drug/medication section: star alleles + PharmGKB + polypharmacy."""
    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'How your genes affect drug processing. '
        'Show this section to every prescribing physician.</p>'
    )

    # Star allele phenotypes with gauges
    if star_alleles_data:
        parts.append("<h3>Your Drug-Processing Enzymes</h3>")
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Enzyme speed determines how fast your liver clears each drug.</p>'
        )
        parts.append('<div class="gauge-row" style="justify-content:center">')
        phenotype_levels = {
            "poor": 0, "intermediate": 0.5, "normal": 1,
            "rapid": 1.5, "ultrarapid": 2, "Unknown": 1,
        }
        phenotype_colors = {
            "poor": C["poor"], "intermediate": C["intermediate"], "normal": C["normal"],
            "rapid": C["rapid"], "ultrarapid": C["ultrarapid"], "Unknown": C["unknown_pheno"],
        }
        for gene, r in star_alleles_data.items():
            level = phenotype_levels.get(r["phenotype"], 1)
            color = phenotype_colors.get(r["phenotype"], C["unknown_pheno"])
            parts.append(svg_metabolism_gauge(gene, level, color))
        parts.append("</div>")

        parts.append(
            '<table><tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th>'
            '<th>SNPs</th><th>What This Means</th></tr>'
        )
        for gene, r in star_alleles_data.items():
            phenotype = r["phenotype"].replace("_", " ").title()
            eli5 = _eli5_for_gene(gene)
            note = eli5 if eli5 else r["clinical_note"][:120]
            parts.append(
                f'<tr><td><strong>{gene}</strong></td>'
                f'<td><code>{r["diplotype"]}</code></td>'
                f'<td>{phenotype}</td>'
                f'<td>{r["snps_found"]}/{r["snps_total"]}</td>'
                f'<td style="font-size:.85em">{note}</td></tr>'
            )
        parts.append("</table>")

    # PharmGKB annotations
    if pharmgkb_findings:
        parts.append("<h3>Drug-Gene Interactions (PharmGKB)</h3>")
        parts.append('<table><tr><th>Gene</th><th>RSID</th><th>Level</th>'
                     '<th>Drugs</th><th>Genotype</th></tr>')
        for p in pharmgkb_findings:
            parts.append(
                f'<tr><td><strong>{_esc(p["gene"])}</strong></td>'
                f'<td><code>{_esc(p["rsid"])}</code> {db_links_html(p["rsid"])}</td>'
                f'<td>{_esc(p["level"])}</td>'
                f'<td>{_esc(p["drugs"])}</td>'
                f'<td><code>{_esc(p["genotype"])}</code></td></tr>'
            )
        parts.append("</table>")

    # Polypharmacy warnings
    if polypharmacy_data and polypharmacy_data.get("warnings"):
        parts.append("<h3>Drug Combination Warnings</h3>")
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Drug combinations that may interfere based on your enzyme profile.</p>'
        )
        severity_colors = {"high": C["red"], "moderate": C["amber"], "low": C["green"]}
        for w in polypharmacy_data["warnings"]:
            color = severity_colors.get(w["severity"], "var(--border)")
            genes_str = ", ".join(f"{g}: {s}" for g, s in w.get("matched_genes", {}).items())
            drugs_str = ", ".join(w.get("drugs_affected", []))
            parts.append(
                f'<details open class="rec-card" style="border-left:4px solid {color}">'
                f'<summary><span class="mag-badge" style="background:{color};color:#fff">'
                f'{w["severity"].upper()}</span> '
                f'<strong>{_esc(w["name"])}</strong></summary>'
            )
            if genes_str:
                parts.append(f'<p><strong>Your genotype:</strong> {_esc(genes_str)}</p>')
            parts.append(f'<p><strong>Drugs affected:</strong> {_esc(drugs_str)}</p>')
            parts.append(f'<p>{_esc(w["warning"])}</p>')
            parts.append(f'<p><strong>Clinical action:</strong> {_esc(w["action"])}</p>')
            parts.append("</details>")

    # Drug card for doctor
    drug_card = (recommendations_data if isinstance(recommendations_data, dict) else {}).get("drug_card", []) if False else []
    parts.append(
        '<div class="doctor-callout">Share this entire drug section with '
        'every prescribing physician. Print it or save as PDF.</div>'
    )

    if not star_alleles_data and not pharmgkb_findings:
        return "<p>No drug-gene interaction data available.</p>"

    return "\n".join(parts)


# =============================================================================
# SECTION 4: DISEASE RISK OVERVIEW
# =============================================================================

def build_disease_risk_overview(prs_results, disease_findings_data, acmg_data):
    """Unified disease risk: PRS gauges + ClinVar pathogenic + ACMG flags."""
    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Genetic risk estimates from polygenic scores and ClinVar variants. '
        'Lifestyle and environment also affect risk significantly.</p>'
    )

    # PRS gauges
    if prs_results:
        any_non_applicable = any(not r["ancestry_applicable"] for r in prs_results.values())
        if any_non_applicable:
            parts.append(
                '<div class="doctor-callout" style="border-color:var(--warn)">'
                'PRS models are calibrated on European-ancestry populations. '
                'Your ancestry profile is substantially non-European — interpret with caution.'
                '</div>'
            )

        parts.append("<h3>Polygenic Risk Scores</h3>")
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Your combined variant score vs. the general population.</p>'
        )
        parts.append('<div class="gauge-row" style="justify-content:center">')
        for cid, r in prs_results.items():
            short_name = r["name"].replace("Age-Related ", "").replace("Macular Degeneration", "AMD")
            parts.append(svg_prs_gauge(short_name, r["percentile"], r["risk_category"]))
        parts.append("</div>")

        parts.append('<table class="sortable"><tr><th>Condition</th><th>Percentile</th>'
                     '<th>Category</th><th>SNPs</th><th>Reference</th></tr>')
        for cid, r in prs_results.items():
            parts.append(
                f'<tr><td><strong>{r["name"]}</strong></td>'
                f'<td>{r["percentile"]:.0f}th</td>'
                f'<td>{r["risk_category"].title()}</td>'
                f'<td>{r["snps_found"]}/{r["snps_total"]}</td>'
                f'<td style="font-size:.8em">{r["reference"]}</td></tr>'
            )
        parts.append("</table>")

        elevated = [r for r in prs_results.values() if r["risk_category"] in ("elevated", "high")]
        if elevated:
            parts.append("<h3>Elevated Risk Details</h3>")
            for r in elevated:
                parts.append(f'<details open><summary><strong>{r["name"]}</strong> — '
                             f'{r["percentile"]:.0f}th percentile ({r["risk_category"]})</summary>')
                if r["contributing_snps"]:
                    parts.append('<table><tr><th>Gene</th><th>rsID</th><th>Copies</th>'
                                 '<th>Effect</th></tr>')
                    for s in r["contributing_snps"][:5]:
                        parts.append(
                            f'<tr><td>{s["gene"]}</td><td><code>{s["rsid"]}</code> '
                            f'{db_links_html(s["rsid"])}</td>'
                            f'<td>{s["copies"]}</td><td>{s["contribution"]:.3f}</td></tr>'
                        )
                    parts.append("</table>")
                parts.append("</details>")

    # ACMG secondary findings
    if acmg_data:
        parts.append("<h3>Medically Actionable Genes (ACMG SF v3.2)</h3>")
        parts.append(
            f'<p style="font-size:.9em;color:var(--accent2)">'
            f'Screened {acmg_data.get("genes_screened", 81)} medically actionable genes (ACMG SF v3.2).</p>'
        )
        acmg_findings = acmg_data.get("acmg_findings", [])
        if not acmg_findings:
            parts.append(
                '<div class="doctor-callout" style="border-color:var(--green)">'
                'No pathogenic/likely pathogenic variants found in ACMG genes.</div>'
            )
        else:
            parts.append(
                f'<div class="doctor-callout" style="border-color:var(--warn)">'
                f'<strong>{len(acmg_findings)} variant(s)</strong> found in '
                f'{acmg_data["genes_with_variants"]} ACMG gene(s). '
                f'Genetic counseling recommended.</div>'
            )
            parts.append('<table><tr><th>Gene</th><th>Condition</th>'
                         '<th>Genotype</th><th>Stars</th><th>Actionability</th></tr>')
            for f in acmg_findings:
                gene = f.get("gene", "Unknown")
                condition = _clean_condition(f.get("traits") or "Unknown")
                genotype = f.get("user_genotype", "")
                stars = f.get("gold_stars", 0)
                action = f.get("acmg_actionability", "")
                eli5 = _eli5_for_gene(gene)
                parts.append(
                    f'<tr><td><strong>{gene}</strong>'
                    f'{"<br><span class=eli5-inline>" + eli5 + "</span>" if eli5 else ""}'
                    f'</td>'
                    f'<td>{condition}</td>'
                    f'<td><code>{genotype}</code></td>'
                    f'<td>{"&#9733;" * stars}{"&#9734;" * (4 - stars)}</td>'
                    f'<td style="font-size:.85em">{action}</td></tr>'
                )
            parts.append("</table>")

    # ClinVar disease findings
    if disease_findings_data:
        def _variant_table(variants, label, color):
            if not variants:
                return ""
            rows = []
            rows.append(
                f'<h3>{label} '
                f'<span class="badge" style="background:{color}">{len(variants)}</span></h3>'
            )
            rows.append(
                '<table><tr><th>Gene</th><th>Condition</th><th>Genotype</th>'
                '<th>Stars</th><th>Zygosity</th></tr>'
            )
            for v in variants[:50]:
                gene = _esc(v.get("gene", "Unknown"))
                condition = _esc(_clean_condition(v.get("traits") or "Unknown"))
                genotype = _esc(v.get("user_genotype", ""))
                stars = v.get("gold_stars", 0)
                star_str = "&#9733;" * stars + "&#9734;" * (4 - stars)
                zyg = _esc(v.get("zygosity", "").replace("_", " ").title())
                rows.append(
                    f'<tr><td><strong>{gene}</strong></td>'
                    f'<td style="max-width:400px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" '
                    f'title="{_esc(v.get("traits", ""))}">{condition}</td>'
                    f'<td><code>{genotype}</code></td>'
                    f'<td>{star_str}</td>'
                    f'<td>{zyg}</td></tr>'
                )
            rows.append("</table>")
            if len(variants) > 50:
                rows.append(f'<p style="font-size:.85em;color:var(--accent2)">'
                            f'Showing 50 of {len(variants)} variants.</p>')
            return "\n".join(rows)

        parts.append(_variant_table(
            disease_findings_data.get("pathogenic", []), "Pathogenic Variants", C["red"]))
        parts.append(_variant_table(
            disease_findings_data.get("likely_pathogenic", []), "Likely Pathogenic Variants", C["orange"]))
        parts.append(_variant_table(
            disease_findings_data.get("risk_factor", []), "Risk Factor Variants", C["amber"]))
        parts.append(_variant_table(
            disease_findings_data.get("drug_response", []), "Drug Response Variants", C["blue"]))

        protective = disease_findings_data.get("protective", [])
        if protective:
            parts.append(
                f'<h3>Protective Variants '
                f'<span class="badge" style="background:var(--green)">{len(protective)}</span></h3>'
            )
            parts.append('<div class="good-news-grid">')
            for v in protective:
                gene = _esc(v.get("gene", "Unknown"))
                condition = _esc(_clean_condition(v.get("traits") or ""))
                parts.append(
                    f'<div class="good-news-card">'
                    f'<strong>{gene}</strong>: {condition}</div>'
                )
            parts.append("</div>")

    parts.append(
        '<p style="font-size:.85em;color:var(--accent2)">'
        'PRS estimates relative genetic risk from common variants. '
        'Lifestyle, environment, and rare variants also affect risk. '
        'Not a clinical diagnosis.</p>'
    )

    if not prs_results and not disease_findings_data and not acmg_data:
        return "<p>No disease risk data available.</p>"

    return "\n".join(parts)


# =============================================================================
# SECTION 5: YOUR BODY PROFILE
# =============================================================================

def build_body_profile(traits_data, blood_type_data, sleep_data, longevity_data, findings):
    """Fun traits, blood type, chronotype, longevity, athletic profile."""
    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Traits, blood type, chronotype, longevity, and athletic profile.</p>'
    )

    # Traits
    if traits_data:
        parts.append("<h3>Predicted Traits</h3>")
        trait_labels = {
            "eye_color": "Eye Color", "hair_color": "Hair Color",
            "earwax_type": "Earwax Type", "freckling": "Freckling / Sun Sensitivity",
        }
        parts.append('<table><tr><th>Trait</th><th>Prediction</th>'
                     '<th>Confidence</th><th>What It Means</th></tr>')
        for trait_id, trait in traits_data.items():
            label = trait_labels.get(trait_id, trait_id.replace("_", " ").title())
            conf_color = {"high": "var(--green)", "moderate": "var(--accent)",
                          "low": "var(--warn)"}.get(trait["confidence"], "inherit")
            parts.append(
                f'<tr><td><strong>{label}</strong></td>'
                f'<td>{trait["prediction"]}</td>'
                f'<td style="color:{conf_color}">{trait["confidence"].title()}</td>'
                f'<td style="font-size:.85em">{trait["description"]}</td></tr>'
            )
        parts.append("</table>")

    # Blood type
    if blood_type_data and blood_type_data.get("blood_type") != "Unknown":
        parts.append("<h3>Blood Type</h3>")
        bt = blood_type_data["blood_type"]
        conf = blood_type_data["confidence"].title()
        parts.append(
            f'<div style="text-align:center;margin:1em 0">'
            f'<span style="font-size:3em;font-weight:bold;color:var(--warn)">{bt}</span>'
            f'<br><span style="font-size:.9em;color:var(--accent2)">'
            f'Confidence: {conf}</span></div>'
        )
        parts.append(
            '<p style="font-size:.85em;color:var(--accent2)">Blood type from SNP data has limitations. '
            'Confirm with clinical blood typing.</p>'
        )

    # Sleep / Chronotype
    if sleep_data:
        parts.append("<h3>Sleep &amp; Chronotype</h3>")
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Genetic influence on circadian rhythm and sleep architecture.</p>'
        )
        parts.append(
            f'<div style="text-align:center;margin:1em 0">'
            f'<p style="font-size:1.5em;font-weight:bold">{_esc(sleep_data.get("chronotype", "Unknown"))}</p>'
            f'<p style="color:var(--accent2)">Chronotype Score: {sleep_data.get("chronotype_score", 50)}/100 '
            f'(0=extreme morning, 100=extreme evening)</p>'
            f'</div>'
        )
        parts.append(
            f'<table>'
            f'<tr><td><strong>Optimal Sleep Window</strong></td><td>{_esc(sleep_data.get("optimal_sleep_window", ""))}</td></tr>'
            f'<tr><td><strong>Peak Alertness</strong></td><td>{_esc(sleep_data.get("peak_alertness", ""))}</td></tr>'
            f'<tr><td><strong>Last Caffeine By</strong></td><td>{_esc(sleep_data.get("caffeine_cutoff", ""))}</td></tr>'
            f'</table>'
        )
        recs = sleep_data.get("recommendations", [])
        if recs:
            parts.append("<p><strong>Sleep Recommendations:</strong></p><ul>")
            for r in recs:
                parts.append(f"<li>{_esc(r)}</li>")
            parts.append("</ul>")

    # Longevity
    if longevity_data:
        parts.append("<h3>Longevity &amp; Healthspan</h3>")
        score = longevity_data.get("longevity_score", 50)
        score_color = C["green"] if score >= 60 else C["amber"] if score >= 40 else C["red"]
        parts.append(
            f'<div style="text-align:center;margin:1em 0">'
            f'<span style="font-size:2.5em;font-weight:bold;color:{score_color}">{score}</span>'
            f'<span style="font-size:1.2em;color:var(--accent2)">/100</span>'
            f'<p style="color:var(--accent2)">Longevity Genetic Score</p>'
            f'<p style="font-size:.9em">{_esc(longevity_data.get("summary", ""))}</p>'
            f'</div>'
        )
        domains = longevity_data.get("healthspan_domains", {})
        if domains:
            parts.append('<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:1em">')
            for domain, info in domains.items():
                s = info["score"]
                c = C["green"] if s >= 60 else C["amber"] if s >= 40 else C["red"]
                parts.append(
                    f'<div style="border:1px solid var(--border);border-radius:8px;padding:1em;text-align:center">'
                    f'<strong>{_esc(domain.replace("_", " ").title())}</strong><br>'
                    f'<span style="font-size:1.5em;color:{c}">{s}</span>/90<br>'
                    f'<span style="font-size:.85em;color:var(--accent2)">{info["rating"]}</span>'
                    f'</div>'
                )
            parts.append("</div>")
        interventions = longevity_data.get("interventions", [])
        if interventions:
            parts.append("<p><strong>Genetically-Supported Interventions:</strong></p><ul>")
            for i in interventions:
                parts.append(f'<li><strong>{_esc(i["intervention"])}</strong> — {_esc(i["why"])}</li>')
            parts.append("</ul>")

    # Athletic profile from ACTN3
    for f in findings:
        if f.get("gene") == "ACTN3":
            parts.append("<h3>Athletic Profile</h3>")
            eli5 = _eli5_for_gene("ACTN3")
            parts.append(f'<p class="eli5">{eli5}</p>')
            status = f.get("status", "").replace("_", " ").title()
            parts.append(f'<p><strong>ACTN3</strong>: {status} — {f.get("description", "")}</p>')
            parts.append(paper_refs_html("rs1815739"))
            break

    return "\n".join(parts)


# =============================================================================
# SECTION 6: MENTAL HEALTH & BRAIN
# =============================================================================

def build_mental_health_section(data):
    """Mental health genetics: risk domains, resilience, treatment matching."""
    if not data:
        return "<p>No mental health genetic data available.</p>"
    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Genetic susceptibility markers. Environment and lifestyle are equally important.</p>'
    )
    parts.append(f'<p><strong>{_esc(data.get("summary", ""))}</strong></p>')

    domains = data.get("domains", {})
    if domains:
        parts.append('<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:1em;margin:1em 0">')
        level_colors = {"elevated": C["red"], "moderate": C["amber"], "low": C["green"]}
        for domain, info in domains.items():
            color = level_colors.get(info["risk_level"], "var(--border)")
            parts.append(
                f'<div style="border:1px solid var(--border);border-radius:8px;padding:1em;text-align:center">'
                f'<strong>{_esc(domain.replace("_", " ").title())}</strong><br>'
                f'<span style="font-size:1.3em;color:{color}">{_esc(info["risk_level"].title())}</span>'
                f'</div>'
            )
        parts.append("</div>")

    risk = data.get("risk_factors", [])
    resilience = data.get("resilience_factors", [])
    if risk:
        parts.append("<h3>Risk Factors</h3><ul>")
        for r in risk:
            parts.append(f"<li>{_esc(r)}</li>")
        parts.append("</ul>")
    if resilience:
        parts.append('<h3>Resilience Factors</h3><ul>')
        for r in resilience:
            parts.append(f"<li>{_esc(r)}</li>")
        parts.append("</ul>")

    notes = data.get("treatment_notes", [])
    if notes:
        parts.append("<h3>Treatment Matching Notes</h3>")
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Share these with your therapist or psychiatrist.</p>'
        )
        parts.append("<ul>")
        for n in notes:
            parts.append(f"<li>{_esc(n)}</li>")
        parts.append("</ul>")

    recs = data.get("recommendations", [])
    if recs:
        parts.append("<h3>Recommendations</h3><ol>")
        for r in recs:
            parts.append(f"<li>{_esc(r)}</li>")
        parts.append("</ol>")
    return "\n".join(parts)


# =============================================================================
# SECTION 7: CLINICAL FINDINGS DETAIL
# =============================================================================

def build_clinical_detail(findings, epistasis_results, carrier_screen_data):
    """Full technical findings table + epistasis + carrier screening."""
    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Complete technical findings for clinical review. '
        'Searchable, sortable, and exportable.</p>'
    )

    # All lifestyle findings grouped by category
    categories = sorted(set(f.get("category", "Other") for f in findings))
    for cat in categories:
        cat_findings = sorted(
            [f for f in findings if f.get("category") == cat],
            key=lambda x: -x.get("magnitude", 0),
        )
        if not cat_findings:
            continue

        parts.append(
            f'<details class="category-section" open>'
            f'<summary><h3 class="inline-h3">{cat} '
            f'<span class="badge">{len(cat_findings)}</span></h3></summary>'
        )

        for f in cat_findings:
            mag = f.get("magnitude", 0)
            mag_class = "high" if mag >= 3 else "mod" if mag == 2 else "low" if mag == 1 else "info"
            rsid = _esc(f.get("rsid", ""))
            gene = _esc(f.get("gene", "Unknown"))
            genotype = _esc(f.get("genotype", ""))
            status = _esc(f.get("status", "").replace("_", " ").title())
            desc = _esc(f.get("description", ""))
            note = _esc(f.get("note", ""))

            parts.append(f'<div class="finding-card {mag_class}">')
            parts.append(
                f'<div class="finding-header">'
                f'<span class="mag-badge mag-{mag_class}">{mag}/6</span> '
                f'<strong>{gene}</strong> '
                f'<code>{rsid}</code> — '
                f'<code>{genotype}</code> — {status}'
                f'</div>'
            )
            eli5 = _eli5_for_gene(f.get("gene", ""))
            if eli5:
                parts.append(f'<p class="eli5-inline" style="font-size:.85em;margin:.2em 0">{eli5}</p>')
            parts.append(f'<p class="finding-desc">{desc}</p>')
            if note:
                parts.append(f'<p class="finding-note">Note: {note}</p>')
            freq = f.get("freq")
            if freq and isinstance(freq, dict):
                freq_labels = {"EUR": "European", "AFR": "African",
                               "EAS": "East Asian", "SAS": "South Asian", "AMR": "American"}
                freq_parts = [f'{freq_labels.get(p, p)}: {v:.0%}'
                              for p, v in sorted(freq.items(), key=lambda x: -x[1])
                              if v > 0.001]
                if freq_parts:
                    parts.append(
                        f'<p class="finding-note" style="font-size:.8em">'
                        f'Population freq: {" &middot; ".join(freq_parts)}</p>'
                    )
            parts.append(db_links_html(rsid))
            parts.append(paper_refs_html(rsid))
            parts.append("</div>")

        parts.append("</details>")

    # Pathway analysis
    from ..clinical_context import PATHWAYS
    parts.append("<h3>Pathway Analysis</h3>")
    gene_map = {f["gene"]: f for f in findings}
    for pathway_name, pathway_genes in PATHWAYS.items():
        pathway_findings = [gene_map[g] for g in pathway_genes if g in gene_map]
        if not pathway_findings:
            continue
        parts.append(f'<details><summary><strong>{_esc(pathway_name)}</strong></summary><ul>')
        for pf in pathway_findings:
            mag = pf.get("magnitude", 0)
            mag_class = "high" if mag >= 3 else "mod" if mag == 2 else "low" if mag == 1 else "info"
            parts.append(
                f'<li><span class="mag-dot mag-{mag_class}"></span> '
                f'<strong>{_esc(pf["gene"])}</strong>: '
                f'{_esc(pf.get("status","").replace("_"," ").title())}</li>'
            )
        parts.append("</ul></details>")

    # Epistasis
    if epistasis_results:
        parts.append("<h3>Gene-Gene Interactions (Epistasis)</h3>")
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Combined effects when multiple gene variants interact.</p>'
        )
        risk_colors = {"high": "var(--warn)", "moderate": "var(--accent)", "low": "var(--green)"}
        for interaction in epistasis_results:
            color = risk_colors.get(interaction["risk_level"], "var(--accent)")
            genes_involved = interaction.get("genes_involved", {})
            genes = ", ".join(genes_involved.keys()) if genes_involved else _esc(interaction.get("name", ""))
            parts.append(
                f'<details open><summary>'
                f'<span class="mag-badge" style="background:{color};color:#fff">'
                f'{interaction["risk_level"].upper()}</span> '
                f'<strong>{interaction["name"]}</strong></summary>'
            )
            parts.append(f'<p><strong>Genes:</strong> {genes}</p>')
            parts.append(f'<p><strong>Effect:</strong> {interaction["effect"]}</p>')
            parts.append(f'<p><strong>Mechanism:</strong> {interaction["mechanism"]}</p>')
            parts.append("<p><strong>Recommended Actions:</strong></p><ul>")
            for action in interaction["actions"]:
                parts.append(f"<li>{action}</li>")
            parts.append("</ul></details>")

    # Carrier screening
    if carrier_screen_data and carrier_screen_data.get("total_carriers", 0) > 0:
        parts.append("<h3>Carrier Screening</h3>")
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Carrier status: you have one copy of a recessive variant. '
            'Relevant for reproductive planning.</p>'
        )
        parts.append(
            f'<p>{carrier_screen_data["total_carriers"]} carrier finding(s) organized by disease system.</p>'
        )
        for system, carriers in sorted(carrier_screen_data.get("by_system", {}).items()):
            parts.append(f'<h4>{system}</h4><ul>')
            for c in carriers:
                note = f' &mdash; <em>{c["reproductive_note"]}</em>' if c.get("reproductive_note") else ""
                parts.append(
                    f'<li><strong>{c["gene"]}</strong> ({c.get("rsid", "")}): '
                    f'{c["condition"]} <span class="badge">{c["inheritance"]}</span>{note}</li>'
                )
            parts.append("</ul>")

        couples = carrier_screen_data.get("couples_relevant", [])
        if couples:
            parts.append('<h4>Couples-Relevant Conditions</h4><ul>')
            for c in couples:
                parts.append(f'<li><strong>{c["gene"]}</strong>: {c["condition"]}</li>')
            parts.append("</ul>")

    return "\n".join(parts)


# =============================================================================
# SECTION 8: ANCESTRY & POPULATION CONTEXT
# =============================================================================

def build_ancestry_section(ancestry_results, mt_haplogroup_data=None):
    """Ancestry proportions + MT haplogroup + population warnings."""
    if not ancestry_results or ancestry_results.get("markers_found", 0) == 0:
        if mt_haplogroup_data and mt_haplogroup_data.get("haplogroup") != "Unknown":
            pass  # still show MT haplogroup
        else:
            return "<p>No ancestry-informative markers found in genome data.</p>"

    parts = []

    if ancestry_results and ancestry_results.get("markers_found", 0) > 0:
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Superpopulation estimates from ~55 ancestry-informative markers.</p>'
        )
        parts.append('<div class="chart-grid">')
        parts.append("<div>")
        parts.append("<h3>Ancestry Proportions</h3>")
        parts.append(svg_ancestry_donut(ancestry_results))
        parts.append("</div>")

        parts.append("<div>")
        parts.append("<h3>Details</h3>")
        parts.append(f'<p><strong>Confidence:</strong> {ancestry_results["confidence"].title()} '
                     f'({ancestry_results["markers_found"]} markers)</p>')
        parts.append('<table><tr><th>Population</th><th>Proportion</th></tr>')
        labels = {"EUR": "European", "AFR": "African", "EAS": "East Asian",
                  "SAS": "South Asian", "AMR": "Admixed American"}
        for pop in sorted(ancestry_results["proportions"],
                          key=lambda p: -ancestry_results["proportions"][p]):
            prop = ancestry_results["proportions"][pop]
            label = labels.get(pop, pop)
            parts.append(f'<tr><td>{label} ({pop})</td><td>{prop:.1%}</td></tr>')
        parts.append("</table>")
        parts.append("</div>")
        parts.append("</div>")

    # MT Haplogroup
    if mt_haplogroup_data and mt_haplogroup_data.get("haplogroup") != "Unknown":
        mt = mt_haplogroup_data
        parts.append("<h3>Maternal Haplogroup (Mitochondrial DNA)</h3>")
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Maternal lineage from mitochondrial DNA.</p>'
        )
        parts.append(
            f'<div style="text-align:center;margin:1em 0">'
            f'<span style="font-size:2.5em;font-weight:bold;color:var(--accent)">{mt["haplogroup"]}</span>'
            f'<br><span style="font-size:.95em">{mt["description"]}</span>'
            f'<br><span style="font-size:.85em;color:var(--accent2)">'
            f'{mt["lineage"]} — {mt["confidence"].title()} confidence '
            f'({mt["markers_found"]}/{mt["markers_tested"]} markers)</span></div>'
        )

    return "\n".join(parts)


# =============================================================================
# SECTION 9: NUTRIGENOMICS DETAIL
# =============================================================================

def build_nutrigenomics_section(nutrigenomics_data, recommendations_data, insights_data):
    """Full nutrient-by-nutrient breakdown + lab tests."""
    parts = []

    # Nutrigenomics data
    if nutrigenomics_data:
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'How your genes affect nutrient needs. Personalized supplement and dietary guidance.</p>'
        )
        parts.append(f'<p>{_esc(nutrigenomics_data.get("summary", ""))}</p>')

        needs = nutrigenomics_data.get("nutrient_needs", [])
        if needs:
            need_colors = {"high": C["red"], "moderate": C["amber"], "low": C["blue"],
                           "normal": C["green"], "caution_excess": C["purple"]}
            for n in needs:
                if n["need_level"] == "normal" and not n["gene_impacts"]:
                    continue
                color = need_colors.get(n["need_level"], "var(--border)")
                parts.append(
                    f'<details class="rec-card" style="border-left:4px solid {color}">'
                    f'<summary><span class="mag-badge" style="background:{color};color:#fff">'
                    f'{_esc(n["need_level"].upper().replace("_", " "))}</span> '
                    f'<strong>{_esc(n["name"])}</strong></summary>'
                )
                parts.append(f'<p>{_esc(n["recommendation"])}</p>')
                parts.append(f'<p><strong>Food sources:</strong> {_esc(n["food_sources"])}</p>')
                for gi in n.get("gene_impacts", []):
                    parts.append(f'<p style="font-size:.9em">— {_esc(gi["gene"])}: {_esc(gi["impact"])}</p>')
                parts.append("</details>")

        supps = nutrigenomics_data.get("supplement_priorities", [])
        if supps:
            parts.append("<h3>Supplement Priority List</h3>")
            parts.append('<table><tr><th>Nutrient</th><th>Form</th><th>Dose</th><th>Priority</th></tr>')
            for s in supps:
                parts.append(
                    f'<tr><td><strong>{_esc(s["nutrient"])}</strong></td>'
                    f'<td>{_esc(s["form"])}</td>'
                    f'<td>{_esc(s["dose"])}</td>'
                    f'<td>{_esc(s["priority"])}</td></tr>'
                )
            parts.append("</table>")

        testing = nutrigenomics_data.get("testing_recommendations", [])
        if testing:
            parts.append("<h3>Recommended Lab Tests</h3><ul>")
            for t in testing:
                parts.append(f'<li><strong>{_esc(t["nutrient"])}:</strong> {_esc(t["test"])}</li>')
            parts.append("</ul>")

    # Nutrition narratives from insights
    narratives = (insights_data or {}).get("narratives", [])
    nutrition_narrative_ids = {"methylation_choline", "vitamin_d_profile", "iron_profile", "caffeine_metabolism"}
    nutrition_narratives = [n for n in narratives if n.get("id", "") in nutrition_narrative_ids]
    if nutrition_narratives:
        parts.append("<h3>Nutritional Gene Stories</h3>")
        for n in nutrition_narratives:
            genes_str = ", ".join(n["matched_genes"])
            parts.append(
                f'<details open class="rec-card" style="border-left:4px solid var(--accent2)">'
                f'<summary><strong>{n["title"]}</strong> '
                f'<span class="badge" style="background:var(--accent2)">'
                f'{len(n["matched_genes"])} genes</span></summary>'
                f'<p><strong>Genes:</strong> {genes_str}</p>'
                f'<p>{n["narrative"]}</p>'
                f'<p><strong>What this means for you:</strong> {n["practical"]}</p>'
            )
            if n.get("references"):
                parts.append('<div class="paper-refs">References: ')
                ref_links = []
                for ref in n["references"]:
                    if "PMID:" in ref:
                        pmid = ref.split("PMID:")[1].strip().rstrip(",.")
                        ref_links.append(
                            f'<a href="https://pubmed.ncbi.nlm.nih.gov/{pmid}/" '
                            f'target="_blank" rel="noopener">{ref}</a>'
                        )
                    else:
                        ref_links.append(ref)
                parts.append(" | ".join(ref_links) + "</div>")
            parts.append("</details>")

    if not parts:
        return "<p>No nutrigenomics data available.</p>"

    return "\n".join(parts)


# =============================================================================
# SECTION 10: DATA QUALITY
# =============================================================================

def build_quality_section(metrics):
    """Call rate, het rate, coverage."""
    if not metrics:
        return "<p>No quality metrics available.</p>"

    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Raw data quality assessment — call rate, heterozygosity, and coverage.</p>'
    )

    call_pct = metrics.get("call_rate", 0) * 100
    if call_pct >= 99:
        badge_color = "var(--green)"
        badge_text = "Excellent"
    elif call_pct >= 97:
        badge_color = C["blue"]
        badge_text = "Good"
    elif call_pct >= 95:
        badge_color = C["amber"]
        badge_text = "Fair"
    else:
        badge_color = "var(--warn)"
        badge_text = "Low"

    parts.append(
        f'<p><strong>Call Rate:</strong> '
        f'<span class="mag-badge" style="background:{badge_color};color:#fff">'
        f'{call_pct:.1f}% — {badge_text}</span></p>'
    )

    parts.append(
        f'<table><tr><th>Metric</th><th>Value</th></tr>'
        f'<tr><td>Total SNPs</td><td>{metrics["total_snps"]:,}</td></tr>'
        f'<tr><td>No-call positions</td><td>{metrics["no_call_count"]:,}</td></tr>'
        f'<tr><td>Autosomal SNPs</td><td>{metrics["autosomal_count"]:,}</td></tr>'
        f'<tr><td>Mitochondrial SNPs</td><td>{metrics["mt_snp_count"]}</td></tr>'
        f'<tr><td>Heterozygosity rate</td><td>{metrics["het_rate"]:.3f}</td></tr>'
        f'<tr><td>Inferred sex</td><td>{"Male (Y detected)" if metrics["has_y"] else "Female"}</td></tr>'
        f'</table>'
    )

    chroms = metrics.get("chromosomes", {})
    autosomal = sorted(
        [(ch, cnt) for ch, cnt in chroms.items() if ch.isdigit()],
        key=lambda x: int(x[0])
    )
    if autosomal:
        max_cnt = max(cnt for _, cnt in autosomal)
        bar_w = 300
        bar_h = 14
        height = len(autosomal) * (bar_h + 4) + 20
        bars = []
        for i, (ch, cnt) in enumerate(autosomal):
            y = 10 + i * (bar_h + 4)
            w = int(cnt / max_cnt * bar_w) if max_cnt else 0
            bars.append(
                f'<text x="30" y="{y + 11}" text-anchor="end" '
                f'fill="currentColor" font-size="10">{ch}</text>'
                f'<rect x="35" y="{y}" width="{w}" height="{bar_h}" '
                f'rx="2" fill="var(--accent)" opacity="0.7"/>'
                f'<text x="{40 + w}" y="{y + 11}" '
                f'fill="currentColor" font-size="9">{cnt:,}</text>'
            )
        parts.append(
            f'<h3>Chromosome Coverage</h3>'
            f'<svg viewBox="0 0 440 {height}" class="chart" role="img" '
            f'aria-label="Chromosome SNP coverage">'
            f'<title>SNPs per Chromosome</title>'
            f'{"".join(bars)}</svg>'
        )

    return "\n".join(parts)


# =============================================================================
# SECTION 11: DOCTOR CARD
# =============================================================================

def build_doctor_card(recommendations_data, star_alleles_data, apoe_data, acmg_data):
    """Print-optimized doctor card with key clinical findings."""
    parts = []
    parts.append('<div style="border:2px solid var(--accent);border-radius:8px;padding:1.5em">')
    parts.append('<h3 style="margin-top:0;color:var(--accent)">Patient Genetic Summary</h3>')

    priorities = (recommendations_data or {}).get("priorities", [])
    high_priorities = [p for p in priorities if p["priority"] == "high"]
    if high_priorities:
        parts.append("<h4>High-Priority Conditions</h4>")
        parts.append('<table><tr><th>Condition</th><th>Doctor Note</th></tr>')
        for p in high_priorities:
            note = p.get("doctor_note", "") or p.get("why", "")
            note = _clean_why(note)
            parts.append(
                f'<tr><td><strong>{p["title"]}</strong></td>'
                f'<td>{note}</td></tr>'
            )
        parts.append("</table>")

    referrals = (recommendations_data or {}).get("specialist_referrals", [])
    if referrals:
        parts.append("<h4>Specialist Referrals</h4>")
        urgency_colors = {"soon": C["red"], "routine": C["amber"]}
        for ref in referrals:
            color = urgency_colors.get(ref.get("urgency", ""), "var(--border)")
            parts.append(
                f'<p><strong>{ref["specialist"]}</strong>: {_clean_why(ref["reason"])} '
                f'<span class="mag-badge" style="background:{color};color:#fff">'
                f'{ref.get("urgency", "routine")}</span></p>'
            )

    if star_alleles_data:
        parts.append("<h4>Pharmacogenomic Profile</h4>")
        parts.append('<table><tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th></tr>')
        for gene, r in star_alleles_data.items():
            phenotype = r["phenotype"].replace("_", " ").title()
            parts.append(
                f'<tr><td><strong>{gene}</strong></td>'
                f'<td><code>{r["diplotype"]}</code></td>'
                f'<td>{phenotype}</td></tr>'
            )
        parts.append("</table>")

    if apoe_data and apoe_data.get("apoe_type") != "Unknown":
        risk_colors = {
            "reduced": "var(--green)", "average": C["blue"],
            "moderate": C["amber"], "elevated": C["orange"], "high": C["red"],
        }
        color = risk_colors.get(apoe_data.get("risk_level", ""), "var(--accent2)")
        parts.append(
            f'<h4>APOE Status</h4>'
            f'<p><strong>{apoe_data["apoe_type"]}</strong> — '
            f'<span class="mag-badge" style="background:{color};color:#fff">'
            f'{apoe_data["risk_level"].title()} Risk</span></p>'
        )

    acmg_findings = (acmg_data or {}).get("acmg_findings", [])
    if acmg_findings:
        parts.append("<h4>ACMG Medically Actionable Findings</h4>")
        parts.append('<table><tr><th>Gene</th><th>Condition</th><th>Actionability</th></tr>')
        for f in acmg_findings:
            gene = f.get("gene", "Unknown")
            condition = _clean_condition(f.get("traits") or "Unknown")
            action = f.get("acmg_actionability", "")
            parts.append(
                f'<tr><td><strong>{gene}</strong></td>'
                f'<td>{condition}</td>'
                f'<td>{action}</td></tr>'
            )
        parts.append("</table>")
        parts.append(
            '<div class="doctor-callout" style="border-color:var(--warn)">'
            'ACMG actionable findings detected. Refer to genetic counselor.</div>'
        )

    schedule = (recommendations_data or {}).get("monitoring_schedule", [])
    if schedule:
        parts.append("<h4>Recommended Monitoring</h4>")
        parts.append('<table><tr><th>Test</th><th>Frequency</th></tr>')
        for m in schedule:
            parts.append(
                f'<tr><td>{m["test"]}</td><td>{m["frequency"]}</td></tr>'
            )
        parts.append("</table>")

    parts.append("</div>")

    if not any([high_priorities, referrals, star_alleles_data, acmg_findings]):
        return "<p>No significant clinical findings for doctor review.</p>"

    return "\n".join(parts)


# =============================================================================
# SECTION 12: REFERENCES & DATABASE LINKS
# =============================================================================

def build_references(findings):
    """All rsID links, paper citations, methodology."""
    parts = []

    # Curated papers
    parts.append("<h3>Key Papers</h3><ul>")
    seen = set()
    for rsid, refs in PAPER_REFS.items():
        for r in refs:
            key = r["pmid"]
            if key in seen:
                continue
            seen.add(key)
            parts.append(
                f'<li><a href="https://pubmed.ncbi.nlm.nih.gov/{r["pmid"]}/" '
                f'target="_blank" rel="noopener">{r["title"]}</a> '
                f'(PMID: {r["pmid"]}, {r["year"]})</li>'
            )
    parts.append("</ul>")

    # All rsID links
    rsids = sorted(set(f.get("rsid", "") for f in findings if f.get("rsid", "").startswith("rs")))
    if rsids:
        parts.append("<h3>Database Links for All Analyzed rsIDs</h3>")
        parts.append('<div class="rsid-grid">')
        for rsid in rsids:
            ref_html = paper_refs_html(rsid)
            parts.append(f"<div><code>{rsid}</code> {db_links_html(rsid)}{ref_html}</div>")
        parts.append("</div>")

    # Methodology
    parts.append("<h3>Methodology &amp; Disclaimers</h3>")
    parts.append(
        "<ul>"
        "<li>Lifestyle findings from curated SNP database (~260 variants)</li>"
        "<li>Disease risk from ClinVar (~341K variants scanned)</li>"
        "<li>Drug interactions from PharmGKB clinical annotations</li>"
        "<li>Polygenic risk scores from published GWAS (8 conditions)</li>"
        "<li>Star allele calling for 6 pharmacogenes (CPIC-style)</li>"
        "<li>Only true SNPs analyzed (indels filtered to prevent false positives)</li>"
        "<li>This report is for <strong>informational purposes only</strong> — not a clinical diagnosis</li>"
        "<li>Genetic associations are probabilistic, not deterministic</li>"
        "<li>Always consult healthcare providers before making medical decisions</li>"
        "</ul>"
    )

    return "\n".join(parts)


# =============================================================================
# BACKWARD-COMPAT WRAPPERS (used by tests)
# =============================================================================

def build_prs_section(prs_results):
    """Build HTML for polygenic risk scores (test-compatible wrapper)."""
    if not prs_results:
        return "<p>No PRS data available.</p>"
    return build_disease_risk_overview(prs_results, None, None)


def build_epistasis_section(epistasis_results):
    """Build HTML for gene-gene interactions (test-compatible wrapper)."""
    if not epistasis_results:
        return "<p>No significant gene-gene interactions detected.</p>"
    parts = []
    risk_colors = {"high": "var(--warn)", "moderate": "var(--accent)", "low": "var(--green)"}
    for interaction in epistasis_results:
        color = risk_colors.get(interaction["risk_level"], "var(--accent)")
        genes_involved = interaction.get("genes_involved", {})
        genes = ", ".join(genes_involved.keys()) if genes_involved else _esc(interaction.get("name", ""))
        parts.append(
            f'<details open><summary>'
            f'<span class="mag-badge" style="background:{color};color:#fff">'
            f'{interaction["risk_level"].upper()}</span> '
            f'<strong>{interaction["name"]}</strong></summary>'
        )
        parts.append(f'<p><strong>Genes:</strong> {genes}</p>')
        parts.append(f'<p><strong>Effect:</strong> {interaction["effect"]}</p>')
        parts.append(f'<p><strong>Mechanism:</strong> {interaction["mechanism"]}</p>')
        parts.append("<p><strong>Recommended Actions:</strong></p><ul>")
        for action in interaction["actions"]:
            parts.append(f"<li>{action}</li>")
        parts.append("</ul></details>")
    return "\n".join(parts)


def build_disease_risk(disease_findings):
    """Build disease risk section from ClinVar findings (test-compatible wrapper)."""
    if not disease_findings:
        return "<p>No ClinVar disease risk data available.</p>"

    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Variants identified by scanning your genome against the ClinVar database.</p>'
    )

    def _variant_table(variants, label, color):
        if not variants:
            return ""
        rows = []
        rows.append(
            f'<h3>{label} '
            f'<span class="badge" style="background:{color}">{len(variants)}</span></h3>'
        )
        rows.append(
            '<table><tr><th>Gene</th><th>Condition</th><th>Genotype</th>'
            '<th>Stars</th><th>Zygosity</th></tr>'
        )
        for v in variants[:50]:
            gene = _esc(v.get("gene", "Unknown"))
            condition = _esc(_clean_condition(v.get("traits") or "Unknown"))
            genotype = _esc(v.get("user_genotype", ""))
            stars = v.get("gold_stars", 0)
            star_str = "&#9733;" * stars + "&#9734;" * (4 - stars)
            zyg = _esc(v.get("zygosity", "").replace("_", " ").title())
            rows.append(
                f'<tr><td><strong>{gene}</strong></td>'
                f'<td style="max-width:400px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap" title="{_esc(v.get("traits", ""))}">{condition}</td>'
                f'<td><code>{genotype}</code></td>'
                f'<td>{star_str}</td>'
                f'<td>{zyg}</td></tr>'
            )
        rows.append("</table>")
        return "\n".join(rows)

    parts.append(_variant_table(
        disease_findings.get("pathogenic", []), "Pathogenic Variants", C["red"]))
    parts.append(_variant_table(
        disease_findings.get("likely_pathogenic", []), "Likely Pathogenic Variants", C["orange"]))
    parts.append(_variant_table(
        disease_findings.get("risk_factor", []), "Risk Factor Variants", C["amber"]))
    parts.append(_variant_table(
        disease_findings.get("drug_response", []), "Drug Response Variants", C["blue"]))

    protective = disease_findings.get("protective", [])
    if protective:
        parts.append(
            f'<h3>Protective Variants '
            f'<span class="badge" style="background:var(--green)">{len(protective)}</span></h3>'
        )
        parts.append('<div class="good-news-grid">')
        for v in protective:
            gene = _esc(v.get("gene", "Unknown"))
            condition = _esc(_clean_condition(v.get("traits") or ""))
            parts.append(
                f'<div class="good-news-card">'
                f'<strong>{gene}</strong>: {condition}</div>'
            )
        parts.append("</div>")

    return "\n".join(parts)


def build_monitoring(recommendations_data):
    """Build monitoring schedule (test-compatible wrapper)."""
    schedule = (recommendations_data or {}).get("monitoring_schedule", [])
    if not schedule:
        return "<p>No monitoring tests recommended based on your genetic profile.</p>"
    parts = []
    freq_colors = {
        "weekly": C["red"], "monthly": C["orange"], "quarterly": C["amber"],
        "semi-annually": C["blue"], "annually": C["blue"], "baseline": C["green"],
    }
    parts.append('<table><tr><th>Test</th><th>Frequency</th><th>Reason</th></tr>')
    for m in schedule:
        freq = m.get("frequency", "").lower()
        color = next((c for key, c in freq_colors.items() if key in freq), "var(--accent2)")
        parts.append(
            f'<tr><td><strong>{m["test"]}</strong></td>'
            f'<td><span class="mag-badge" style="background:{color};color:#fff">'
            f'{m["frequency"]}</span></td>'
            f'<td>{m["reason"]}</td></tr>'
        )
    parts.append("</table>")
    parts.append(
        '<div class="doctor-callout" style="border-color:var(--accent)">'
        'Print this monitoring schedule and bring it to your next doctor visit.</div>'
    )
    return "\n".join(parts)


def build_nutrition_section(recommendations_data, insights_data):
    """Build nutrition section (test-compatible wrapper)."""
    if not recommendations_data and not insights_data:
        return "<p>No nutrition data available.</p>"
    parts = []
    nutrition_groups = {"methylation", "iron_metabolism", "caffeine", "metabolic_diabetes",
                        "nutrition", "vitamin_d", "iron"}
    nutrition_narrative_ids = {"methylation_choline", "vitamin_d_profile", "iron_profile",
                               "caffeine_metabolism"}
    priorities = (recommendations_data or {}).get("priorities", [])
    nutrition_priorities = [
        p for p in priorities
        if any(g in p.get("id", "").lower() for g in nutrition_groups)
    ]
    if nutrition_priorities:
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Supplement and dietary guidance based on your genetic variants.</p>'
        )
        priority_colors = {"high": C["red"], "moderate": C["amber"], "low": C["green"]}
        for p in nutrition_priorities:
            color = priority_colors.get(p["priority"], "var(--border)")
            parts.append(
                f'<details open class="rec-card" style="border-left:4px solid {color}">'
                f'<summary>'
                f'<span class="mag-badge" style="background:{color};color:#fff">'
                f'{p["priority"].upper()}</span> '
                f'<strong>{p["title"]}</strong></summary>'
            )
            parts.append(f'<p><strong>Why:</strong> {_clean_why(p["why"])}</p>')
            parts.append("<p><strong>Actions:</strong></p><ol>")
            for action in p["actions"]:
                parts.append(f"<li>{action}</li>")
            parts.append("</ol>")
            if p.get("clinical_actions"):
                parts.append("<p><strong>Clinical guidance:</strong></p><ul>")
                for ca in p["clinical_actions"]:
                    parts.append(f"<li>{ca}</li>")
                parts.append("</ul>")
            parts.append("</details>")
    narratives = (insights_data or {}).get("narratives", [])
    nutrition_narratives = [n for n in narratives if n.get("id", "") in nutrition_narrative_ids]
    if nutrition_narratives:
        parts.append("<h3>Nutritional Gene Stories</h3>")
        for n in nutrition_narratives:
            genes_str = ", ".join(n["matched_genes"])
            parts.append(
                f'<details open class="rec-card" style="border-left:4px solid var(--accent2)">'
                f'<summary><strong>{n["title"]}</strong> '
                f'<span class="badge" style="background:var(--accent2)">'
                f'{len(n["matched_genes"])} genes</span></summary>'
                f'<p><strong>Genes:</strong> {genes_str}</p>'
                f'<p>{n["narrative"]}</p>'
                f'<p><strong>What this means for you:</strong> {n["practical"]}</p>'
                f'</details>'
            )
    if not parts:
        parts.append("<p>No nutrition-specific genetic findings identified.</p>")
    return "\n".join(parts)


def build_protective(recommendations_data, insights_data):
    """Build protective variants section (test-compatible wrapper)."""
    good_news = (recommendations_data or {}).get("good_news", [])
    protective = (insights_data or {}).get("protective_findings", [])
    if not good_news and not protective:
        return "<p>No protective variants identified.</p>"
    parts = []
    if good_news:
        parts.append('<div class="good-news-grid">')
        seen_genes = set()
        for g in good_news:
            if g["gene"] in seen_genes:
                continue
            seen_genes.add(g["gene"])
            desc = _clean_condition(g["description"]) if "|" in g.get("description", "") else g["description"]
            parts.append(
                f'<div class="good-news-card">'
                f'<strong>{g["gene"]}</strong>: {desc}</div>'
            )
        parts.append("</div>")
    if protective:
        parts.append("<h3>Research-Backed Protective Findings</h3>")
        for p in protective:
            status = p["status"].replace("_", " ").title()
            parts.append(
                f'<details class="rec-card" style="border-left:4px solid var(--green)">'
                f'<summary><strong>{p["gene"]}</strong> ({status}): {p["title"]}</summary>'
                f'<p>{p["finding"]}</p>'
                f'<div class="paper-refs">Ref: {p["reference"]}</div>'
                f'</details>'
            )
    return "\n".join(parts)


# =============================================================================
# HTML TEMPLATE — 12 sections, practical typography
# =============================================================================

HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Genetic Health Report{subject_title}</title>
<style>
:root {{
  --bg: #faf9f6; --fg: #2a2a2a; --accent: #2b6777; --accent2: #7b5ea7;
  --border: #d4d0c8; --code-bg: #f0ede6; --card-bg: #f5f4f0;
  --table-stripe: #f5f3ee; --warn: #bf4040; --green: #3a8a5c;
  --shadow: 0 1px 4px rgba(0,0,0,0.06);
  --body-font: "Charter", "Bitstream Charter", "Sitka Text", Cambria, serif;
  --heading-font: "Concourse", -apple-system, "Segoe UI", Roboto, "Helvetica Neue", sans-serif;
  --code-font: "Triplicate", "Source Code Pro", Menlo, Consolas, monospace;
}}
@media (prefers-color-scheme: dark) {{
  :root {{
    --bg: #1a1b1e; --fg: #e5e2db; --accent: #6db3c4; --accent2: #b094d0;
    --border: #3a3a42; --code-bg: #252528; --card-bg: #222225;
    --table-stripe: #252528; --warn: #d06050; --green: #5cc47a;
    --shadow: 0 1px 4px rgba(0,0,0,0.25);
  }}
}}
*, *::before, *::after {{ box-sizing: border-box; }}
body {{
  font-family: var(--body-font);
  font-size: 1.05rem;
  line-height: 1.65; color: var(--fg); background: var(--bg);
  max-width: 42em; margin: 0 auto; padding: 2.5em 1.5em;
  -webkit-font-smoothing: antialiased;
  text-rendering: optimizeLegibility;
}}
h1 {{
  font-family: var(--heading-font);
  font-size: 1.9em; font-weight: 700; letter-spacing: -0.01em;
  border-bottom: 3px solid var(--accent); padding-bottom: .3em; margin-top: 0;
}}
h2 {{
  font-family: var(--heading-font);
  font-size: 1.35em; font-weight: 600; letter-spacing: -0.01em;
  border-bottom: 2px solid var(--border); padding-bottom: .2em;
  margin-top: 2.5em; color: var(--accent);
}}
h3 {{
  font-family: var(--heading-font);
  font-size: 1.1em; font-weight: 600;
  margin-top: 1.5em;
}}
h4 {{
  font-family: var(--heading-font);
  font-size: 1em; font-weight: 600; margin-top: 1.2em;
}}
.inline-h3 {{ display: inline; margin: 0; }}
hr {{ border: none; border-top: 1px solid var(--border); margin: 2em 0; }}
a {{ color: var(--accent); text-decoration: underline; text-decoration-color: rgba(26,82,118,0.3); text-underline-offset: 2px; }}
a:hover {{ text-decoration-color: var(--accent); }}
code {{
  font-family: var(--code-font);
  background: var(--code-bg); padding: .12em .35em; border-radius: 3px;
  font-size: .85em;
}}
table {{
  border-collapse: collapse; width: 100%; margin: 1em 0;
  font-size: 0.9em;
}}
th, td {{ border: 1px solid var(--border); padding: .5em .7em; text-align: left; }}
th {{ background: var(--code-bg); font-weight: 600; font-family: var(--heading-font); font-size: .9em; }}
tr:nth-child(even) {{ background: var(--table-stripe); }}
ul, ol {{ padding-left: 1.3em; }}
li {{ margin: .35em 0; }}
p {{ margin: .6em 0; }}

/* ELI5 styling */
.eli5 {{
  background: var(--code-bg); border-left: 3px solid var(--accent2);
  padding: .6em 1em; margin: .8em 0; font-size: .92em;
  border-radius: 0 6px 6px 0; color: var(--fg);
}}
.eli5-inline {{ font-size: .85em; color: var(--accent2); font-style: italic; }}

/* Part headers */
.part-header {{
  font-family: var(--heading-font);
  font-size: .85em; font-weight: 700; text-transform: uppercase;
  letter-spacing: .12em; color: var(--accent2);
  margin-top: 3em; padding: .3em 0;
  border-bottom: 1px solid var(--accent2);
}}

/* Sections */
.section {{
  background: var(--card-bg); border: 1px solid var(--border);
  border-radius: 6px; padding: 1.8em 1.6em; margin: 1.8em 0;
  box-shadow: var(--shadow);
}}
.section:first-of-type {{
  border-top: 3px solid var(--accent);
}}
.section-number {{
  display: inline-block; background: var(--accent); color: #fff;
  width: 24px; height: 24px; border-radius: 50%; text-align: center;
  line-height: 24px; font-size: 12px; font-weight: 700;
  margin-right: .4em; vertical-align: middle;
  font-family: var(--heading-font);
}}

/* Finding cards */
.finding-card {{
  border-left: 4px solid var(--border); padding: .75em 1em;
  margin: .75em 0; background: var(--card-bg); border-radius: 0 6px 6px 0;
}}
.finding-card.high {{ border-left-color: #bf4040; }}
.finding-card.mod {{ border-left-color: #c49a20; }}
.finding-card.low {{ border-left-color: #3a8a5c; }}
.finding-card.info {{ border-left-color: #8993a0; }}
.finding-header {{ font-size: .92em; }}
.finding-desc {{ margin: .3em 0; font-size: .9em; }}
.finding-note {{ font-size: .85em; color: var(--accent2); }}
.mag-badge {{
  display: inline-block; padding: .1em .5em; border-radius: 12px;
  font-size: .78em; font-weight: bold; color: #fff;
  font-family: var(--heading-font);
}}
.mag-high {{ background: #bf4040; }}
.mag-mod {{ background: #c49a20; }}
.mag-low {{ background: #3a8a5c; }}
.mag-info {{ background: #8993a0; }}
.mag-dot {{
  display: inline-block; width: 10px; height: 10px; border-radius: 50%;
  margin-right: .3em; vertical-align: middle;
}}
.mag-dot.mag-high {{ background: #bf4040; }}
.mag-dot.mag-mod {{ background: #c49a20; }}
.mag-dot.mag-low {{ background: #3a8a5c; }}
.mag-dot.mag-info {{ background: #8993a0; }}

/* Key findings bullets */
.key-findings {{ list-style: none; padding-left: 0; }}
.key-findings li {{
  padding: .7em 1.1em; margin: .4em 0; border-radius: 4px;
  border-left: 4px solid var(--border);
  line-height: 1.55;
}}
.key-findings li strong {{ font-family: var(--heading-font); }}
.kf-red {{ border-left-color: #bf4040; background: rgba(191,64,64,0.05); }}
.kf-yellow {{ border-left-color: #c49a20; background: rgba(196,154,32,0.05); }}
.kf-green {{ border-left-color: #3a8a5c; background: rgba(58,138,92,0.04); }}
@media (prefers-color-scheme: dark) {{
  .kf-red {{ background: rgba(208,96,80,0.1); }}
  .kf-yellow {{ background: rgba(196,154,32,0.1); }}
  .kf-green {{ background: rgba(92,196,122,0.08); }}
}}

/* Badge */
.badge {{
  display: inline-block; background: var(--accent); color: #fff;
  padding: .1em .55em; border-radius: 12px; font-size: .75em;
  vertical-align: middle; margin-left: .3em;
  font-family: var(--heading-font);
}}

/* Charts */
.chart-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 1.5em; }}
@media (max-width: 700px) {{ .chart-grid {{ grid-template-columns: 1fr; }} }}
.chart {{ width: 100%; height: auto; color: var(--fg); }}
.gauge {{ width: 140px; height: 100px; display: inline-block; }}
.gauge-row {{ display: flex; gap: 1.5em; flex-wrap: wrap; }}

/* DB links */
.db-links {{ font-size: .78em; }}
.db-links a {{ margin: 0 .2em; }}
.paper-refs {{ font-size: .78em; margin-top: .3em; color: var(--accent2); }}
.rsid-grid {{ display: grid; grid-template-columns: repeat(auto-fill, minmax(320px, 1fr)); gap: .5em; font-size: .85em; }}

/* Doctor callout */
.doctor-callout {{
  background: var(--code-bg); border: 2px solid var(--accent);
  border-radius: 8px; padding: 1em; margin: 1em 0;
  font-weight: 600; text-align: center;
}}

/* Collapsible */
details {{ margin: .5em 0; }}
summary {{
  cursor: pointer; padding: .4em 0; user-select: none;
  list-style: none;
}}
summary::-webkit-details-marker {{ display: none; }}
summary::before {{
  content: "\\25B6"; display: inline-block; margin-right: .5em;
  transition: transform .2s; font-size: .8em;
}}
details[open] > summary::before {{ transform: rotate(90deg); }}
.category-section {{ margin: .75em 0; }}

/* TOC */
.toc {{ columns: 2; column-gap: 2em; font-size: .88em; font-family: var(--heading-font); }}
.toc a {{ display: block; padding: .25em 0; text-decoration: none; }}
.toc a:hover {{ text-decoration: underline; }}
@media (max-width: 600px) {{ .toc {{ columns: 1; }} }}

/* Doctor card print */
.doctor-card {{ page-break-before: always; }}
@media print {{
  body {{ max-width: 100%; font-size: 10pt; padding: 0.5em; }}
  .chart, .gauge, .chart-grid, .gauge-row, .toc,
  .no-print {{ display: none !important; }}
  .section {{ break-inside: avoid; box-shadow: none; border: 1px solid #ccc; }}
  .doctor-card {{ page-break-before: always; }}
  details {{ display: block; }}
  details > summary {{ display: none; }}
  a {{ color: inherit; }}
}}

/* Skip link */
.skip-link {{ position: absolute; left: -9999px; }}
.skip-link:focus {{
  left: 1em; top: 1em; background: var(--accent); color: #fff;
  padding: .5em 1em; z-index: 999; border-radius: 4px;
}}

/* Search & toolbar */
.toolbar {{
  display: flex; gap: .75em; align-items: center;
  margin: 1em 0; flex-wrap: wrap;
}}
#search-box {{
  flex: 1; min-width: 200px; padding: .5em .75em;
  border: 2px solid var(--border); border-radius: 6px;
  background: var(--bg); color: var(--fg); font-size: .92em;
  font-family: var(--body-font);
}}
#search-box:focus {{ border-color: var(--accent); outline: none; }}
.toolbar button {{
  padding: .45em 1em; border: 1px solid var(--border);
  border-radius: 6px; background: var(--code-bg); color: var(--fg);
  cursor: pointer; font-size: .82em; white-space: nowrap;
  font-family: var(--heading-font);
}}
.toolbar button:hover {{ background: var(--accent); color: #fff; }}

/* Sortable tables */
th.sortable {{ cursor: pointer; user-select: none; }}
th.sortable::after {{ content: " \\2195"; opacity: 0.3; }}
th.sortable.asc::after {{ content: " \\2191"; opacity: 1; }}
th.sortable.desc::after {{ content: " \\2193"; opacity: 1; }}

/* PRS gauges */
.prs-gauge {{ width: 160px; height: 105px; display: inline-block; }}

/* Recommendation cards */
.rec-card {{
  background: var(--card-bg); padding: .75em 1em; margin: .75em 0;
  border-radius: 0 8px 8px 0;
}}
.good-news-grid {{
  display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
  gap: .6em;
}}
.good-news-card {{
  background: var(--card-bg); border: 1px solid rgba(58,138,92,0.3);
  border-left: 3px solid var(--green); border-radius: 0 4px 4px 0;
  padding: .5em .9em; font-size: .88em; line-height: 1.5;
}}
.good-news-card strong {{ font-family: var(--heading-font); }}
</style>
</head>
<body>
<a href="#main" class="skip-link">Skip to content</a>

<header style="margin-bottom:2em">
<h1>Genetic Health Report{subject_title}</h1>
<div style="display:flex;flex-wrap:wrap;gap:.5em 1.5em;font-size:.85em;color:var(--accent2);font-family:var(--heading-font)">
<span><strong>{generated_date}</strong></span>
<span>{total_snps:,} SNPs analyzed</span>
<span>{num_findings} lifestyle findings</span>
<span>{num_pharmgkb} drug interactions</span>
</div>
</header>

<nav class="no-print">
<h2 id="toc-heading">Contents</h2>
<p class="part-header">Part 1 &mdash; What You Should Do</p>
<div class="toc">
<a href="#key-findings">1. Your Key Findings</a>
<a href="#action-plan">2. Your Action Plan</a>
<a href="#drug-guide">3. Drug &amp; Medication Guide</a>
<a href="#disease-risk">4. Disease Risk Overview</a>
<a href="#body-profile">5. Your Body Profile</a>
<a href="#mental-health">6. Mental Health &amp; Brain</a>
</div>
<p class="part-header">Part 2 &mdash; Details for Your Doctor</p>
<div class="toc">
<a href="#clinical-detail">7. Clinical Findings Detail</a>
<a href="#ancestry">8. Ancestry &amp; Population</a>
<a href="#nutrigenomics">9. Nutrigenomics Detail</a>
<a href="#quality">10. Data Quality</a>
<a href="#doctor-card">11. Doctor Card</a>
<a href="#references">12. References &amp; Links</a>
</div>
</nav>

<div class="toolbar no-print">
<input type="text" id="search-box" placeholder="Search findings (gene, rsID, keyword)..." aria-label="Search findings">
<button id="export-csv" title="Copy findings as CSV to clipboard">Export CSV</button>
</div>

<main id="main">

<p class="part-header no-print">Part 1 &mdash; What You Should Do</p>

<div class="section" id="key-findings">
<h2><span class="section-number">1</span> Your Key Findings</h2>
{key_findings_content}
</div>

<div class="section" id="action-plan">
<h2><span class="section-number">2</span> Your Action Plan</h2>
{action_plan_content}
</div>

<div class="section" id="drug-guide">
<h2><span class="section-number">3</span> Drug &amp; Medication Guide</h2>
{drug_guide_content}
</div>

<div class="section" id="disease-risk">
<h2><span class="section-number">4</span> Disease Risk Overview</h2>
{disease_risk_content}
</div>

<div class="section" id="body-profile">
<h2><span class="section-number">5</span> Your Body Profile</h2>
{body_profile_content}
</div>

<div class="section" id="mental-health">
<h2><span class="section-number">6</span> Mental Health &amp; Brain</h2>
{mental_health_content}
</div>

<p class="part-header no-print">Part 2 &mdash; Details for Your Doctor</p>

<div class="section" id="clinical-detail">
<h2><span class="section-number">7</span> Clinical Findings Detail</h2>
{clinical_detail_content}
</div>

<div class="section" id="ancestry">
<h2><span class="section-number">8</span> Ancestry &amp; Population Context</h2>
{ancestry_content}
</div>

<div class="section" id="nutrigenomics">
<h2><span class="section-number">9</span> Nutrigenomics Detail</h2>
{nutrigenomics_content}
</div>

<div class="section" id="quality">
<h2><span class="section-number">10</span> Data Quality</h2>
{quality_content}
</div>

<div class="section doctor-card" id="doctor-card">
<h2><span class="section-number">11</span> Doctor Card</h2>
<p><em>Print this page &mdash; this section is optimized for printing.</em></p>
{doctor_card_content}
</div>

<div class="section" id="references">
<h2><span class="section-number">12</span> References &amp; Links</h2>
{references_content}
</div>

</main>

<footer style="text-align:center;margin-top:3em;padding-top:1.5em;border-top:1px solid var(--border);font-size:.82em;color:var(--accent2);">
Generated by Genetic Health Analysis Pipeline &mdash; {generated_date}<br>
For informational purposes only. Not a clinical diagnosis.
</footer>

<script>
// Keyboard accessibility for summary elements
document.querySelectorAll('summary').forEach(function(s) {{
  s.setAttribute('tabindex', '0');
}});

// Smooth scroll for TOC links
document.querySelectorAll('.toc a').forEach(function(a) {{
  a.addEventListener('click', function(e) {{
    var target = document.querySelector(this.getAttribute('href'));
    if (target) {{
      e.preventDefault();
      target.scrollIntoView({{ behavior: 'smooth', block: 'start' }});
      history.pushState(null, null, this.getAttribute('href'));
    }}
  }});
}});

// --- Search / Filter ---
(function() {{
  var searchBox = document.getElementById('search-box');
  if (!searchBox) return;
  searchBox.addEventListener('input', function() {{
    var q = this.value.toLowerCase().trim();
    document.querySelectorAll('.finding-card').forEach(function(card) {{
      card.style.display = (!q || card.textContent.toLowerCase().indexOf(q) !== -1)
        ? '' : 'none';
    }});
    document.querySelectorAll('table').forEach(function(table) {{
      var rows = table.querySelectorAll('tr');
      for (var i = 1; i < rows.length; i++) {{
        rows[i].style.display = (!q || rows[i].textContent.toLowerCase().indexOf(q) !== -1)
          ? '' : 'none';
      }}
    }});
  }});
}})();

// --- Sortable Tables ---
(function() {{
  document.querySelectorAll('table').forEach(function(table) {{
    var headers = table.querySelectorAll('th');
    if (headers.length < 2) return;
    headers.forEach(function(th, colIdx) {{
      th.classList.add('sortable');
      th.addEventListener('click', function() {{
        var isAsc = th.classList.contains('asc');
        headers.forEach(function(h) {{ h.classList.remove('asc', 'desc'); }});
        th.classList.add(isAsc ? 'desc' : 'asc');
        var rows = Array.from(table.querySelectorAll('tr')).slice(1);
        rows.sort(function(a, b) {{
          var aText = (a.children[colIdx] || {{}}).textContent || '';
          var bText = (b.children[colIdx] || {{}}).textContent || '';
          var aNum = parseFloat(aText.replace(/[^0-9.\\-]/g, ''));
          var bNum = parseFloat(bText.replace(/[^0-9.\\-]/g, ''));
          if (!isNaN(aNum) && !isNaN(bNum)) {{
            return isAsc ? bNum - aNum : aNum - bNum;
          }}
          return isAsc ? bText.localeCompare(aText) : aText.localeCompare(bText);
        }});
        var tbody = table.querySelector('tbody') || table;
        rows.forEach(function(r) {{ tbody.appendChild(r); }});
      }});
    }});
  }});
}})();

// --- CSV Export ---
(function() {{
  var btn = document.getElementById('export-csv');
  if (!btn) return;
  btn.addEventListener('click', function() {{
    var lines = ['Gene,rsID,Genotype,Status,Category,Magnitude'];
    document.querySelectorAll('.finding-card').forEach(function(card) {{
      var header = card.querySelector('.finding-header');
      if (!header) return;
      var text = header.textContent;
      var parts = text.split(/\\s+—\\s+/);
      var gene = '', rsid = '', genotype = '', status = '';
      if (parts.length >= 1) {{
        var m = parts[0].match(/([A-Z0-9]+)\\s+(rs\\d+)\\s+(.+)/);
        if (m) {{ gene = m[1]; rsid = m[2]; genotype = m[3]; }}
      }}
      if (parts.length >= 2) status = parts[parts.length - 1].trim();
      var mag = '';
      var magEl = card.querySelector('.mag-badge');
      if (magEl) mag = magEl.textContent.trim();
      var cat = '';
      var section = card.closest('.category-section');
      if (section) {{
        var sum = section.querySelector('summary');
        if (sum) cat = sum.textContent.replace(/\\d+$/, '').trim();
      }}
      lines.push([gene, rsid, genotype, status, cat, mag].join(','));
    }});
    var csv = lines.join('\\n');
    if (navigator.clipboard) {{
      navigator.clipboard.writeText(csv).then(function() {{
        btn.textContent = 'Copied!';
        setTimeout(function() {{ btn.textContent = 'Export CSV'; }}, 2000);
      }});
    }}
  }});
}})();
</script>
</body>
</html>
"""


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 60)
    print("Genetic Health Report (Redesigned)")
    print("=" * 60)

    results_path = REPORTS_DIR / "comprehensive_results.json"

    if not results_path.exists():
        print(f"ERROR: {results_path} not found. Run run_full_analysis.py first.")
        sys.exit(1)

    print(f"\n>>> Loading {results_path.name}")
    data = load_json(results_path)
    findings = data.get("findings", [])
    pharmgkb_findings = data.get("pharmgkb_findings", [])
    summary = data.get("summary", {})
    ancestry_data = data.get("ancestry", {})
    prs_data = data.get("prs", {})
    epistasis_data = data.get("epistasis", [])
    recommendations_data = data.get("recommendations", {})
    quality_data = data.get("quality_metrics", {})
    blood_type_data = data.get("blood_type", {})
    mt_haplogroup_data = data.get("mt_haplogroup", {})
    star_alleles_data = data.get("star_alleles", {})
    apoe_data = data.get("apoe", {})
    acmg_data = data.get("acmg", {})
    carrier_screen_data = data.get("carrier_screen", {})
    traits_data = data.get("traits", {})
    insights_data = data.get("insights", {})
    disease_findings_data = data.get("disease_findings", {})
    polypharmacy_data = data.get("polypharmacy", {})
    longevity_data = data.get("longevity", {})
    sleep_data = data.get("sleep_profile", {})
    nutrigenomics_data = data.get("nutrigenomics", {})
    mental_health_data = data.get("mental_health", {})

    subject_name = data.get("subject_name", "")

    # Build each section
    print("\n>>> Building key findings")
    key_findings = build_key_findings(
        findings, recommendations_data, apoe_data, acmg_data,
        prs_data, star_alleles_data, longevity_data,
        traits_data=traits_data, blood_type_data=blood_type_data,
        sleep_data=sleep_data, mental_health_data=mental_health_data,
        ancestry_data=ancestry_data, nutrigenomics_data=nutrigenomics_data,
        disease_findings_data=disease_findings_data)

    print(">>> Building action plan")
    action_plan = build_action_plan(recommendations_data, insights_data)

    print(">>> Building drug guide")
    drug_guide = build_drug_guide(
        star_alleles_data, findings, pharmgkb_findings, polypharmacy_data)

    print(">>> Building disease risk overview")
    disease_risk = build_disease_risk_overview(prs_data, disease_findings_data, acmg_data)

    print(">>> Building body profile")
    body_profile = build_body_profile(
        traits_data, blood_type_data, sleep_data, longevity_data, findings)

    print(">>> Building mental health section")
    mental_health_html = build_mental_health_section(mental_health_data)

    print(">>> Building clinical detail")
    clinical_detail = build_clinical_detail(findings, epistasis_data, carrier_screen_data)

    print(">>> Building ancestry section")
    ancestry_html = build_ancestry_section(ancestry_data, mt_haplogroup_data)

    print(">>> Building nutrigenomics section")
    nutrigenomics_html = build_nutrigenomics_section(
        nutrigenomics_data, recommendations_data, insights_data)

    print(">>> Building quality section")
    quality_html = build_quality_section(quality_data)

    print(">>> Building doctor card")
    doctor_card = build_doctor_card(
        recommendations_data, star_alleles_data, apoe_data, acmg_data)

    print(">>> Building references")
    references = build_references(findings)

    # Assemble HTML
    print("\n>>> Assembling HTML report")
    subject_title = f" — {_esc(subject_name)}" if subject_name else ""
    html = HTML_TEMPLATE.format(
        generated_date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        total_snps=summary.get("total_snps", 0),
        num_findings=len(findings),
        num_pharmgkb=len(pharmgkb_findings),
        subject_title=subject_title,
        key_findings_content=key_findings,
        action_plan_content=action_plan,
        drug_guide_content=drug_guide,
        disease_risk_content=disease_risk,
        body_profile_content=body_profile,
        mental_health_content=mental_health_html,
        clinical_detail_content=clinical_detail,
        ancestry_content=ancestry_html,
        nutrigenomics_content=nutrigenomics_html,
        quality_content=quality_html,
        doctor_card_content=doctor_card,
        references_content=references,
    )

    output_path = REPORTS_DIR / "GENETIC_HEALTH_REPORT.html"
    output_path.write_text(html, encoding="utf-8")

    print(f"\n{'='*60}")
    print(f"Report generated: {output_path}")
    print(f"Size: {len(html):,} characters")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
