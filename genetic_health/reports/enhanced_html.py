"""
Enhanced All-in-One Genetic Health Report Generator

Produces a single self-contained HTML file merging:
  - Lifestyle/health findings (comprehensive_results.json)
  - Disease risk analysis (EXHAUSTIVE_DISEASE_RISK_REPORT.md)
  - Actionable protocol (ACTIONABLE_HEALTH_PROTOCOL_V3.md)
  - Personal health summary (PERSONAL_HEALTH_SUMMARY.md)

Output: reports/ENHANCED_HEALTH_REPORT.html
  - Zero external dependencies (all CSS/JS/SVG inline)
  - Dark/light mode via prefers-color-scheme
  - Collapsible sections (vanilla JS)
  - Print-optimized doctor card
  - Database links for every rsID
  - Curated paper references
"""

import json
import re
import sys
from datetime import datetime
from collections import defaultdict
from pathlib import Path



from ..config import REPORTS_DIR

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
    return '<span class="db-links">' + " · ".join(links) + "</span>"


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


def load_text(path):
    if path.exists():
        return path.read_text(encoding="utf-8")
    return ""


def parse_disease_report(md_text):
    """Parse the disease risk Markdown report into structured data."""
    sections = {
        "pathogenic": [],
        "unclear_inheritance": [],
        "risk_factors": [],
        "drug_response": [],
        "protective": [],
    }
    current = None

    for line in md_text.split("\n"):
        stripped = line.strip()

        if "Pathogenic Variants - Affected Status" in stripped:
            current = "pathogenic"
            continue
        elif "Pathogenic/Likely Pathogenic - Inheritance Unclear" in stripped:
            current = "unclear_inheritance"
            continue
        elif "Risk Factor Variants" in stripped:
            current = "risk_factors"
            continue
        elif "Drug Response Variants" in stripped:
            current = "drug_response"
            continue
        elif "Protective Variants" in stripped and "Good News" not in stripped:
            current = "protective"
            continue
        elif stripped.startswith("## Disclaimer") or stripped.startswith("## Executive Summary"):
            current = None
            continue

        if current is None:
            continue

        # Parse ### headers as variant entries (pathogenic / unclear)
        if current in ("pathogenic", "unclear_inheritance"):
            m = re.match(r"^###\s+(.+?)\s*-\s*(.+)", stripped)
            if m:
                gene = m.group(1).strip()
                condition = m.group(2).strip()
                sections[current].append({"gene": gene, "condition": condition, "details": {}})
                continue
            if sections[current]:
                entry = sections[current][-1]
                m2 = re.match(r"^-\s+\*\*(.+?)\*\*:?\s*(.*)", stripped)
                if m2:
                    key = m2.group(1).strip()
                    val = m2.group(2).strip()
                    entry["details"][key] = val

        # Parse bullet list entries (risk, drug, protective)
        if current in ("risk_factors", "drug_response", "protective"):
            m = re.match(r"^-\s+\*\*(.+?)\*\*\s*\((.+?)\):\s*`?([^`]*)`?\s*-\s*(.*)", stripped)
            if m:
                sections[current].append({
                    "gene": m.group(1).strip(),
                    "rsid": m.group(2).strip(),
                    "genotype": m.group(3).strip(),
                    "description": m.group(4).strip(),
                })

    return sections


def parse_protocol_summary(md_text):
    """Parse the actionable protocol for risk factor summary table."""
    risk_table = []
    in_risk_table = False
    for line in md_text.split("\n"):
        stripped = line.strip()
        if "| Condition |" in stripped or "| Area |" in stripped:
            in_risk_table = True
            continue
        if in_risk_table:
            if stripped.startswith("|---"):
                continue
            if not stripped.startswith("|"):
                in_risk_table = False
                continue
            parts = [p.strip() for p in stripped.split("|") if p.strip()]
            if len(parts) >= 2:
                risk_table.append({"condition": parts[0], "genes": parts[1],
                                   "notes": parts[2] if len(parts) > 2 else ""})
    return risk_table


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
        bar(10, high, "#ef4444", "High"),
        bar(48, mod, "#f59e0b", "Moderate"),
        bar(86, low, "#22c55e", "Low"),
        bar(124, info, "#94a3b8", "Info"),
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

    colors = [
        "#ef4444", "#f59e0b", "#22c55e", "#3b82f6", "#8b5cf6",
        "#ec4899", "#14b8a6", "#f97316", "#6366f1", "#84cc16",
        "#06b6d4", "#e11d48",
    ]
    total = sum(counts.values())
    cx, cy, r = 120, 120, 90
    inner_r = 55
    angle = -90  # start at top
    paths = []
    legend_items = []
    import math

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


def svg_risk_cards(personal_summary_text):
    """Colored risk overview cards based on personal summary risk table."""
    risks = [
        ("Blood Pressure", "#ef4444", "AGTR1+AGT+ADRB1"),
        ("Diabetes", "#f59e0b", "KCNJ11, HNF1A"),
        ("Thrombosis", "#f97316", "F13B, FGA"),
        ("Eye/Macular", "#8b5cf6", "C3"),
        ("IBD", "#3b82f6", "ATG16L1, NOD2"),
        ("Thyroid", "#14b8a6", "TG, CTLA4"),
    ]
    cards = []
    for i, (label, color, genes) in enumerate(risks):
        x = (i % 3) * 148
        y = (i // 3) * 80
        cards.append(
            f'<rect x="{x}" y="{y}" width="140" height="70" rx="8" '
            f'fill="{color}" opacity="0.15" stroke="{color}" stroke-width="1.5"/>'
            f'<text x="{x+70}" y="{y+28}" text-anchor="middle" '
            f'fill="{color}" font-size="14" font-weight="bold">{label}</text>'
            f'<text x="{x+70}" y="{y+48}" text-anchor="middle" '
            f'fill="currentColor" font-size="10">{genes}</text>'
        )
    return (
        '<svg viewBox="0 0 444 170" class="chart" role="img" '
        'aria-label="Risk area cards">'
        '<title>Risk Overview</title>' + "".join(cards) + "</svg>"
    )


def svg_metabolism_gauge(label, level, color):
    """Simple gauge for drug metabolism speed. level: 0-2 (slow/intermediate/fast)."""
    import math
    cx, cy, r = 70, 65, 50
    # Arc sweeps from 150° to 30° (going counterclockwise through the top)
    # In SVG coords (y-down): use negative angles for upper half
    # Points on the arc: (cx + r*cos(a), cy - r*sin(a))  [standard math angles]
    start_math = 150  # lower-left
    end_math = 30     # lower-right
    range_math = start_math - end_math  # 120° sweep

    # Pointer angle: 150° (slow) -> 90° (intermediate) -> 30° (fast)
    ptr_math = start_math - (level / 2) * range_math
    ptr_rad = math.radians(ptr_math)
    ptr_r = r - 10
    px = cx + ptr_r * math.cos(ptr_rad)
    py = cy - ptr_r * math.sin(ptr_rad)

    # Arc endpoints
    s_rad = math.radians(start_math)
    e_rad = math.radians(end_math)
    sx = cx + r * math.cos(s_rad)
    sy = cy - r * math.sin(s_rad)
    ex = cx + r * math.cos(e_rad)
    ey = cy - r * math.sin(e_rad)

    speed_labels = {0: "Slow", 1: "Intermediate", 2: "Fast"}
    speed_text = speed_labels.get(level, "Unknown")

    # SVG arc: sweep-flag=0 means counterclockwise (going through the top)
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

    import math
    proportions = ancestry_results["proportions"]
    colors = {"EUR": "#3b82f6", "AFR": "#f59e0b", "EAS": "#ef4444",
              "SAS": "#8b5cf6", "AMR": "#22c55e"}
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
        color = colors.get(pop, "#94a3b8")
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

    # Center text
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
    import math
    cx, cy, r = 80, 70, 55

    # Arc from 180° to 0° (left to right through top)
    # Zone boundaries: low <20th, avg 20-80, elevated 80-95, high >95
    zones = [
        (0.0, 0.2, "#22c55e"),    # low - green
        (0.2, 0.8, "#3b82f6"),    # average - blue
        (0.8, 0.95, "#f59e0b"),   # elevated - orange
        (0.95, 1.0, "#ef4444"),   # high - red
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

    # Pointer
    frac = max(0.0, min(1.0, percentile / 100.0))
    ptr_angle = math.radians(180 - frac * 180)
    ptr_r = r - 18
    px = cx + ptr_r * math.cos(ptr_angle)
    py = cy - ptr_r * math.sin(ptr_angle)

    cat_colors = {"low": "#22c55e", "average": "#3b82f6",
                  "elevated": "#f59e0b", "high": "#ef4444"}
    ptr_color = cat_colors.get(category, "#94a3b8")

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


def build_ancestry_section(ancestry_results):
    """Build HTML for ancestry estimation section."""
    if not ancestry_results or ancestry_results.get("markers_found", 0) == 0:
        return "<p>No ancestry-informative markers found in genome data.</p>"

    parts = []
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
    parts.append('<p style="font-size:.85em;color:var(--accent2)">'
                 'Based on ~55 ancestry-informative markers. '
                 'This is a rough superpopulation estimate.</p>')
    parts.append("</div>")
    parts.append("</div>")

    return "\n".join(parts)


def build_prs_section(prs_results):
    """Build HTML for polygenic risk scores section."""
    if not prs_results:
        return "<p>No PRS data available.</p>"

    parts = []

    # Ancestry disclaimer
    any_non_applicable = any(not r["ancestry_applicable"] for r in prs_results.values())
    if any_non_applicable:
        parts.append(
            '<div class="doctor-callout" style="border-color:var(--warn)">'
            'PRS models are calibrated on European-ancestry populations. '
            'Your ancestry profile is substantially non-European — interpret with caution.'
            '</div>'
        )

    # Gauge row
    parts.append('<div class="gauge-row" style="justify-content:center">')
    for cid, r in prs_results.items():
        short_name = r["name"].replace("Age-Related ", "").replace("Macular Degeneration", "AMD")
        parts.append(svg_prs_gauge(short_name, r["percentile"], r["risk_category"]))
    parts.append("</div>")

    # Summary table
    parts.append("<h3>Summary</h3>")
    parts.append('<table class="sortable"><tr><th>Condition</th><th>Percentile</th>'
                 '<th>Category</th><th>SNPs</th><th>Reference</th></tr>')
    for cid, r in prs_results.items():
        cat_class = {"low": "green", "average": "accent", "elevated": "warn", "high": "warn"}
        color = cat_class.get(r["risk_category"], "")
        parts.append(
            f'<tr><td><strong>{r["name"]}</strong></td>'
            f'<td>{r["percentile"]:.0f}th</td>'
            f'<td>{r["risk_category"].title()}</td>'
            f'<td>{r["snps_found"]}/{r["snps_total"]}</td>'
            f'<td style="font-size:.8em">{r["reference"]}</td></tr>'
        )
    parts.append("</table>")

    # Elevated conditions detail
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
                        f'<tr><td>{s["gene"]}</td><td><code>{s["rsid"]}</code></td>'
                        f'<td>{s["copies"]}</td><td>{s["contribution"]:.3f}</td></tr>'
                    )
                parts.append("</table>")
            parts.append("</details>")

    parts.append(
        '<p style="font-size:.85em;color:var(--accent2)">'
        'PRS estimates relative genetic risk from common variants. '
        'Lifestyle, environment, and rare variants also affect risk. '
        'Not a clinical diagnosis.</p>'
    )

    return "\n".join(parts)


# =============================================================================
# SECTION BUILDERS
# =============================================================================

def build_eli5(findings, disease_sections, personal_text):
    """Build ELI5 (explain like I'm 5) summary."""
    high_impact = [f for f in findings if f.get("magnitude", 0) >= 3]

    # Detect key topics from data
    has_bp = any(f["gene"] in ("AGTR1", "AGT", "ADRB1", "ACE", "GNB3")
                 for f in findings if f.get("magnitude", 0) >= 1)
    has_drug = any(f["gene"].startswith("CYP") for f in findings)
    has_lung = "SERPINA1" in personal_text
    has_caffeine = any(f.get("category") == "Caffeine Response" for f in findings)

    parts = []
    parts.append(
        "<p>Your body has a recipe book called DNA. We read yours and "
        "found some interesting things. Here's the simple version:</p>"
    )

    items = []
    if has_lung:
        items.append(
            "Your lungs have a little less of a protective shield than most people "
            "(a protein called AAT). This means taking extra care of your lungs "
            "is really important — no smoke, clean air, and regular check-ups."
        )
    if has_bp:
        items.append(
            "Several of your genes point to blood pressure wanting to run a bit "
            "high. The good news: eating less salt, exercising, and monitoring "
            "your blood pressure can make a big difference."
        )
    if has_drug:
        items.append(
            "Your body processes some medicines differently than average. If a "
            "doctor ever prescribes you warfarin (a blood thinner) or certain "
            "other drugs, they should know you may need a smaller dose."
        )
    if has_caffeine:
        items.append(
            "Coffee and caffeine stick around in your body longer than most "
            "people, and you're also more likely to feel jittery from it. "
            "Morning-only coffee is a good rule for you."
        )

    # Protective good news
    items.append(
        "Good news: you carry protective variants that lower your risk for "
        "Alzheimer's disease, nicotine addiction, and coronary heart disease."
    )

    parts.append("<ul>" + "".join(f"<li>{i}</li>" for i in items) + "</ul>")
    parts.append(
        '<p class="doctor-callout">Show this report to your doctor! '
        "They can help you use this information wisely.</p>"
    )
    return "\n".join(parts)


def build_dashboard(findings, personal_text):
    """Build dashboard section with SVG charts."""
    parts = []

    parts.append('<div class="chart-grid">')
    parts.append("<div>")
    parts.append("<h3>Impact Distribution</h3>")
    parts.append(svg_impact_bar(findings))
    parts.append("</div>")

    parts.append("<div>")
    parts.append("<h3>Findings by Category</h3>")
    parts.append(svg_category_donut(findings))
    parts.append("</div>")
    parts.append("</div>")  # end chart-grid

    parts.append("<h3>Risk Overview</h3>")
    parts.append(svg_risk_cards(personal_text))

    # Drug metabolism gauges
    parts.append("<h3>Drug Metabolism</h3>")
    parts.append('<div class="gauge-row">')

    # Determine CYP2C9 level from findings
    cyp2c9_level = 1  # default intermediate
    for f in findings:
        if f["gene"] == "CYP2C9":
            if f["status"] in ("poor", "significantly_reduced"):
                cyp2c9_level = 0
            elif f["status"] in ("intermediate",):
                cyp2c9_level = 1
            else:
                cyp2c9_level = 2

    cyp1a2_level = 1
    for f in findings:
        if f["gene"] == "CYP1A2":
            if f["status"] in ("slow",):
                cyp1a2_level = 0
            elif f["status"] in ("intermediate",):
                cyp1a2_level = 1
            else:
                cyp1a2_level = 2

    parts.append(svg_metabolism_gauge("CYP2C9", cyp2c9_level, "#ef4444"))
    parts.append(svg_metabolism_gauge("CYP1A2", cyp1a2_level, "#3b82f6"))
    parts.append("</div>")

    return "\n".join(parts)


def build_critical_findings(personal_text):
    """Build critical findings from personal summary."""
    # Extract ## Critical Findings and ## Asthma-Specific sections
    parts = []

    # Parse sections from personal summary
    current_section = None
    section_lines = []

    for line in personal_text.split("\n"):
        if line.startswith("## Critical Findings"):
            current_section = "critical"
            section_lines = []
            continue
        elif line.startswith("## Asthma") or line.startswith("## Nutrition") or line.startswith("## Exercise"):
            if current_section == "critical" and section_lines:
                parts.append(_personal_section_to_html(section_lines))
            current_section = None
            section_lines = []
            continue

        if current_section == "critical":
            section_lines.append(line)

    if current_section == "critical" and section_lines:
        parts.append(_personal_section_to_html(section_lines))

    return "\n".join(parts) if parts else "<p>See lifestyle findings below.</p>"


def build_asthma_section(personal_text):
    """Extract asthma-specific findings from personal summary."""
    parts = []
    in_section = False
    lines = []

    for line in personal_text.split("\n"):
        if line.startswith("## Asthma"):
            in_section = True
            lines = []
            continue
        elif in_section and line.startswith("## "):
            break
        if in_section:
            lines.append(line)

    if lines:
        return _personal_section_to_html(lines)
    return "<p>No asthma-specific findings in personal summary.</p>"


def _personal_section_to_html(lines):
    """Convert personal summary markdown lines to HTML."""
    html = []
    in_list = False
    in_table = False
    table_rows = []

    for line in lines:
        stripped = line.strip()
        if not stripped:
            if in_list:
                html.append("</ul>")
                in_list = False
            if in_table:
                html.append(_render_table(table_rows))
                in_table = False
                table_rows = []
            continue

        if stripped.startswith("---"):
            continue

        # Table
        if stripped.startswith("|"):
            in_table = True
            if not stripped.startswith("|---"):
                table_rows.append(stripped)
            continue

        # Headers
        m = re.match(r"^(#{1,6})\s+(.*)", stripped)
        if m:
            level = min(len(m.group(1)) + 1, 6)  # shift down one level
            text = _inline_md(m.group(2))
            html.append(f"<h{level}>{text}</h{level}>")
            continue

        # List items
        if stripped.startswith("- "):
            if not in_list:
                html.append("<ul>")
                in_list = True
            html.append(f"<li>{_inline_md(stripped[2:])}</li>")
            continue

        # Paragraph
        html.append(f"<p>{_inline_md(stripped)}</p>")

    if in_list:
        html.append("</ul>")
    if in_table and table_rows:
        html.append(_render_table(table_rows))

    return "\n".join(html)


def _render_table(rows):
    """Render markdown table rows to HTML."""
    if not rows:
        return ""
    parts = ['<table>']
    for i, row in enumerate(rows):
        cells = [c.strip() for c in row.split("|") if c.strip()]
        tag = "th" if i == 0 else "td"
        tr = "".join(f"<{tag}>{_inline_md(c)}</{tag}>" for c in cells)
        parts.append(f"<tr>{tr}</tr>")
    parts.append("</table>")
    return "\n".join(parts)


def _inline_md(text):
    """Convert inline markdown (bold, code, links) to HTML."""
    text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)", r'<a href="\2" target="_blank">\1</a>', text)
    return text


def build_lifestyle_findings(findings):
    """Build all lifestyle findings grouped by category with collapsible sections."""
    categories = sorted(set(f.get("category", "Other") for f in findings))
    parts = []

    for cat in categories:
        cat_findings = sorted(
            [f for f in findings if f.get("category") == cat],
            key=lambda x: -x.get("magnitude", 0),
        )
        if not cat_findings:
            continue

        cat_id = cat.lower().replace(" ", "-").replace("/", "-")
        parts.append(
            f'<details class="category-section" open>'
            f'<summary><h3 class="inline-h3">{cat} '
            f'<span class="badge">{len(cat_findings)}</span></h3></summary>'
        )

        for f in cat_findings:
            mag = f.get("magnitude", 0)
            mag_class = "high" if mag >= 3 else "mod" if mag == 2 else "low" if mag == 1 else "info"
            rsid = f.get("rsid", "")
            gene = f.get("gene", "Unknown")
            genotype = f.get("genotype", "")
            status = f.get("status", "").replace("_", " ").title()
            desc = f.get("description", "")
            note = f.get("note", "")

            parts.append(f'<div class="finding-card {mag_class}">')
            parts.append(
                f'<div class="finding-header">'
                f'<span class="mag-badge mag-{mag_class}">{mag}/6</span> '
                f'<strong>{gene}</strong> '
                f'<code>{rsid}</code> — '
                f'<code>{genotype}</code> — {status}'
                f'</div>'
            )
            parts.append(f'<p class="finding-desc">{desc}</p>')
            if note:
                parts.append(f'<p class="finding-note">Note: {note}</p>')
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

        parts.append(f'<details><summary><strong>{pathway_name}</strong></summary><ul>')
        for pf in pathway_findings:
            mag = pf.get("magnitude", 0)
            mag_class = "high" if mag >= 3 else "mod" if mag == 2 else "low" if mag == 1 else "info"
            parts.append(
                f'<li><span class="mag-dot mag-{mag_class}"></span> '
                f'<strong>{pf["gene"]}</strong>: '
                f'{pf.get("status","").replace("_"," ").title()}</li>'
            )
        parts.append("</ul></details>")

    return "\n".join(parts)


def build_disease_risk(disease_sections):
    """Build disease risk section from parsed data."""
    parts = []

    # Pathogenic
    path = disease_sections.get("pathogenic", [])
    unclear = disease_sections.get("unclear_inheritance", [])
    parts.append(
        f'<p><strong>{len(path)}</strong> pathogenic (affected), '
        f'<strong>{len(unclear)}</strong> unclear inheritance</p>'
    )

    if path:
        parts.append('<details open><summary><strong>Pathogenic — Affected Status '
                     f'({len(path)})</strong></summary>')
        parts.append('<table><tr><th>Gene</th><th>Condition</th>'
                     '<th>Genotype</th><th>Confidence</th></tr>')
        for v in path:
            d = v.get("details", {})
            gt = d.get("Genotype", "")
            conf = d.get("Confidence", "")
            parts.append(
                f'<tr><td><strong>{v["gene"]}</strong></td>'
                f'<td>{v["condition"]}</td>'
                f'<td><code>{gt}</code></td>'
                f'<td>{conf}</td></tr>'
            )
        parts.append("</table></details>")

    if unclear:
        parts.append(f'<details><summary><strong>Inheritance Unclear ({len(unclear)})</strong></summary>')
        parts.append('<table><tr><th>Gene</th><th>Condition</th>'
                     '<th>Genotype</th><th>Confidence</th></tr>')
        for v in unclear:
            d = v.get("details", {})
            gt = d.get("Genotype", "")
            conf = d.get("Confidence", "")
            parts.append(
                f'<tr><td><strong>{v["gene"]}</strong></td>'
                f'<td>{v["condition"]}</td>'
                f'<td><code>{gt}</code></td>'
                f'<td>{conf}</td></tr>'
            )
        parts.append("</table></details>")

    # Risk factors
    risk = disease_sections.get("risk_factors", [])
    if risk:
        parts.append(
            f'<details><summary><strong>Risk Factor Variants ({len(risk)})</strong></summary>'
        )
        parts.append('<table><tr><th>Gene</th><th>RSID</th><th>Genotype</th>'
                     '<th>Description</th></tr>')
        for v in risk:
            parts.append(
                f'<tr><td><strong>{v["gene"]}</strong></td>'
                f'<td><code>{v["rsid"]}</code></td>'
                f'<td><code>{v["genotype"]}</code></td>'
                f'<td>{v["description"][:120]}...</td></tr>'
            )
        parts.append("</table></details>")

    # Drug response
    drug = disease_sections.get("drug_response", [])
    if drug:
        parts.append(
            f'<details><summary><strong>Drug Response Variants ({len(drug)})</strong></summary>'
        )
        parts.append('<table><tr><th>Gene</th><th>RSID</th><th>Genotype</th>'
                     '<th>Response</th></tr>')
        for v in drug:
            parts.append(
                f'<tr><td><strong>{v["gene"]}</strong></td>'
                f'<td><code>{v["rsid"]}</code></td>'
                f'<td><code>{v["genotype"]}</code></td>'
                f'<td>{v["description"][:120]}...</td></tr>'
            )
        parts.append("</table></details>")

    return "\n".join(parts)


def build_drug_gene(findings, pharmgkb_findings, disease_drug):
    """Build drug-gene interaction section."""
    parts = []

    if pharmgkb_findings:
        parts.append("<h3>PharmGKB Annotations</h3>")
        parts.append('<table><tr><th>Gene</th><th>RSID</th><th>Level</th>'
                     '<th>Drugs</th><th>Genotype</th></tr>')
        for p in pharmgkb_findings:
            parts.append(
                f'<tr><td><strong>{p["gene"]}</strong></td>'
                f'<td><code>{p["rsid"]}</code></td>'
                f'<td>{p["level"]}</td>'
                f'<td>{p["drugs"]}</td>'
                f'<td><code>{p["genotype"]}</code></td></tr>'
            )
        parts.append("</table>")

    # ClinVar drug response from disease sections
    if disease_drug:
        parts.append("<h3>ClinVar Drug Response Variants</h3>")
        parts.append('<table><tr><th>Gene</th><th>RSID</th><th>Genotype</th>'
                     '<th>Response</th></tr>')
        for v in disease_drug[:20]:  # show first 20
            parts.append(
                f'<tr><td><strong>{v["gene"]}</strong></td>'
                f'<td><code>{v["rsid"]}</code></td>'
                f'<td><code>{v["genotype"]}</code></td>'
                f'<td>{v["description"][:100]}...</td></tr>'
            )
        if len(disease_drug) > 20:
            parts.append(
                f'<tr><td colspan="4"><em>...and {len(disease_drug)-20} more</em></td></tr>'
            )
        parts.append("</table>")

    parts.append(
        '<div class="doctor-callout">Share this entire drug-gene section with '
        "every prescribing physician. Print it or save as PDF.</div>"
    )

    return "\n".join(parts)


def build_nutrition_section(personal_text):
    """Extract nutrition, supplement, and exercise info from personal summary."""
    sections_to_extract = ["Nutrition & Supplements", "Exercise", "Skin"]
    parts = []
    current = None
    lines = []

    for line in personal_text.split("\n"):
        if any(line.startswith(f"## {s}") for s in sections_to_extract):
            if current and lines:
                parts.append(_personal_section_to_html(lines))
            current = line
            lines = []
            parts.append(f"<h3>{line.replace('## ', '')}</h3>")
            continue
        elif current and line.startswith("## "):
            if lines:
                parts.append(_personal_section_to_html(lines))
            current = None
            lines = []
            continue
        if current:
            lines.append(line)

    if current and lines:
        parts.append(_personal_section_to_html(lines))

    return "\n".join(parts) if parts else "<p>See personal summary for details.</p>"


def build_monitoring(personal_text):
    """Extract monitoring schedule from personal summary."""
    in_section = False
    lines = []
    for line in personal_text.split("\n"):
        if line.startswith("## Monitoring"):
            in_section = True
            lines = []
            continue
        elif in_section and line.startswith("## "):
            break
        if in_section:
            lines.append(line)

    if lines:
        return _personal_section_to_html(lines)
    return "<p>See actionable protocol for monitoring details.</p>"


def build_protective(disease_sections, personal_text):
    """Build protective variants section."""
    prot = disease_sections.get("protective", [])
    parts = []

    # Also pull from personal summary
    in_prot = False
    personal_prot_lines = []
    for line in personal_text.split("\n"):
        if "Protective Variants" in line and line.startswith("## "):
            in_prot = True
            continue
        elif in_prot and line.startswith("## "):
            break
        if in_prot:
            personal_prot_lines.append(line)

    if personal_prot_lines:
        parts.append(_personal_section_to_html(personal_prot_lines))

    if prot:
        parts.append(f'<h3>All {len(prot)} ClinVar Protective Variants</h3>')
        parts.append('<table><tr><th>Gene</th><th>RSID</th><th>Genotype</th>'
                     '<th>Protection</th></tr>')
        for v in prot:
            parts.append(
                f'<tr><td><strong>{v["gene"]}</strong></td>'
                f'<td><code>{v["rsid"]}</code></td>'
                f'<td><code>{v["genotype"]}</code></td>'
                f'<td>{v["description"][:150]}</td></tr>'
            )
        parts.append("</table>")

    return "\n".join(parts)


def build_doctor_card(personal_text):
    """Build print-optimized doctor card."""
    parts = []
    in_doctor = False
    lines = []

    for line in personal_text.split("\n"):
        if "What to Share with Your Doctors" in line:
            in_doctor = True
            lines = []
            continue
        elif in_doctor and line.startswith("## Important"):
            break
        elif in_doctor and line.startswith("## ") and "Disclaim" not in line:
            break
        if in_doctor:
            lines.append(line)

    if lines:
        parts.append(_personal_section_to_html(lines))
    else:
        parts.append("<p>See personal summary for doctor-specific notes.</p>")

    return "\n".join(parts)


def build_references(findings):
    """Build references section with all paper links and database links for every rsID."""
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
            parts.append(f"<div><code>{rsid}</code> {db_links_html(rsid)}</div>")
        parts.append("</div>")

    # Methodology
    parts.append("<h3>Methodology &amp; Disclaimers</h3>")
    parts.append(
        "<ul>"
        "<li>Lifestyle findings from curated SNP database (~200 variants)</li>"
        "<li>Disease risk from ClinVar (~341K variants scanned)</li>"
        "<li>Drug interactions from PharmGKB clinical annotations</li>"
        "<li>Only true SNPs analyzed (indels filtered to prevent false positives)</li>"
        "<li>This report is for <strong>informational purposes only</strong> — not a clinical diagnosis</li>"
        "<li>Genetic associations are probabilistic, not deterministic</li>"
        "<li>Always consult healthcare providers before making medical decisions</li>"
        "</ul>"
    )

    return "\n".join(parts)


# =============================================================================
# HTML TEMPLATE
# =============================================================================

HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Enhanced Genetic Health Report</title>
<style>
:root {{
  --bg: #ffffff; --fg: #1a1a2e; --accent: #2563eb; --accent2: #7c3aed;
  --border: #e2e8f0; --code-bg: #f1f5f9; --card-bg: #f8fafc;
  --table-stripe: #f8fafc; --warn: #dc2626; --green: #16a34a;
  --shadow: 0 1px 3px rgba(0,0,0,0.08);
}}
@media (prefers-color-scheme: dark) {{
  :root {{
    --bg: #0f172a; --fg: #e2e8f0; --accent: #60a5fa; --accent2: #a78bfa;
    --border: #334155; --code-bg: #1e293b; --card-bg: #1e293b;
    --table-stripe: #1e293b; --warn: #f87171; --green: #4ade80;
    --shadow: 0 1px 3px rgba(0,0,0,0.3);
  }}
}}
*, *::before, *::after {{ box-sizing: border-box; }}
body {{
  font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
  line-height: 1.7; color: var(--fg); background: var(--bg);
  max-width: 62em; margin: 0 auto; padding: 2em 1.5em;
}}
h1 {{ border-bottom: 3px solid var(--accent); padding-bottom: .3em; margin-top: 0; }}
h2 {{
  border-bottom: 2px solid var(--border); padding-bottom: .2em;
  margin-top: 2.5em; color: var(--accent);
}}
h3 {{ margin-top: 1.5em; }}
.inline-h3 {{ display: inline; margin: 0; }}
hr {{ border: none; border-top: 1px solid var(--border); margin: 2em 0; }}
a {{ color: var(--accent); text-decoration: none; }}
a:hover {{ text-decoration: underline; }}
code {{
  background: var(--code-bg); padding: .15em .35em; border-radius: 3px;
  font-size: .88em;
}}
table {{
  border-collapse: collapse; width: 100%; margin: 1em 0;
  font-size: 0.92em;
}}
th, td {{ border: 1px solid var(--border); padding: .45em .65em; text-align: left; }}
th {{ background: var(--code-bg); font-weight: 600; }}
tr:nth-child(even) {{ background: var(--table-stripe); }}
ul, ol {{ padding-left: 1.5em; }}
li {{ margin: .3em 0; }}

/* Sections */
.section {{
  background: var(--card-bg); border: 1px solid var(--border);
  border-radius: 10px; padding: 1.5em; margin: 1.5em 0;
  box-shadow: var(--shadow);
}}
.section-number {{
  display: inline-block; background: var(--accent); color: #fff;
  width: 28px; height: 28px; border-radius: 50%; text-align: center;
  line-height: 28px; font-size: 14px; font-weight: bold;
  margin-right: .5em; vertical-align: middle;
}}

/* Finding cards */
.finding-card {{
  border-left: 4px solid var(--border); padding: .75em 1em;
  margin: .75em 0; background: var(--card-bg); border-radius: 0 6px 6px 0;
}}
.finding-card.high {{ border-left-color: #ef4444; }}
.finding-card.mod {{ border-left-color: #f59e0b; }}
.finding-card.low {{ border-left-color: #22c55e; }}
.finding-card.info {{ border-left-color: #94a3b8; }}
.finding-header {{ font-size: .95em; }}
.finding-desc {{ margin: .3em 0; font-size: .9em; }}
.finding-note {{ font-size: .85em; color: var(--accent2); }}
.mag-badge {{
  display: inline-block; padding: .1em .5em; border-radius: 12px;
  font-size: .8em; font-weight: bold; color: #fff;
}}
.mag-high {{ background: #ef4444; }}
.mag-mod {{ background: #f59e0b; }}
.mag-low {{ background: #22c55e; }}
.mag-info {{ background: #94a3b8; }}
.mag-dot {{
  display: inline-block; width: 10px; height: 10px; border-radius: 50%;
  margin-right: .3em; vertical-align: middle;
}}
.mag-dot.mag-high {{ background: #ef4444; }}
.mag-dot.mag-mod {{ background: #f59e0b; }}
.mag-dot.mag-low {{ background: #22c55e; }}
.mag-dot.mag-info {{ background: #94a3b8; }}

/* Badge */
.badge {{
  display: inline-block; background: var(--accent); color: #fff;
  padding: .1em .55em; border-radius: 12px; font-size: .75em;
  vertical-align: middle; margin-left: .3em;
}}

/* Charts */
.chart-grid {{ display: grid; grid-template-columns: 1fr 1fr; gap: 1.5em; }}
@media (max-width: 700px) {{ .chart-grid {{ grid-template-columns: 1fr; }} }}
.chart {{ width: 100%; height: auto; color: var(--fg); }}
.gauge {{ width: 140px; height: 100px; display: inline-block; }}
.gauge-row {{ display: flex; gap: 1.5em; flex-wrap: wrap; }}

/* DB links */
.db-links {{ font-size: .8em; }}
.db-links a {{ margin: 0 .2em; }}
.paper-refs {{ font-size: .8em; margin-top: .3em; color: var(--accent2); }}
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
.toc {{ columns: 2; column-gap: 2em; font-size: .9em; }}
.toc a {{ display: block; padding: .2em 0; }}
@media (max-width: 600px) {{ .toc {{ columns: 1; }} }}

/* Doctor card print */
.doctor-card {{
  page-break-before: always;
}}
@media print {{
  body {{ max-width: 100%; font-size: 10pt; padding: 0.5em; }}
  .chart, .gauge, .chart-grid, .gauge-row, .toc,
  .no-print {{ display: none !important; }}
  .section {{ break-inside: avoid; box-shadow: none; border: 1px solid #ccc; }}
  .doctor-card {{ page-break-before: always; }}
  details {{ display: block; }}
  details > summary {{ display: none; }}
  a {{ color: inherit; text-decoration: underline; }}
}}

/* Skip link for a11y */
.skip-link {{
  position: absolute; left: -9999px;
}}
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
  background: var(--bg); color: var(--fg); font-size: .95em;
}}
#search-box:focus {{ border-color: var(--accent); outline: none; }}
.toolbar button {{
  padding: .45em 1em; border: 1px solid var(--border);
  border-radius: 6px; background: var(--code-bg); color: var(--fg);
  cursor: pointer; font-size: .85em; white-space: nowrap;
}}
.toolbar button:hover {{ background: var(--accent); color: #fff; }}

/* Sortable tables */
th.sortable {{ cursor: pointer; user-select: none; }}
th.sortable::after {{ content: " \\2195"; opacity: 0.3; }}
th.sortable.asc::after {{ content: " \\2191"; opacity: 1; }}
th.sortable.desc::after {{ content: " \\2193"; opacity: 1; }}

/* PRS gauges */
.prs-gauge {{ width: 160px; height: 105px; display: inline-block; }}
</style>
</head>
<body>
<a href="#main" class="skip-link">Skip to content</a>

<h1>Enhanced Genetic Health Report</h1>
<p><strong>Generated:</strong> {generated_date} &nbsp;|&nbsp;
<strong>SNPs analyzed:</strong> {total_snps:,} &nbsp;|&nbsp;
<strong>Lifestyle findings:</strong> {num_findings} &nbsp;|&nbsp;
<strong>Drug interactions:</strong> {num_pharmgkb}</p>

<nav class="no-print">
<h2 id="toc-heading">Table of Contents</h2>
<div class="toc">
<a href="#eli5">1. Simple Summary</a>
<a href="#dashboard">2. Dashboard</a>
<a href="#ancestry">3. Ancestry</a>
<a href="#prs">4. Polygenic Risk Scores</a>
<a href="#critical">5. Critical Findings</a>
<a href="#asthma">6. Asthma &amp; Medications</a>
<a href="#lifestyle">7. Lifestyle Findings</a>
<a href="#disease">8. Disease Risk</a>
<a href="#drugs">9. Drug-Gene Interactions</a>
<a href="#nutrition">10. Nutrition &amp; Lifestyle</a>
<a href="#monitoring">11. Monitoring Schedule</a>
<a href="#protective">12. Protective Variants</a>
<a href="#doctor-card">13. Doctor Card</a>
<a href="#references">14. References &amp; Links</a>
</div>
</nav>

<div class="toolbar no-print">
<input type="text" id="search-box" placeholder="Search findings (gene, rsID, keyword)..." aria-label="Search findings">
<button id="export-csv" title="Copy findings as CSV to clipboard">Export CSV</button>
</div>

<main id="main">

<div class="section" id="eli5">
<h2><span class="section-number">1</span> The Simple Version</h2>
{eli5_content}
</div>

<div class="section no-print" id="dashboard">
<h2><span class="section-number">2</span> Dashboard</h2>
{dashboard_content}
</div>

<div class="section" id="ancestry">
<h2><span class="section-number">3</span> Ancestry Estimation</h2>
{ancestry_content}
</div>

<div class="section" id="prs">
<h2><span class="section-number">4</span> Polygenic Risk Scores</h2>
{prs_content}
</div>

<div class="section" id="critical">
<h2><span class="section-number">5</span> Critical Findings</h2>
{critical_content}
</div>

<div class="section" id="asthma">
<h2><span class="section-number">6</span> Asthma &amp; Medications</h2>
{asthma_content}
</div>

<div class="section" id="lifestyle">
<h2><span class="section-number">7</span> All Lifestyle Findings</h2>
{lifestyle_content}
</div>

<div class="section" id="disease">
<h2><span class="section-number">8</span> Disease Risk</h2>
{disease_content}
</div>

<div class="section" id="drugs">
<h2><span class="section-number">9</span> Drug-Gene Interactions</h2>
{drugs_content}
</div>

<div class="section" id="nutrition">
<h2><span class="section-number">10</span> Nutrition, Supplements &amp; Lifestyle</h2>
{nutrition_content}
</div>

<div class="section" id="monitoring">
<h2><span class="section-number">11</span> Monitoring Schedule</h2>
{monitoring_content}
</div>

<div class="section" id="protective">
<h2><span class="section-number">12</span> Protective Variants (Good News)</h2>
{protective_content}
</div>

<div class="section doctor-card" id="doctor-card">
<h2><span class="section-number">13</span> Doctor Card</h2>
<p><em>Print this page — only this section will appear in the printout.</em></p>
{doctor_card_content}
</div>

<div class="section" id="references">
<h2><span class="section-number">14</span> References &amp; Links</h2>
{references_content}
</div>

</main>

<footer style="text-align:center;margin-top:3em;font-size:.85em;color:var(--border);">
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
    // Filter finding cards
    document.querySelectorAll('.finding-card').forEach(function(card) {{
      card.style.display = (!q || card.textContent.toLowerCase().indexOf(q) !== -1)
        ? '' : 'none';
    }});
    // Filter table rows (skip header rows)
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
          var aNum = parseFloat(aText.replace(/[^0-9.\-]/g, ''));
          var bNum = parseFloat(bText.replace(/[^0-9.\-]/g, ''));
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
      var parts = text.split(/\s+—\s+/);
      var gene = '', rsid = '', genotype = '', status = '';
      if (parts.length >= 1) {{
        var m = parts[0].match(/([A-Z0-9]+)\s+(rs\d+)\s+(.+)/);
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
    print("Enhanced All-in-One Genetic Health Report")
    print("=" * 60)

    # Load data sources
    results_path = REPORTS_DIR / "comprehensive_results.json"
    disease_path = REPORTS_DIR / "EXHAUSTIVE_DISEASE_RISK_REPORT.md"
    protocol_path = REPORTS_DIR / "ACTIONABLE_HEALTH_PROTOCOL_V3.md"
    personal_path = REPORTS_DIR / "PERSONAL_HEALTH_SUMMARY.md"

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

    print(f">>> Loading {disease_path.name}")
    disease_text = load_text(disease_path)
    disease_sections = parse_disease_report(disease_text)

    print(f">>> Loading {protocol_path.name}")
    protocol_text = load_text(protocol_path)

    print(f">>> Loading {personal_path.name}")
    personal_text = load_text(personal_path)

    # Build each section
    print("\n>>> Building ELI5 summary")
    eli5 = build_eli5(findings, disease_sections, personal_text)

    print(">>> Building dashboard")
    dashboard = build_dashboard(findings, personal_text)

    print(">>> Building ancestry section")
    ancestry = build_ancestry_section(ancestry_data)

    print(">>> Building PRS section")
    prs = build_prs_section(prs_data)

    print(">>> Building critical findings")
    critical = build_critical_findings(personal_text)

    print(">>> Building asthma section")
    asthma = build_asthma_section(personal_text)

    print(">>> Building lifestyle findings")
    lifestyle = build_lifestyle_findings(findings)

    print(">>> Building disease risk")
    disease = build_disease_risk(disease_sections)

    print(">>> Building drug-gene interactions")
    drugs = build_drug_gene(
        findings, pharmgkb_findings,
        disease_sections.get("drug_response", []),
    )

    print(">>> Building nutrition section")
    nutrition = build_nutrition_section(personal_text)

    print(">>> Building monitoring schedule")
    monitoring = build_monitoring(personal_text)

    print(">>> Building protective variants")
    protective = build_protective(disease_sections, personal_text)

    print(">>> Building doctor card")
    doctor_card = build_doctor_card(personal_text)

    print(">>> Building references")
    references = build_references(findings)

    # Assemble HTML
    print("\n>>> Assembling HTML report")
    html = HTML_TEMPLATE.format(
        generated_date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        total_snps=summary.get("total_snps", 0),
        num_findings=len(findings),
        num_pharmgkb=len(pharmgkb_findings),
        eli5_content=eli5,
        dashboard_content=dashboard,
        ancestry_content=ancestry,
        prs_content=prs,
        critical_content=critical,
        asthma_content=asthma,
        lifestyle_content=lifestyle,
        disease_content=disease,
        drugs_content=drugs,
        nutrition_content=nutrition,
        monitoring_content=monitoring,
        protective_content=protective,
        doctor_card_content=doctor_card,
        references_content=references,
    )

    output_path = REPORTS_DIR / "ENHANCED_HEALTH_REPORT.html"
    output_path.write_text(html, encoding="utf-8")

    print(f"\n{'='*60}")
    print(f"Report generated: {output_path}")
    print(f"Size: {len(html):,} characters")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
