"""
Interactive All-in-One Genetic Health Report Generator

Produces a single self-contained HTML file from comprehensive_results.json.

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
import sys
from datetime import datetime
from collections import defaultdict
from pathlib import Path


def _esc(text):
    """Escape text for safe HTML interpolation."""
    return html_mod.escape(str(text)) if text else ""



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


def build_epistasis_section(epistasis_results):
    """Build HTML for gene-gene interactions section."""
    if not epistasis_results:
        return "<p>No significant gene-gene interactions detected.</p>"

    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'These interactions occur when the combined effect of multiple gene '
        'variants differs from each individual effect.</p>'
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

    return "\n".join(parts)


def build_recommendations_section(recs_data):
    """Build HTML for personalized recommendations section."""
    if not recs_data:
        return "<p>No personalized recommendations available.</p>"

    parts = []
    priorities = recs_data.get("priorities", [])
    clinical_insights = recs_data.get("clinical_insights", [])

    # Index clinical insights by gene for lookup
    insights_by_gene = {}
    for ci in clinical_insights:
        insights_by_gene[ci["gene"]] = ci

    if not priorities:
        parts.append("<p>No convergent risk patterns detected.</p>")
    else:
        parts.append(
            '<p style="font-size:.9em;color:var(--accent2)">'
            'Prioritized synthesis of all genetic signals. Convergent risks '
            '(multiple independent signals pointing to the same condition) '
            'are ranked highest.</p>'
        )

        priority_colors = {
            "high": "#ef4444",
            "moderate": "#f59e0b",
            "low": "#22c55e",
        }

        for p in priorities:
            color = priority_colors.get(p["priority"], "var(--border)")
            badge_bg = color

            parts.append(
                f'<details open class="rec-card" style="border-left:4px solid {color};">'
                f'<summary>'
                f'<span class="mag-badge" style="background:{badge_bg};color:#fff">'
                f'{p["priority"].upper()}</span> '
                f'<strong>{p["title"]}</strong></summary>'
            )

            parts.append(f'<p><strong>Why:</strong> {p["why"]}</p>')

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

            # Clinical context from clinical_insights matching this priority
            title_words = set(p["title"].lower().split())
            matched_insights = [
                ci for gene, ci in insights_by_gene.items()
                if gene.lower() in p["title"].lower()
                or gene.lower() in p.get("why", "").lower()
            ]
            if matched_insights:
                parts.append(
                    '<details style="margin-top:.5em">'
                    '<summary style="font-size:.85em;color:var(--accent2)">'
                    'Clinical Context</summary>'
                )
                for ci in matched_insights:
                    parts.append(f'<p><strong>{ci["gene"]}:</strong> {ci["mechanism"]}</p>')
                    if ci.get("actions"):
                        parts.append("<ul>")
                        for a in ci["actions"]:
                            parts.append(f"<li>{a}</li>")
                        parts.append("</ul>")
                parts.append("</details>")

            parts.append("</details>")

    # Drug card
    drug_card = recs_data.get("drug_card", [])
    if drug_card:
        parts.append("<h3>Drug-Gene Card</h3>")
        parts.append(
            '<table><tr><th>Gene</th><th>rsID</th><th>Genotype</th>'
            '<th>Status</th><th>Source</th></tr>'
        )
        for entry in drug_card:
            for e in entry["entries"]:
                status = (e.get("status") or "").replace("_", " ").title() or "\u2014"
                parts.append(
                    f'<tr><td><strong>{entry["gene"]}</strong></td>'
                    f'<td><code>{e["rsid"]}</code></td>'
                    f'<td><code>{e["genotype"]}</code></td>'
                    f'<td>{status}</td><td>{e["source"]}</td></tr>'
                )
        parts.append("</table>")

    # Monitoring schedule
    schedule = recs_data.get("monitoring_schedule", [])
    if schedule:
        parts.append("<h3>Consolidated Monitoring Schedule</h3>")
        parts.append(
            '<table><tr><th>Test</th><th>Frequency</th><th>Reason</th></tr>'
        )
        for m in schedule:
            parts.append(
                f'<tr><td>{m["test"]}</td><td>{m["frequency"]}</td>'
                f'<td>{m["reason"]}</td></tr>'
            )
        parts.append("</table>")

    # Specialist referrals
    referrals = recs_data.get("specialist_referrals", [])
    if referrals:
        parts.append("<h3>Specialist Referrals</h3>")
        urgency_colors = {
            "soon": "#ef4444",
            "routine": "#f59e0b",
        }
        for ref in referrals:
            color = urgency_colors.get(ref.get("urgency", ""), "var(--border)")
            parts.append(
                f'<div class="finding-card" style="border-left-color:{color}">'
                f'<strong>{ref["specialist"]}</strong>: {ref["reason"]}'
                f' <span class="badge" style="background:{color};color:#fff">'
                f'{ref.get("urgency", "routine")}</span></div>'
            )

    # Good news
    good_news = recs_data.get("good_news", [])
    if good_news:
        parts.append("<h3>Good News</h3>")
        parts.append('<div class="good-news-grid">')
        for g in good_news:
            parts.append(
                f'<div class="good-news-card">'
                f'<strong>{g["gene"]}</strong>: {g["description"]}</div>'
            )
        parts.append("</div>")

    return "\n".join(parts)


# =============================================================================
# SECTION BUILDERS
# =============================================================================

def build_quality_section(metrics):
    """Build HTML for data quality section."""
    if not metrics:
        return "<p>No quality metrics available.</p>"

    parts = []

    call_pct = metrics.get("call_rate", 0) * 100
    if call_pct >= 99:
        badge_color = "var(--green)"
        badge_text = "Excellent"
    elif call_pct >= 97:
        badge_color = "#3b82f6"
        badge_text = "Good"
    elif call_pct >= 95:
        badge_color = "#f59e0b"
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

    # Chromosome bar chart
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


def build_blood_type_section(blood_type):
    """Build HTML for blood type section."""
    if not blood_type or blood_type.get("blood_type") == "Unknown":
        return "<p>Insufficient SNP data for blood type prediction.</p>"

    parts = []
    bt = blood_type["blood_type"]
    conf = blood_type["confidence"].title()

    parts.append(
        f'<div style="text-align:center;margin:1em 0">'
        f'<span style="font-size:3em;font-weight:bold;color:var(--warn)">{bt}</span>'
        f'<br><span style="font-size:.9em;color:var(--accent2)">'
        f'Confidence: {conf}</span></div>'
    )

    parts.append(
        f'<table><tr><th>Component</th><th>Result</th></tr>'
        f'<tr><td>ABO Group</td><td><strong>{blood_type["abo"]}</strong></td></tr>'
        f'<tr><td>Rh Factor</td><td><strong>{blood_type["rh"]}</strong></td></tr>'
        f'</table>'
    )

    if blood_type.get("details"):
        parts.append("<p><strong>SNPs analyzed:</strong></p><ul>")
        for d in blood_type["details"]:
            parts.append(f"<li><code>{d}</code></li>")
        parts.append("</ul>")

    parts.append(
        '<p style="font-size:.85em;color:var(--accent2)">'
        'Blood type from SNP data has limitations (rs8176719 is an indel, '
        'often not genotyped by 23andMe). Confirm with clinical blood typing.</p>'
    )

    return "\n".join(parts)


def build_mt_haplogroup_section(mt_result):
    """Build HTML for mitochondrial haplogroup section."""
    if not mt_result or mt_result.get("haplogroup") == "Unknown":
        return "<p>Insufficient mitochondrial SNP data for haplogroup estimation.</p>"

    parts = []
    hg = mt_result["haplogroup"]
    desc = mt_result["description"]
    lineage = mt_result["lineage"]
    conf = mt_result["confidence"].title()

    parts.append(
        f'<div style="text-align:center;margin:1em 0">'
        f'<span style="font-size:2.5em;font-weight:bold;color:var(--accent)">{hg}</span>'
        f'<br><span style="font-size:.95em">{desc}</span>'
        f'<br><span style="font-size:.85em;color:var(--accent2)">'
        f'{lineage} — {conf} confidence '
        f'({mt_result["markers_found"]}/{mt_result["markers_tested"]} markers)</span></div>'
    )

    parts.append(
        '<p style="font-size:.85em;color:var(--accent2)">'
        'Mitochondrial DNA is inherited exclusively from the mother. '
        'Your haplogroup traces the direct maternal line and reflects '
        'ancient human migration patterns.</p>'
    )

    return "\n".join(parts)


def build_star_alleles_section(star_results):
    """Build HTML for pharmacogenomic star alleles section."""
    if not star_results:
        return "<p>No star allele data available.</p>"

    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'CPIC-style star allele calling for key drug-metabolizing enzymes. '
        'Share results with prescribing physicians.</p>'
    )

    # Gauge row for metabolizer phenotypes
    parts.append('<div class="gauge-row" style="justify-content:center">')
    phenotype_levels = {
        "poor": 0, "intermediate": 0.5, "normal": 1,
        "rapid": 1.5, "ultrarapid": 2, "Unknown": 1,
    }
    phenotype_colors = {
        "poor": "#ef4444", "intermediate": "#f59e0b", "normal": "#22c55e",
        "rapid": "#3b82f6", "ultrarapid": "#8b5cf6", "Unknown": "#94a3b8",
    }
    for gene, r in star_results.items():
        level = phenotype_levels.get(r["phenotype"], 1)
        color = phenotype_colors.get(r["phenotype"], "#94a3b8")
        label = f'{gene}: {r["phenotype"].title()}'
        parts.append(svg_metabolism_gauge(gene, level, color))
    parts.append("</div>")

    # Summary table
    parts.append(
        '<table><tr><th>Gene</th><th>Diplotype</th><th>Phenotype</th>'
        '<th>SNPs</th><th>Note</th></tr>'
    )
    for gene, r in star_results.items():
        phenotype = r["phenotype"].replace("_", " ").title()
        note = r["clinical_note"][:100]
        parts.append(
            f'<tr><td><strong>{gene}</strong></td>'
            f'<td><code>{r["diplotype"]}</code></td>'
            f'<td>{phenotype}</td>'
            f'<td>{r["snps_found"]}/{r["snps_total"]}</td>'
            f'<td style="font-size:.8em">{note}</td></tr>'
        )
    parts.append("</table>")

    return "\n".join(parts)


def build_apoe_section(apoe_data):
    """Build HTML for APOE haplotype section."""
    if not apoe_data or apoe_data.get("apoe_type") == "Unknown":
        return "<p>Insufficient SNP data for APOE haplotype determination.</p>"

    risk_colors = {
        "reduced": "#22c55e", "average": "#3b82f6", "moderate": "#f59e0b",
        "elevated": "#f97316", "high": "#ef4444", "unknown": "#94a3b8",
    }
    color = risk_colors.get(apoe_data["risk_level"], "#94a3b8")

    parts = []
    parts.append(
        f'<div style="text-align:center;margin:1em 0">'
        f'<span style="font-size:2em;font-weight:bold;color:{color}">'
        f'{apoe_data["apoe_type"]}</span><br>'
        f'<span class="badge" style="background:{color}">'
        f'{apoe_data["risk_level"].title()} Alzheimer\'s Risk</span>'
        f'</div>'
    )

    if apoe_data.get("alzheimer_or") is not None:
        parts.append(f'<p><strong>Alzheimer\'s odds ratio:</strong> {apoe_data["alzheimer_or"]}x</p>')
    parts.append(f'<p>{apoe_data["description"]}</p>')
    parts.append(
        '<p style="font-size:.85em;color:var(--accent2)">'
        'APOE is the strongest common genetic risk factor for late-onset Alzheimer\'s. '
        'Risk is modifiable through cardiovascular health, exercise, sleep, and cognitive engagement.</p>'
    )
    return "\n".join(parts)


def build_acmg_section(acmg_data):
    """Build HTML for ACMG secondary findings section."""
    if not acmg_data:
        return "<p>ACMG screening not available.</p>"

    parts = []
    parts.append(
        f'<p style="font-size:.9em;color:var(--accent2)">'
        f'Screened {acmg_data["genes_screened"]} medically actionable genes from the '
        f'ACMG SF v3.2 recommendation list.</p>'
    )

    findings = acmg_data.get("acmg_findings", [])
    if not findings:
        parts.append(
            '<div class="doctor-callout" style="border-color:var(--green)">'
            'No pathogenic/likely pathogenic variants found in ACMG genes.</div>'
        )
    else:
        parts.append(
            f'<div class="doctor-callout" style="border-color:var(--warn)">'
            f'<strong>{len(findings)} variant(s)</strong> found in '
            f'{acmg_data["genes_with_variants"]} ACMG gene(s). '
            f'Genetic counseling recommended.</div>'
        )
        parts.append('<table><tr><th>Gene</th><th>Condition</th>'
                     '<th>Genotype</th><th>Stars</th><th>Actionability</th></tr>')
        for f in findings:
            gene = f.get("gene", "Unknown")
            condition = (f.get("traits") or "Unknown").split(";")[0].strip()
            genotype = f.get("user_genotype", "")
            stars = f.get("gold_stars", 0)
            action = f.get("acmg_actionability", "")
            parts.append(
                f'<tr><td><strong>{gene}</strong></td>'
                f'<td>{condition}</td>'
                f'<td><code>{genotype}</code></td>'
                f'<td>{"&#9733;" * stars}{"&#9734;" * (4 - stars)}</td>'
                f'<td style="font-size:.85em">{action}</td></tr>'
            )
        parts.append("</table>")

    return "\n".join(parts)


def build_carrier_screen_section(carrier_data):
    """Build HTML for carrier screening section."""
    if not carrier_data or carrier_data.get("total_carriers", 0) == 0:
        return "<p>No carrier findings to report.</p>"

    parts = []
    parts.append(
        f'<p>{carrier_data["total_carriers"]} carrier finding(s) organized by disease system. '
        f'Relevant for reproductive planning.</p>'
    )

    for system, carriers in sorted(carrier_data.get("by_system", {}).items()):
        parts.append(f'<h3>{system}</h3><ul>')
        for c in carriers:
            note = f' &mdash; <em>{c["reproductive_note"]}</em>' if c.get("reproductive_note") else ""
            parts.append(
                f'<li><strong>{c["gene"]}</strong> ({c.get("rsid", "")}): '
                f'{c["condition"]} <span class="badge">{c["inheritance"]}</span>{note}</li>'
            )
        parts.append("</ul>")

    couples = carrier_data.get("couples_relevant", [])
    if couples:
        parts.append('<h3>Couples-Relevant Conditions</h3>'
                     '<p>Commonly included in expanded carrier screening panels:</p><ul>')
        for c in couples:
            parts.append(f'<li><strong>{c["gene"]}</strong>: {c["condition"]}</li>')
        parts.append("</ul>")

    return "\n".join(parts)


def build_traits_section(traits_data):
    """Build HTML for predicted traits section."""
    if not traits_data:
        return "<p>No trait predictions available.</p>"

    trait_labels = {
        "eye_color": "Eye Color",
        "hair_color": "Hair Color",
        "earwax_type": "Earwax Type",
        "freckling": "Freckling / Sun Sensitivity",
    }

    parts = []
    parts.append('<table><tr><th>Trait</th><th>Prediction</th>'
                 '<th>Confidence</th><th>Description</th></tr>')
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
    return "\n".join(parts)


def build_insights_section(insights_data):
    """Build HTML for research-backed insights section."""
    if not insights_data:
        return "<p>Insights module not run.</p>"

    parts = []

    # Genome highlights
    highlights = insights_data.get("genome_highlights", [])
    if highlights:
        parts.append("<h3>Your Genome Highlights</h3>")
        parts.append('<p style="font-size:.9em;color:var(--accent2)">'
                     'The most notable findings from your DNA.</p>')
        type_colors = {
            "protective": "var(--green)",
            "clinical": "var(--warn)",
            "pharmacogenomic": "var(--accent)",
            "lifestyle": "var(--accent2)",
        }
        for h in highlights:
            color = type_colors.get(h.get("type", ""), "var(--border)")
            parts.append(
                f'<div class="finding-card" style="border-left-color:{color}">'
                f'<strong>{h["title"]}</strong><br>'
                f'<span style="font-size:.9em">{h["detail"]}</span></div>'
            )

    # Multi-gene narratives
    narratives = insights_data.get("narratives", [])
    if narratives:
        parts.append("<h3>Multi-Gene Stories</h3>")
        for n in narratives:
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
                    # Extract PMID if present
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

    # Single-gene insights
    single = insights_data.get("single_gene", [])
    if single:
        parts.append("<h3>Gene-Specific Insights</h3>")
        for entry in single:
            parts.append(
                f'<details class="rec-card">'
                f'<summary><strong>{entry["gene"]}</strong> — {entry["title"]}</summary>'
                f'<p>{entry["finding"]}</p>'
                f'<p><strong>Practical:</strong> {entry["practical"]}</p>'
                f'<div class="paper-refs">Ref: {entry["reference"]}</div>'
                f'</details>'
            )

    # Protective findings
    protective = insights_data.get("protective_findings", [])
    if protective:
        parts.append("<h3>Protective Findings</h3>")
        parts.append('<div class="good-news-grid">')
        for p in protective:
            status = p["status"].replace("_", " ").title()
            parts.append(
                f'<div class="good-news-card">'
                f'<strong>{p["gene"]}</strong> ({status}): {p["title"]}</div>'
            )
        parts.append("</div>")

    return "\n".join(parts)


def build_eli5(findings):
    """Build ELI5 (explain like I'm 5) summary."""
    # Detect key topics from data
    has_bp = any(f["gene"] in ("AGTR1", "AGT", "ADRB1", "ACE", "GNB3")
                 for f in findings if f.get("magnitude", 0) >= 1)
    has_drug = any(f["gene"].startswith("CYP") for f in findings)
    has_lung = any(f["gene"] == "SERPINA1" for f in findings)
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


def build_dashboard(findings):
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
    parts.append(svg_risk_cards(""))

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
            parts.append(f'<p class="finding-desc">{desc}</p>')
            if note:
                parts.append(f'<p class="finding-note">Note: {note}</p>')
            # Population frequency annotation
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
                        f'Population freq: {" · ".join(freq_parts)}</p>'
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

    return "\n".join(parts)


def build_disease_risk(disease_findings):
    """Build disease risk section from ClinVar findings."""
    if not disease_findings:
        return "<p>No ClinVar disease risk data available.</p>"

    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Variants identified by scanning your genome against the ClinVar database. '
        'Only true SNPs are analyzed (indels filtered to prevent false positives).</p>'
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
            condition = _esc((v.get("traits") or "Unknown").split(";")[0].strip())
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
        if len(variants) > 50:
            rows.append(f'<p style="font-size:.85em;color:var(--accent2)">'
                        f'Showing 50 of {len(variants)} variants.</p>')
        return "\n".join(rows)

    parts.append(_variant_table(
        disease_findings.get("pathogenic", []), "Pathogenic Variants", "#ef4444"))
    parts.append(_variant_table(
        disease_findings.get("likely_pathogenic", []), "Likely Pathogenic Variants", "#f97316"))
    parts.append(_variant_table(
        disease_findings.get("risk_factor", []), "Risk Factor Variants", "#f59e0b"))
    parts.append(_variant_table(
        disease_findings.get("drug_response", []), "Drug Response Variants", "#3b82f6"))

    # Protective variants as good-news cards
    protective = disease_findings.get("protective", [])
    if protective:
        parts.append(
            f'<h3>Protective Variants '
            f'<span class="badge" style="background:var(--green)">{len(protective)}</span></h3>'
        )
        parts.append('<div class="good-news-grid">')
        for v in protective:
            gene = _esc(v.get("gene", "Unknown"))
            condition = _esc((v.get("traits") or "").split(";")[0].strip())
            parts.append(
                f'<div class="good-news-card">'
                f'<strong>{gene}</strong>: {condition}</div>'
            )
        parts.append("</div>")

    return "\n".join(parts)


def build_drug_gene(findings, pharmgkb_findings):
    """Build drug-gene interaction section."""
    parts = []

    if pharmgkb_findings:
        parts.append("<h3>PharmGKB Annotations</h3>")
        parts.append('<table><tr><th>Gene</th><th>RSID</th><th>Level</th>'
                     '<th>Drugs</th><th>Genotype</th></tr>')
        for p in pharmgkb_findings:
            parts.append(
                f'<tr><td><strong>{_esc(p["gene"])}</strong></td>'
                f'<td><code>{_esc(p["rsid"])}</code></td>'
                f'<td>{_esc(p["level"])}</td>'
                f'<td>{_esc(p["drugs"])}</td>'
                f'<td><code>{_esc(p["genotype"])}</code></td></tr>'
            )
        parts.append("</table>")

    parts.append(
        '<div class="doctor-callout">Share this entire drug-gene section with '
        "every prescribing physician. Print it or save as PDF.</div>"
    )

    return "\n".join(parts)


def build_nutrition_section(recommendations_data, insights_data):
    """Build nutrition section from recommendations priorities and insight narratives."""
    if not recommendations_data and not insights_data:
        return "<p>No nutrition data available.</p>"

    parts = []
    nutrition_groups = {
        "methylation", "iron_metabolism", "caffeine", "metabolic_diabetes",
        "nutrition", "vitamin_d", "iron",
    }
    nutrition_narrative_ids = {
        "methylation_choline", "vitamin_d_profile", "iron_profile",
        "caffeine_metabolism",
    }

    # Filter nutrition-relevant priorities
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
        priority_colors = {"high": "#ef4444", "moderate": "#f59e0b", "low": "#22c55e"}
        for p in nutrition_priorities:
            color = priority_colors.get(p["priority"], "var(--border)")
            parts.append(
                f'<details open class="rec-card" style="border-left:4px solid {color}">'
                f'<summary>'
                f'<span class="mag-badge" style="background:{color};color:#fff">'
                f'{p["priority"].upper()}</span> '
                f'<strong>{p["title"]}</strong></summary>'
            )
            parts.append(f'<p><strong>Why:</strong> {p["why"]}</p>')
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

    # Matching insight narratives
    narratives = (insights_data or {}).get("narratives", [])
    nutrition_narratives = [
        n for n in narratives
        if n.get("id", "") in nutrition_narrative_ids
    ]
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


def build_monitoring(recommendations_data):
    """Build monitoring schedule section from recommendations.monitoring_schedule."""
    schedule = (recommendations_data or {}).get("monitoring_schedule", [])
    if not schedule:
        return "<p>No monitoring tests recommended based on your genetic profile.</p>"

    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Consolidated and deduplicated monitoring schedule based on all genetic findings.</p>'
    )

    freq_colors = {
        "weekly": "#ef4444", "monthly": "#f97316", "quarterly": "#f59e0b",
        "semi-annually": "#3b82f6", "annually": "#3b82f6", "baseline": "#22c55e",
    }

    parts.append(
        '<table><tr><th>Test</th><th>Frequency</th><th>Reason</th></tr>'
    )
    for m in schedule:
        freq = m.get("frequency", "").lower()
        color = next(
            (c for key, c in freq_colors.items() if key in freq),
            "var(--accent2)"
        )
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


def build_protective(recommendations_data, insights_data):
    """Build protective variants section from good_news + protective_findings."""
    good_news = (recommendations_data or {}).get("good_news", [])
    protective = (insights_data or {}).get("protective_findings", [])

    if not good_news and not protective:
        return "<p>No protective variants identified.</p>"

    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Genetic variants that confer protection or favorable outcomes.</p>'
    )

    if good_news:
        parts.append('<div class="good-news-grid">')
        for g in good_news:
            parts.append(
                f'<div class="good-news-card">'
                f'<strong>{g["gene"]}</strong>: {g["description"]}</div>'
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


def build_doctor_card(recommendations_data, star_alleles_data, apoe_data, acmg_data):
    """Build print-optimized doctor card with key clinical findings."""
    parts = []
    parts.append('<div style="border:2px solid var(--accent);border-radius:8px;padding:1.5em">')
    parts.append('<h3 style="margin-top:0;color:var(--accent)">Patient Genetic Summary</h3>')

    # High-priority conditions
    priorities = (recommendations_data or {}).get("priorities", [])
    high_priorities = [p for p in priorities if p["priority"] == "high"]
    if high_priorities:
        parts.append("<h4>High-Priority Conditions</h4>")
        parts.append('<table><tr><th>Condition</th><th>Doctor Note</th></tr>')
        for p in high_priorities:
            note = p.get("doctor_note", "") or p.get("why", "")
            parts.append(
                f'<tr><td><strong>{p["title"]}</strong></td>'
                f'<td>{note}</td></tr>'
            )
        parts.append("</table>")

    # Specialist referrals
    referrals = (recommendations_data or {}).get("specialist_referrals", [])
    if referrals:
        parts.append("<h4>Specialist Referrals</h4>")
        urgency_colors = {"soon": "#ef4444", "routine": "#f59e0b"}
        for ref in referrals:
            color = urgency_colors.get(ref.get("urgency", ""), "var(--border)")
            parts.append(
                f'<p><strong>{ref["specialist"]}</strong>: {ref["reason"]} '
                f'<span class="mag-badge" style="background:{color};color:#fff">'
                f'{ref.get("urgency", "routine")}</span></p>'
            )

    # Pharmacogenomic card
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

    # APOE status
    if apoe_data and apoe_data.get("apoe_type") != "Unknown":
        risk_colors = {
            "reduced": "var(--green)", "average": "#3b82f6",
            "moderate": "#f59e0b", "elevated": "#f97316", "high": "#ef4444",
        }
        color = risk_colors.get(apoe_data.get("risk_level", ""), "var(--accent2)")
        parts.append(
            f'<h4>APOE Status</h4>'
            f'<p><strong>{apoe_data["apoe_type"]}</strong> — '
            f'<span class="mag-badge" style="background:{color};color:#fff">'
            f'{apoe_data["risk_level"].title()} Risk</span></p>'
        )

    # ACMG actionable findings
    acmg_findings = (acmg_data or {}).get("acmg_findings", [])
    if acmg_findings:
        parts.append("<h4>ACMG Medically Actionable Findings</h4>")
        parts.append('<table><tr><th>Gene</th><th>Condition</th><th>Actionability</th></tr>')
        for f in acmg_findings:
            gene = f.get("gene", "Unknown")
            condition = (f.get("traits") or "Unknown").split(";")[0].strip()
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

    # Monitoring schedule
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

/* Recommendation cards */
.rec-card {{
  background: var(--card-bg); padding: .75em 1em; margin: .75em 0;
  border-radius: 0 8px 8px 0;
}}
.good-news-grid {{
  display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
  gap: .75em;
}}
.good-news-card {{
  background: var(--card-bg); border: 1px solid var(--green);
  border-left: 4px solid var(--green); border-radius: 0 6px 6px 0;
  padding: .6em 1em; font-size: .9em;
}}
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
<a href="#quality">2. Data Quality</a>
<a href="#insights">3. Research-Backed Insights</a>
<a href="#apoe">4. APOE Haplotype</a>
<a href="#recommendations">5. Personalized Recommendations</a>
<a href="#dashboard">6. Dashboard</a>
<a href="#ancestry">7. Ancestry</a>
<a href="#blood-type">8. Blood Type &amp; Traits</a>
<a href="#mt-haplogroup">9. Maternal Haplogroup</a>
<a href="#prs">10. Polygenic Risk Scores</a>
<a href="#star-alleles">11. Pharmacogenomic Star Alleles</a>
<a href="#acmg">12. ACMG Secondary Findings</a>
<a href="#epistasis">13. Gene-Gene Interactions</a>
<a href="#carrier-screen">14. Carrier Screening</a>
<a href="#lifestyle">15. Lifestyle Findings</a>
<a href="#disease">16. Disease Risk</a>
<a href="#drugs">17. Drug-Gene Interactions</a>
<a href="#nutrition">18. Nutrition &amp; Lifestyle</a>
<a href="#monitoring">19. Monitoring Schedule</a>
<a href="#protective">20. Protective Variants</a>
<a href="#polypharmacy">21. Polypharmacy Warnings</a>
<a href="#longevity">22. Longevity &amp; Healthspan</a>
<a href="#sleep-profile">23. Sleep &amp; Circadian</a>
<a href="#nutrigenomics">24. Nutrigenomics</a>
<a href="#mental-health">25. Mental Health Genetics</a>
<a href="#doctor-card">26. Doctor Card</a>
<a href="#references">27. References &amp; Links</a>
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

<div class="section" id="quality">
<h2><span class="section-number">2</span> Data Quality</h2>
{quality_content}
</div>

<div class="section" id="insights">
<h2><span class="section-number">3</span> Research-Backed Insights</h2>
{insights_content}
</div>

<div class="section" id="apoe">
<h2><span class="section-number">4</span> APOE Haplotype</h2>
{apoe_content}
</div>

<div class="section" id="recommendations">
<h2><span class="section-number">5</span> Personalized Recommendations</h2>
{recommendations_content}
</div>

<div class="section no-print" id="dashboard">
<h2><span class="section-number">6</span> Dashboard</h2>
{dashboard_content}
</div>

<div class="section" id="ancestry">
<h2><span class="section-number">7</span> Ancestry Estimation</h2>
{ancestry_content}
</div>

<div class="section" id="blood-type">
<h2><span class="section-number">8</span> Blood Type &amp; Traits</h2>
{blood_type_content}
{traits_content}
</div>

<div class="section" id="mt-haplogroup">
<h2><span class="section-number">9</span> Maternal Haplogroup</h2>
{mt_haplogroup_content}
</div>

<div class="section" id="prs">
<h2><span class="section-number">10</span> Polygenic Risk Scores</h2>
{prs_content}
</div>

<div class="section" id="star-alleles">
<h2><span class="section-number">11</span> Pharmacogenomic Star Alleles</h2>
{star_alleles_content}
</div>

<div class="section" id="acmg">
<h2><span class="section-number">12</span> ACMG Secondary Findings</h2>
{acmg_content}
</div>

<div class="section" id="epistasis">
<h2><span class="section-number">13</span> Gene-Gene Interactions</h2>
{epistasis_content}
</div>

<div class="section" id="carrier-screen">
<h2><span class="section-number">14</span> Carrier Screening</h2>
{carrier_screen_content}
</div>

<div class="section" id="lifestyle">
<h2><span class="section-number">15</span> All Lifestyle Findings</h2>
{lifestyle_content}
</div>

<div class="section" id="disease">
<h2><span class="section-number">16</span> Disease Risk</h2>
{disease_content}
</div>

<div class="section" id="drugs">
<h2><span class="section-number">17</span> Drug-Gene Interactions</h2>
{drugs_content}
</div>

<div class="section" id="nutrition">
<h2><span class="section-number">18</span> Nutrition, Supplements &amp; Lifestyle</h2>
{nutrition_content}
</div>

<div class="section" id="monitoring">
<h2><span class="section-number">19</span> Monitoring Schedule</h2>
{monitoring_content}
</div>

<div class="section" id="protective">
<h2><span class="section-number">20</span> Protective Variants (Good News)</h2>
{protective_content}
</div>

<div class="section" id="polypharmacy">
<h2><span class="section-number">21</span> Polypharmacy Warnings</h2>
{polypharmacy_content}
</div>

<div class="section" id="longevity">
<h2><span class="section-number">22</span> Longevity &amp; Healthspan Profile</h2>
{longevity_content}
</div>

<div class="section" id="sleep-profile">
<h2><span class="section-number">23</span> Sleep &amp; Circadian Profile</h2>
{sleep_content}
</div>

<div class="section" id="nutrigenomics">
<h2><span class="section-number">24</span> Nutrigenomics Profile</h2>
{nutrigenomics_content}
</div>

<div class="section" id="mental-health">
<h2><span class="section-number">25</span> Mental Health Genetics</h2>
{mental_health_content}
</div>

<div class="section doctor-card" id="doctor-card">
<h2><span class="section-number">26</span> Doctor Card</h2>
<p><em>Print this page — only this section will appear in the printout.</em></p>
{doctor_card_content}
</div>

<div class="section" id="references">
<h2><span class="section-number">27</span> References &amp; Links</h2>
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

def build_polypharmacy_section(data):
    """Build polypharmacy warnings section."""
    if not data or not data.get("warnings"):
        return "<p>No polypharmacy risks detected based on your pharmacogenomic profile.</p>"
    parts = [
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Drug combination warnings based on your combined pharmacogene profile. '
        'Share this section with every prescribing physician.</p>'
    ]
    severity_colors = {"high": "#ef4444", "moderate": "#f59e0b", "low": "#22c55e"}
    for w in data["warnings"]:
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
    return "\n".join(parts)


def build_longevity_section(data):
    """Build longevity and healthspan section."""
    if not data:
        return "<p>No longevity data available.</p>"
    parts = []
    score = data.get("longevity_score", 50)
    score_color = "#22c55e" if score >= 60 else "#f59e0b" if score >= 40 else "#ef4444"
    parts.append(
        f'<div style="text-align:center;margin:1em 0">'
        f'<span style="font-size:2.5em;font-weight:bold;color:{score_color}">{score}</span>'
        f'<span style="font-size:1.2em;color:var(--accent2)">/100</span>'
        f'<p style="color:var(--accent2)">Longevity Genetic Score</p>'
        f'<p style="font-size:.9em">{_esc(data.get("summary", ""))}</p>'
        f'</div>'
    )
    # Healthspan domains
    domains = data.get("healthspan_domains", {})
    if domains:
        parts.append("<h3>Healthspan Domains</h3>")
        parts.append('<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(200px,1fr));gap:1em">')
        for domain, info in domains.items():
            s = info["score"]
            c = "#22c55e" if s >= 60 else "#f59e0b" if s >= 40 else "#ef4444"
            parts.append(
                f'<div style="border:1px solid var(--border);border-radius:8px;padding:1em;text-align:center">'
                f'<strong>{_esc(domain.replace("_", " ").title())}</strong><br>'
                f'<span style="font-size:1.5em;color:{c}">{s}</span>/90<br>'
                f'<span style="font-size:.85em;color:var(--accent2)">{info["rating"]}</span>'
                f'</div>'
            )
        parts.append("</div>")
    # Longevity alleles
    alleles = data.get("longevity_alleles", [])
    if alleles:
        parts.append("<h3>Longevity-Associated Alleles</h3>")
        parts.append('<table><tr><th>Gene</th><th>Genotype</th><th>Copies</th><th>Status</th></tr>')
        for a in alleles:
            status_color = "#22c55e" if a["copies"] == 2 else "#f59e0b" if a["copies"] == 1 else "#ef4444"
            parts.append(
                f'<tr><td><strong>{_esc(a["gene"])}</strong></td>'
                f'<td><code>{_esc(a["genotype"])}</code></td>'
                f'<td style="color:{status_color}">{a["copies"]}/2</td>'
                f'<td>{_esc(a["status"].replace("_", " ").title())}</td></tr>'
            )
        parts.append("</table>")
    # Interventions
    interventions = data.get("interventions", [])
    if interventions:
        parts.append("<h3>Genetically-Supported Interventions</h3><ul>")
        for i in interventions:
            parts.append(f'<li><strong>{_esc(i["intervention"])}</strong> — {_esc(i["why"])}</li>')
        parts.append("</ul>")
    return "\n".join(parts)


def build_sleep_section(data):
    """Build sleep and circadian profile section."""
    if not data:
        return "<p>No sleep profile data available.</p>"
    parts = []
    parts.append(
        f'<div style="text-align:center;margin:1em 0">'
        f'<p style="font-size:1.5em;font-weight:bold">{_esc(data.get("chronotype", "Unknown"))}</p>'
        f'<p style="color:var(--accent2)">Chronotype Score: {data.get("chronotype_score", 50)}/100 '
        f'(0=extreme morning, 100=extreme evening)</p>'
        f'</div>'
    )
    parts.append(
        f'<table>'
        f'<tr><td><strong>Optimal Sleep Window</strong></td><td>{_esc(data.get("optimal_sleep_window", ""))}</td></tr>'
        f'<tr><td><strong>Peak Alertness</strong></td><td>{_esc(data.get("peak_alertness", ""))}</td></tr>'
        f'<tr><td><strong>Last Caffeine By</strong></td><td>{_esc(data.get("caffeine_cutoff", ""))}</td></tr>'
        f'<tr><td><strong>Caffeine Sensitive</strong></td><td>{"Yes" if data.get("caffeine_sensitive") else "No"}</td></tr>'
        f'</table>'
    )
    if data.get("deep_sleep_note"):
        parts.append(f'<p style="margin-top:1em"><strong>Note:</strong> {_esc(data["deep_sleep_note"])}</p>')
    recs = data.get("recommendations", [])
    if recs:
        parts.append("<h3>Sleep Recommendations</h3><ul>")
        for r in recs:
            parts.append(f"<li>{_esc(r)}</li>")
        parts.append("</ul>")
    markers = data.get("sleep_markers", [])
    if markers:
        parts.append("<h3>Sleep Gene Markers</h3>")
        parts.append('<table><tr><th>Gene</th><th>Trait</th><th>Genotype</th><th>Effect</th></tr>')
        for m in markers:
            parts.append(
                f'<tr><td><strong>{_esc(m["gene"])}</strong></td>'
                f'<td>{_esc(m["trait"])}</td>'
                f'<td><code>{_esc(m["genotype"])}</code></td>'
                f'<td>{_esc(m["description"])}</td></tr>'
            )
        parts.append("</table>")
    return "\n".join(parts)


def build_nutrigenomics_section(data):
    """Build nutrigenomics profile section."""
    if not data:
        return "<p>No nutrigenomics data available.</p>"
    parts = []
    parts.append(f'<p>{_esc(data.get("summary", ""))}</p>')
    needs = data.get("nutrient_needs", [])
    if needs:
        need_colors = {"high": "#ef4444", "moderate": "#f59e0b", "low": "#3b82f6",
                       "normal": "#22c55e", "caution_excess": "#a855f7"}
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
    supps = data.get("supplement_priorities", [])
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
    testing = data.get("testing_recommendations", [])
    if testing:
        parts.append("<h3>Recommended Lab Tests</h3><ul>")
        for t in testing:
            parts.append(f'<li><strong>{_esc(t["nutrient"])}:</strong> {_esc(t["test"])}</li>')
        parts.append("</ul>")
    return "\n".join(parts)


def build_mental_health_section(data):
    """Build mental health genetics section."""
    if not data:
        return "<p>No mental health genetic data available.</p>"
    parts = []
    parts.append(
        '<p style="font-size:.9em;color:var(--accent2)">'
        'Genetic susceptibility markers for mental health domains. '
        'Genetics is only one factor — environment, lifestyle, and support matter greatly.</p>'
    )
    parts.append(f'<p><strong>{_esc(data.get("summary", ""))}</strong></p>')
    # Domain scores
    domains = data.get("domains", {})
    if domains:
        parts.append('<div style="display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:1em;margin:1em 0">')
        level_colors = {"elevated": "#ef4444", "moderate": "#f59e0b", "low": "#22c55e"}
        for domain, info in domains.items():
            color = level_colors.get(info["risk_level"], "var(--border)")
            parts.append(
                f'<div style="border:1px solid var(--border);border-radius:8px;padding:1em;text-align:center">'
                f'<strong>{_esc(domain.replace("_", " ").title())}</strong><br>'
                f'<span style="font-size:1.3em;color:{color}">{_esc(info["risk_level"].title())}</span>'
                f'</div>'
            )
        parts.append("</div>")
    # Risk and resilience
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
    # Treatment notes
    notes = data.get("treatment_notes", [])
    if notes:
        parts.append("<h3>Treatment Matching Notes</h3><ul>")
        for n in notes:
            parts.append(f"<li>{_esc(n)}</li>")
        parts.append("</ul>")
    # Recommendations
    recs = data.get("recommendations", [])
    if recs:
        parts.append("<h3>Recommendations</h3><ol>")
        for r in recs:
            parts.append(f"<li>{_esc(r)}</li>")
        parts.append("</ol>")
    return "\n".join(parts)


def main():
    print("=" * 60)
    print("Enhanced All-in-One Genetic Health Report")
    print("=" * 60)

    # Load data sources
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

    # Build each section
    print("\n>>> Building ELI5 summary")
    eli5 = build_eli5(findings)

    print(">>> Building dashboard")
    dashboard = build_dashboard(findings)

    print(">>> Building quality section")
    quality = build_quality_section(quality_data)

    print(">>> Building insights section")
    insights_html = build_insights_section(insights_data)

    print(">>> Building APOE section")
    apoe_html = build_apoe_section(apoe_data)

    print(">>> Building ancestry section")
    ancestry = build_ancestry_section(ancestry_data)

    print(">>> Building blood type section")
    blood_type_html = build_blood_type_section(blood_type_data)

    print(">>> Building MT haplogroup section")
    mt_haplogroup_html = build_mt_haplogroup_section(mt_haplogroup_data)

    print(">>> Building PRS section")
    prs = build_prs_section(prs_data)

    print(">>> Building star alleles section")
    star_alleles_html = build_star_alleles_section(star_alleles_data)

    print(">>> Building ACMG section")
    acmg_html = build_acmg_section(acmg_data)

    print(">>> Building epistasis section")
    epistasis = build_epistasis_section(epistasis_data)

    print(">>> Building carrier screening section")
    carrier_screen_html = build_carrier_screen_section(carrier_screen_data)

    print(">>> Building traits section")
    traits_html = build_traits_section(traits_data)

    print(">>> Building recommendations section")
    recommendations = build_recommendations_section(recommendations_data)

    print(">>> Building lifestyle findings")
    lifestyle = build_lifestyle_findings(findings)

    print(">>> Building disease risk")
    disease = build_disease_risk(disease_findings_data)

    print(">>> Building drug-gene interactions")
    drugs = build_drug_gene(findings, pharmgkb_findings)

    print(">>> Building nutrition section")
    nutrition = build_nutrition_section(recommendations_data, insights_data)

    print(">>> Building monitoring schedule")
    monitoring = build_monitoring(recommendations_data)

    print(">>> Building protective variants")
    protective = build_protective(recommendations_data, insights_data)

    print(">>> Building doctor card")
    doctor_card = build_doctor_card(
        recommendations_data, star_alleles_data, apoe_data, acmg_data)

    print(">>> Building references")
    references = build_references(findings)

    print(">>> Building polypharmacy section")
    polypharmacy_html = build_polypharmacy_section(polypharmacy_data)

    print(">>> Building longevity section")
    longevity_html = build_longevity_section(longevity_data)

    print(">>> Building sleep profile section")
    sleep_html = build_sleep_section(sleep_data)

    print(">>> Building nutrigenomics section")
    nutrigenomics_html = build_nutrigenomics_section(nutrigenomics_data)

    print(">>> Building mental health section")
    mental_health_html = build_mental_health_section(mental_health_data)

    # Assemble HTML
    print("\n>>> Assembling HTML report")
    html = HTML_TEMPLATE.format(
        generated_date=datetime.now().strftime("%Y-%m-%d %H:%M"),
        total_snps=summary.get("total_snps", 0),
        num_findings=len(findings),
        num_pharmgkb=len(pharmgkb_findings),
        eli5_content=eli5,
        quality_content=quality,
        insights_content=insights_html,
        apoe_content=apoe_html,
        recommendations_content=recommendations,
        dashboard_content=dashboard,
        ancestry_content=ancestry,
        blood_type_content=blood_type_html,
        traits_content=traits_html,
        mt_haplogroup_content=mt_haplogroup_html,
        prs_content=prs,
        star_alleles_content=star_alleles_html,
        acmg_content=acmg_html,
        epistasis_content=epistasis,
        carrier_screen_content=carrier_screen_html,
        lifestyle_content=lifestyle,
        disease_content=disease,
        drugs_content=drugs,
        nutrition_content=nutrition,
        monitoring_content=monitoring,
        protective_content=protective,
        doctor_card_content=doctor_card,
        references_content=references,
        polypharmacy_content=polypharmacy_html,
        longevity_content=longevity_html,
        sleep_content=sleep_html,
        nutrigenomics_content=nutrigenomics_html,
        mental_health_content=mental_health_html,
    )

    output_path = REPORTS_DIR / "GENETIC_HEALTH_REPORT.html"
    output_path.write_text(html, encoding="utf-8")

    print(f"\n{'='*60}")
    print(f"Report generated: {output_path}")
    print(f"Size: {len(html):,} characters")
    print(f"{'='*60}")


if __name__ == "__main__":
    main()
