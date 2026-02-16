"""
Report generators for the genetic health analysis pipeline.

Contains three report generators (each writes .md + .html):
  - generate_exhaustive_genetic_report() -> EXHAUSTIVE_GENETIC_REPORT.{md,html}
  - generate_disease_risk_report() -> EXHAUSTIVE_DISEASE_RISK_REPORT.{md,html}
  - generate_actionable_protocol() -> ACTIONABLE_HEALTH_PROTOCOL_V3.{md,html}
"""

import re
from pathlib import Path
from datetime import datetime
from collections import defaultdict


def _print_step(text):
    print(f"\n>>> {text}")


# =============================================================================
# Markdown -> HTML converter (no external dependencies)
# =============================================================================

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  :root {{ --bg: #fff; --fg: #1a1a2e; --accent: #2563eb; --border: #e2e8f0;
           --code-bg: #f1f5f9; --table-stripe: #f8fafc; --warn: #dc2626; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ --bg: #0f172a; --fg: #e2e8f0; --accent: #60a5fa; --border: #334155;
             --code-bg: #1e293b; --table-stripe: #1e293b; --warn: #f87171; }}
  }}
  *, *::before, *::after {{ box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
          line-height: 1.7; color: var(--fg); background: var(--bg);
          max-width: 54em; margin: 0 auto; padding: 2em 1.5em; }}
  h1 {{ border-bottom: 2px solid var(--accent); padding-bottom: .3em; }}
  h2 {{ border-bottom: 1px solid var(--border); padding-bottom: .2em; margin-top: 2em; }}
  h3 {{ margin-top: 1.5em; }}
  hr {{ border: none; border-top: 1px solid var(--border); margin: 2em 0; }}
  code {{ background: var(--code-bg); padding: .15em .35em; border-radius: 3px;
          font-size: .9em; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
  th, td {{ border: 1px solid var(--border); padding: .5em .75em; text-align: left; }}
  th {{ background: var(--code-bg); font-weight: 600; }}
  tr:nth-child(even) {{ background: var(--table-stripe); }}
  ul, ol {{ padding-left: 1.5em; }}
  li {{ margin: .3em 0; }}
  strong {{ color: var(--accent); }}
  .warn {{ color: var(--warn); font-weight: bold; }}
  @media print {{ body {{ max-width: 100%; font-size: 11pt; }} }}
</style>
</head>
<body>
{body}
</body>
</html>
"""


def _md_to_html(md: str) -> str:
    """Convert report Markdown to HTML. Handles the subset used by our reports."""
    lines = md.split("\n")
    html_parts = []
    i = 0

    while i < len(lines):
        line = lines[i]

        # Blank line
        if not line.strip():
            i += 1
            continue

        # Horizontal rule
        if re.match(r"^-{3,}\s*$", line) or re.match(r"^\*{3,}\s*$", line):
            html_parts.append("<hr>")
            i += 1
            continue

        # Headers
        m = re.match(r"^(#{1,6})\s+(.*)", line)
        if m:
            level = len(m.group(1))
            text = _inline(m.group(2))
            html_parts.append(f"<h{level}>{text}</h{level}>")
            i += 1
            continue

        # Table (starts with |)
        if line.strip().startswith("|"):
            table_lines = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                table_lines.append(lines[i])
                i += 1
            html_parts.append(_table_to_html(table_lines))
            continue

        # Unordered list
        if re.match(r"^\s*[-*]\s+", line):
            items = []
            while i < len(lines) and re.match(r"^\s*[-*]\s+", lines[i]):
                items.append(re.sub(r"^\s*[-*]\s+", "", lines[i]))
                i += 1
            html_parts.append("<ul>" + "".join(f"<li>{_inline(it)}</li>" for it in items) + "</ul>")
            continue

        # Paragraph (collect consecutive non-empty, non-special lines)
        para = []
        while i < len(lines) and lines[i].strip() and not re.match(r"^(#{1,6}\s|[-*]{3,}\s*$|\||\s*[-*]\s+)", lines[i]):
            para.append(lines[i])
            i += 1
        if para:
            html_parts.append(f"<p>{_inline(' '.join(para))}</p>")

    return "\n".join(html_parts)


def _inline(text: str) -> str:
    """Convert inline Markdown: bold, code, italic."""
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", text)
    return text


def _table_to_html(lines: list[str]) -> str:
    """Convert Markdown table lines to HTML table."""
    rows = []
    for line in lines:
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        rows.append(cells)

    # Skip separator row (e.g. |---|---|)
    data_rows = [r for r in rows if not all(re.match(r"^[-:]+$", c) for c in r)]
    if not data_rows:
        return ""

    html = "<table>\n<thead><tr>"
    for cell in data_rows[0]:
        html += f"<th>{_inline(cell)}</th>"
    html += "</tr></thead>\n<tbody>\n"
    for row in data_rows[1:]:
        html += "<tr>" + "".join(f"<td>{_inline(c)}</td>" for c in row) + "</tr>\n"
    html += "</tbody></table>"
    return html


def _write_html(md_content: str, md_path: Path):
    """Write an HTML version alongside the Markdown file."""
    html_path = md_path.with_suffix(".html")
    title = "Genetic Health Report"
    # Extract title from first H1
    m = re.search(r"^#\s+(.+)", md_content, re.MULTILINE)
    if m:
        title = m.group(1)
    body = _md_to_html(md_content)
    html = _HTML_TEMPLATE.format(title=title, body=body)
    with open(html_path, "w") as f:
        f.write(html)
    print(f"    HTML:     {html_path}")


def classify_zygosity(finding):
    """Classify zygosity impact."""
    inheritance = finding['inheritance'].lower() if finding['inheritance'] else ''

    if finding['is_homozygous']:
        return 'AFFECTED', 'Homozygous for variant'
    elif finding['is_heterozygous']:
        if 'recessive' in inheritance:
            return 'CARRIER', 'Heterozygous carrier (recessive)'
        elif 'dominant' in inheritance:
            return 'AFFECTED', 'Heterozygous (dominant)'
        else:
            return 'HETEROZYGOUS', 'Heterozygous (inheritance unclear)'
    return 'UNKNOWN', 'Zygosity unclear'


def generate_exhaustive_genetic_report(results: dict, output_path: Path, subject_name: str = None):
    """Generate the exhaustive lifestyle/health genetic report."""
    _print_step("Generating exhaustive genetic report")

    from generate_exhaustive_report import (
        generate_executive_summary, generate_priority_findings,
        generate_pathway_analysis, generate_full_findings,
        generate_pharmgkb_report, generate_action_summary, generate_disclaimer,
    )

    data = {
        'findings': results['findings'],
        'pharmgkb_findings': results['pharmgkb_findings'],
        'summary': results['summary']
    }

    report_parts = []
    report_parts.append(generate_executive_summary(data))
    report_parts.append(generate_priority_findings(results['findings']))
    report_parts.append(generate_pathway_analysis(results['findings']))
    report_parts.append(generate_full_findings(results['findings']))
    report_parts.append(generate_pharmgkb_report(results['pharmgkb_findings']))
    report_parts.append(generate_action_summary(results['findings']))
    report_parts.append(generate_disclaimer())

    full_report = "\n".join(report_parts)

    if subject_name:
        full_report = full_report.replace(
            "# Exhaustive Genetic Health Report",
            f"# Exhaustive Genetic Health Report\n\n**Subject:** {subject_name}"
        )

    with open(output_path, 'w') as f:
        f.write(full_report)

    print(f"    Markdown: {output_path}")
    _write_html(full_report, output_path)


def generate_disease_risk_report(findings: dict, stats: dict, genome_count: int,
                                  output_path: Path, subject_name: str = None):
    """Generate the exhaustive disease risk report."""
    _print_step("Generating disease risk report")

    now = datetime.now().strftime("%Y-%m-%d %H:%M")

    affected = []
    carriers = []
    het_unknown = []

    for f in findings['pathogenic'] + findings['likely_pathogenic']:
        status, desc = classify_zygosity(f)
        f['zygosity_status'] = status
        f['zygosity_description'] = desc

        if status == 'AFFECTED':
            affected.append(f)
        elif status == 'CARRIER':
            carriers.append(f)
        else:
            het_unknown.append(f)

    for lst in [affected, carriers, het_unknown]:
        lst.sort(key=lambda x: (-x['gold_stars'], x['gene']))
    findings['risk_factor'].sort(key=lambda x: (-x['gold_stars'], x['gene']))
    findings['drug_response'].sort(key=lambda x: (-x['gold_stars'], x['gene']))
    findings['protective'].sort(key=lambda x: (-x['gold_stars'], x['gene']))

    subject_line = f"\n**Subject:** {subject_name}" if subject_name else ""

    report = f"""# Exhaustive Disease Risk Report
{subject_line}
**Generated:** {now}

---

## Executive Summary

### Genome Overview
- **Total SNPs in Raw Data:** {genome_count:,}
- **ClinVar Variants Scanned:** {stats['total_clinvar']:,}

### Clinical Findings Summary

| Category | Count | Description |
|----------|-------|-------------|
| **Pathogenic (Affected)** | {len(affected)} | Homozygous or dominant |
| **Pathogenic (Carrier)** | {len(carriers)} | Heterozygous carrier for recessive |
| **Likely Pathogenic** | {len(het_unknown)} | Heterozygous, inheritance unclear |
| **Risk Factors** | {len(findings['risk_factor'])} | Increased susceptibility |
| **Drug Response** | {len(findings['drug_response'])} | Pharmacogenomic variants |
| **Protective** | {len(findings['protective'])} | Reduced disease risk |

---

"""

    if affected:
        report += "## Pathogenic Variants - Affected Status\n\n"
        for f in affected:
            stars = '*' * f['gold_stars'] + '.' * (4 - f['gold_stars'])
            report += f"""### {f['gene']} - {f['traits'].split(';')[0] if f['traits'] else 'Unknown'}

- **Position:** chr{f['chromosome']}:{f['position']}
- **RSID:** {f['rsid']}
- **Genotype:** `{f['user_genotype']}` ({'Homozygous' if f['is_homozygous'] else 'Heterozygous'})
- **Variant:** {f['ref']} -> {f['alt']}
- **Confidence:** {stars} ({f['gold_stars']}/4)
- **Condition:** {f['traits'] if f['traits'] else 'Not specified'}

---

"""

    if carriers:
        report += "## Carrier Status - Recessive Conditions\n\n"
        report += "**You are a carrier - no personal symptoms expected, but reproductive implications.**\n\n"
        for f in carriers:
            stars = '*' * f['gold_stars'] + '.' * (4 - f['gold_stars'])
            report += f"""### {f['gene']} - {f['traits'].split(';')[0] if f['traits'] else 'Unknown'}

- **RSID:** {f['rsid']}
- **Genotype:** `{f['user_genotype']}` (Carrier)
- **Confidence:** {stars} ({f['gold_stars']}/4)
- **Condition:** {f['traits'] if f['traits'] else 'Not specified'}

---

"""

    if het_unknown:
        report += "## Pathogenic/Likely Pathogenic - Inheritance Unclear\n\n"
        for f in het_unknown:
            stars = '*' * f['gold_stars'] + '.' * (4 - f['gold_stars'])
            report += f"""### {f['gene']} - {f['traits'].split(';')[0] if f['traits'] else 'Unknown'}

- **RSID:** {f['rsid']}
- **Genotype:** `{f['user_genotype']}`
- **Confidence:** {stars} ({f['gold_stars']}/4)
- **Condition:** {f['traits'] if f['traits'] else 'Not specified'}

---

"""

    if findings['risk_factor']:
        report += "## Risk Factor Variants\n\n"
        for f in findings['risk_factor'][:30]:
            report += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` - {f['traits'][:80] if f['traits'] else 'Risk factor'}...\n"
        report += "\n---\n\n"

    if findings['drug_response']:
        report += "## Drug Response Variants\n\n"
        for f in findings['drug_response'][:30]:
            report += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` - {f['traits'][:80] if f['traits'] else 'Drug response'}...\n"
        report += "\n---\n\n"

    if findings['protective']:
        report += "## Protective Variants\n\n"
        for f in findings['protective']:
            report += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` - {f['traits'][:80] if f['traits'] else 'Protective'}...\n"
        report += "\n---\n\n"

    report += """## Disclaimer

This report is for **informational purposes only**. It is NOT a clinical diagnosis.

- Consult a genetic counselor or physician for clinical interpretation
- Variant classifications may change as research progresses
- Carrier status has reproductive implications

---

*Generated using ClinVar database*
"""

    with open(output_path, 'w') as f:
        f.write(report)

    print(f"    Markdown: {output_path}")
    _write_html(report, output_path)


def generate_actionable_protocol(health_results: dict, disease_findings: dict,
                                  output_path: Path, subject_name: str = None):
    """Generate comprehensive actionable health protocol combining ALL sources."""
    _print_step("Generating actionable health protocol (comprehensive)")

    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    subject_line = f"\n**Subject:** {subject_name}" if subject_name else ""

    findings_dict = {f['gene']: f for f in health_results['findings']}

    affected = []
    carriers = []
    het_unknown = []

    if disease_findings:
        for f in disease_findings.get('pathogenic', []) + disease_findings.get('likely_pathogenic', []):
            inheritance = f.get('inheritance', '').lower()
            if f.get('is_homozygous'):
                affected.append(f)
            elif f.get('is_heterozygous'):
                if 'recessive' in inheritance:
                    carriers.append(f)
                elif 'dominant' in inheritance:
                    affected.append(f)
                else:
                    het_unknown.append(f)

    total_lifestyle = len(health_results['findings'])
    total_pharmgkb = len(health_results['pharmgkb_findings'])
    total_risk_factors = len(disease_findings.get('risk_factor', [])) if disease_findings else 0
    total_drug_response = len(disease_findings.get('drug_response', [])) if disease_findings else 0
    total_protective = len(disease_findings.get('protective', [])) if disease_findings else 0

    report = f"""# Actionable Health Protocol (V3)
{subject_line}
**Generated:** {now}

This protocol synthesizes ALL genetic findings into concrete recommendations:
- Lifestyle/health genetics ({total_lifestyle} findings)
- PharmGKB drug interactions ({total_pharmgkb} interactions)
- Pathogenic/likely pathogenic variants ({len(affected)} affected, {len(carriers)} carrier, {len(het_unknown)} unclear)
- Risk factors ({total_risk_factors} variants)
- ClinVar drug response ({total_drug_response} variants)
- Protective variants ({total_protective} variants)

---

## Executive Summary

### High-Impact Lifestyle Findings (Magnitude >= 3)

"""

    high_impact = [f for f in health_results['findings'] if f['magnitude'] >= 3]
    if high_impact:
        for f in high_impact:
            report += f"- **{f['gene']}** ({f['category']}): {f['description']}\n"
    else:
        report += "None detected.\n"

    report += "\n### Pathogenic/Likely Pathogenic Variants\n\n"

    if affected:
        report += "**Affected Status:**\n"
        for f in affected:
            condition = f['traits'].split(';')[0] if f['traits'] else 'Unknown condition'
            stars = f['gold_stars']
            confidence = f"({stars}/4 stars)" if stars > 0 else "(low confidence)"
            report += f"- **{f['gene']}**: {condition} {confidence}\n"
        report += "\n"

    if carriers:
        report += "**Carrier Status (Recessive):**\n"
        for f in carriers:
            condition = f['traits'].split(';')[0] if f['traits'] else 'Unknown condition'
            stars = f['gold_stars']
            confidence = f"({stars}/4 stars)" if stars > 0 else "(low confidence)"
            report += f"- **{f['gene']}**: {condition} {confidence}\n"
        report += "\n"

    if het_unknown:
        report += "**Heterozygous (Inheritance Unclear):**\n"
        for f in het_unknown:
            condition = f['traits'].split(';')[0] if f['traits'] else 'Unknown condition'
            stars = f['gold_stars']
            confidence = f"({stars}/4 stars)" if stars > 0 else "(low confidence)"
            report += f"- **{f['gene']}**: {condition} {confidence}\n"
        report += "\n"

    if not affected and not carriers and not het_unknown:
        report += "None detected.\n\n"

    report += "### Protective Variants\n\n"
    if disease_findings and disease_findings.get('protective'):
        for f in disease_findings['protective']:
            condition = f['traits'].split(';')[0] if f['traits'] else 'Protective effect'
            report += f"- **{f['gene']}**: {condition}\n"
    else:
        report += "None detected.\n"

    report += """

---

## Supplement Recommendations

*Discuss with healthcare provider before starting any supplements*

"""

    supplements = []

    if 'MTHFR' in findings_dict and findings_dict['MTHFR']['magnitude'] >= 2:
        supplements.append({
            'name': 'Methylfolate (L-5-MTHF)',
            'dose': '400-800mcg daily',
            'reason': 'MTHFR variant reduces folic acid conversion',
            'source': 'MTHFR (lifestyle)',
            'notes': 'Avoid synthetic folic acid. Start low if slow COMT.'
        })
        supplements.append({
            'name': 'Methylcobalamin (B12)',
            'dose': '1000mcg sublingual',
            'reason': 'Supports methylation cycle',
            'source': 'MTHFR (lifestyle)',
            'notes': 'Prefer methylcobalamin over cyanocobalamin'
        })

    if 'MTRR' in findings_dict and findings_dict['MTRR']['magnitude'] >= 2:
        if not any('B12' in s['name'] for s in supplements):
            supplements.append({
                'name': 'Methylcobalamin (B12)',
                'dose': '1000-5000mcg sublingual',
                'reason': 'MTRR variant impairs B12 recycling',
                'source': 'MTRR (lifestyle)',
                'notes': 'May need higher doses than typical'
            })

    if 'GC' in findings_dict and findings_dict['GC'].get('status') == 'low':
        supplements.append({
            'name': 'Vitamin D3',
            'dose': '2000-5000 IU daily',
            'reason': 'Genetically low vitamin D binding protein',
            'source': 'GC (lifestyle)',
            'notes': 'Take with fat. Test 25-OH-D after 2-3 months. Target 40-60 ng/mL.'
        })
        supplements.append({
            'name': 'Vitamin K2 (MK-7)',
            'dose': '100-200mcg daily',
            'reason': 'Synergistic with D3 for calcium metabolism',
            'source': 'GC (lifestyle)',
            'notes': 'Optional but recommended with high-dose D3'
        })

    if 'FADS1' in findings_dict and findings_dict['FADS1'].get('status') == 'low_conversion':
        supplements.append({
            'name': 'Fish Oil or Algae Oil (EPA/DHA)',
            'dose': '1-2g EPA+DHA daily',
            'reason': 'Poor conversion from plant omega-3s (ALA)',
            'source': 'FADS1 (lifestyle)',
            'notes': 'Direct marine source required. Flax/chia insufficient.'
        })

    if 'COMT' in findings_dict and findings_dict['COMT'].get('status') == 'slow':
        supplements.append({
            'name': 'Magnesium Glycinate',
            'dose': '300-400mg evening',
            'reason': 'Supports COMT function, calming effect',
            'source': 'COMT (lifestyle)',
            'notes': 'Glycinate form preferred for bioavailability and sleep'
        })

    if 'PEMT' in findings_dict:
        supplements.append({
            'name': 'Choline (Phosphatidylcholine or CDP-Choline)',
            'dose': '250-500mg daily',
            'reason': 'PEMT variant increases dietary choline requirement',
            'source': 'PEMT (lifestyle)',
            'notes': 'Eggs are excellent food source (2 eggs = ~300mg)'
        })

    if 'BCMO1' in findings_dict and findings_dict['BCMO1'].get('status') == 'reduced':
        supplements.append({
            'name': 'Preformed Vitamin A or Cod Liver Oil',
            'dose': '2500-5000 IU (as retinol)',
            'reason': 'Poor conversion from beta-carotene',
            'source': 'BCMO1 (lifestyle)',
            'notes': 'Get from food (liver, eggs) or supplement. Avoid excess.'
        })

    if 'IL6' in findings_dict and findings_dict['IL6'].get('status') == 'high':
        supplements.append({
            'name': 'Omega-3 (EPA/DHA)',
            'dose': '2-3g daily',
            'reason': 'Higher baseline inflammation (IL-6)',
            'source': 'IL6 (lifestyle)',
            'notes': 'Anti-inflammatory. Consider curcumin as well.'
        })

    if supplements:
        report += "| Supplement | Dose | Reason | Notes |\n"
        report += "|------------|------|--------|-------|\n"
        for s in supplements:
            report += f"| {s['name']} | {s['dose']} | {s['reason']} | {s['notes']} |\n"
    else:
        report += "No specific supplements indicated by genetic profile.\n"

    report += """

---

## Dietary Recommendations

"""

    diet_recs = []

    if 'APOA2' in findings_dict and findings_dict['APOA2'].get('status') == 'sensitive':
        diet_recs.append("**Limit saturated fat (<7% calories)**: APOA2 variant links sat fat intake to weight gain. Minimize butter, fatty red meat, full-fat dairy, coconut oil. Prefer olive oil, nuts, avocado.")

    if 'MTHFR' in findings_dict and findings_dict['MTHFR']['magnitude'] >= 2:
        diet_recs.append("**Emphasize folate-rich foods**: Leafy greens, legumes, liver. Avoid folic acid-fortified processed foods when possible (UMFA accumulation risk).")

    if 'IL6' in findings_dict:
        diet_recs.append("**Anti-inflammatory diet**: Omega-3 rich fish, colorful vegetables, minimize processed foods. Sleep deprivation spikes IL-6.")

    if 'MCM6/LCT' in findings_dict and 'intolerant' in findings_dict['MCM6/LCT'].get('status', ''):
        diet_recs.append("**Lactose intolerance**: May tolerate small amounts or fermented dairy (yogurt, aged cheese). Lactase supplements available. Ensure calcium from other sources.")

    if 'HLA-DQA1' in findings_dict:
        diet_recs.append("**Celiac risk (HLA-DQ2.5)**: No preventive gluten-free diet needed. If GI symptoms arise, get celiac antibody testing (tTG-IgA) *while still eating gluten*.")

    caffeine_issues = []
    if 'CYP1A2' in findings_dict and findings_dict['CYP1A2'].get('status') in ['slow', 'intermediate']:
        caffeine_issues.append("slow metabolizer")
    if 'ADORA2A' in findings_dict and findings_dict['ADORA2A'].get('status') == 'anxiety_prone':
        caffeine_issues.append("anxiety-prone")
    if 'COMT' in findings_dict and findings_dict['COMT'].get('status') == 'slow':
        caffeine_issues.append("slow COMT")

    if caffeine_issues:
        diet_recs.append(f"**Caffeine caution** ({', '.join(caffeine_issues)}): Limit to morning only (before 10am). Consider lower doses, green tea (L-theanine), or alternatives.")

    if 'HFE' in findings_dict:
        diet_recs.append("**Iron awareness (HFE carrier)**: Don't supplement iron unless deficiency confirmed. Blood donation helps regulate if ferritin runs high.")

    if diet_recs:
        for rec in diet_recs:
            report += f"- {rec}\n\n"
    else:
        report += "No specific dietary modifications beyond general healthy eating.\n"

    report += """
---

## Lifestyle Recommendations

"""

    lifestyle_recs = []

    if 'COMT' in findings_dict and findings_dict['COMT'].get('status') == 'slow':
        lifestyle_recs.append("**Stress management is critical**: Slow COMT means catecholamines (dopamine, norepinephrine) build up under stress. Daily meditation, breathwork, adequate sleep. Avoid combining multiple stimulants.")

    if 'BDNF' in findings_dict and findings_dict['BDNF']['magnitude'] >= 2:
        lifestyle_recs.append("**Exercise is essential**: BDNF variant reduces activity-dependent brain growth factor. Physical activity is one of the strongest natural BDNF boosters.")

    if 'ACTN3' in findings_dict:
        status = findings_dict['ACTN3'].get('status', '')
        if status == 'endurance':
            lifestyle_recs.append("**Training style (ACTN3 endurance)**: Genetics favor endurance/aerobic training. Can still build strength but may excel at higher volume, aerobic work.")
        elif status == 'power':
            lifestyle_recs.append("**Training style (ACTN3 power)**: Genetics favor explosive/strength training. May recover faster from power-based work.")
        else:
            lifestyle_recs.append("**Training style (ACTN3 mixed)**: Versatile profile - respond well to both power and endurance training.")

    if 'ARNTL' in findings_dict:
        lifestyle_recs.append("**Circadian rhythm support (ARNTL)**: May have weaker internal clock. Strong morning light exposure, consistent sleep/wake times even weekends, blue light reduction in evening.")

    bp_genes = ['AGTR1', 'ACE', 'AGT', 'GNB3']
    bp_findings = [findings_dict[g] for g in bp_genes if g in findings_dict]
    if len(bp_findings) >= 2:
        lifestyle_recs.append("**Blood pressure focus**: Multiple BP-related variants. Regular monitoring, sodium restriction, DASH diet pattern, 150+ min/week aerobic exercise.")

    if 'MC1R' in findings_dict:
        lifestyle_recs.append("**Sun protection (MC1R)**: Accelerated skin aging variant. Daily SPF 30+, topical retinoids, antioxidant serums. Avoid excessive sun exposure.")

    if lifestyle_recs:
        for rec in lifestyle_recs:
            report += f"- {rec}\n\n"
    else:
        report += "Standard healthy lifestyle recommendations apply.\n"

    report += """
---

## Monitoring Recommendations

"""

    monitoring = []

    if 'MTHFR' in findings_dict and findings_dict['MTHFR']['magnitude'] >= 2:
        monitoring.append("**Homocysteine**: Annually. Target <10 \u03bcmol/L. MTHFR variant affects metabolism.")

    if 'MTRR' in findings_dict and findings_dict['MTRR']['magnitude'] >= 2:
        monitoring.append("**B12 + Methylmalonic acid (MMA)**: For functional B12 status. MTRR affects recycling.")

    if 'GC' in findings_dict:
        monitoring.append("**25-OH Vitamin D**: After 2-3 months supplementation, then annually. Target 40-60 ng/mL.")

    if any(g in findings_dict for g in ['AGTR1', 'ACE', 'AGT', 'GNB3']):
        monitoring.append("**Blood pressure**: Home monitoring recommended. Multiple BP-related variants.")

    if 'HFE' in findings_dict:
        monitoring.append("**Ferritin/iron panel**: Every 1-2 years. HFE carrier status.")

    if 'TCF7L2' in findings_dict and findings_dict['TCF7L2']['magnitude'] >= 2:
        monitoring.append("**Fasting glucose or HbA1c**: Annually. TCF7L2 diabetes risk variant.")

    if disease_findings:
        risk_conditions = set()
        for f in disease_findings.get('risk_factor', []):
            traits = f.get('traits', '').lower()
            if 'macular degeneration' in traits:
                risk_conditions.add('macular_degeneration')
            if 'diabetes' in traits:
                risk_conditions.add('diabetes')
            if 'hypertension' in traits:
                risk_conditions.add('hypertension')
            if 'thrombosis' in traits or 'thromboembolism' in traits:
                risk_conditions.add('thrombosis')

        if 'macular_degeneration' in risk_conditions:
            monitoring.append("**Eye exams**: Regular ophthalmology. Multiple age-related macular degeneration risk variants (CFH, C3, ERCC6).")
        if 'diabetes' in risk_conditions and 'TCF7L2' not in findings_dict:
            monitoring.append("**Glucose monitoring**: Multiple diabetes susceptibility variants detected.")
        if 'thrombosis' in risk_conditions:
            monitoring.append("**Clotting awareness**: Risk variants for venous thrombosis (F13B, FGA). Stay hydrated, move on long flights, know DVT symptoms.")

    if monitoring:
        for m in monitoring:
            report += f"- {m}\n"
    else:
        report += "Standard health monitoring appropriate for age.\n"

    report += """

---

## Drug-Gene Interactions

**Share this section with prescribing physicians.**

### PharmGKB Level 1 (Clinical Guidelines Exist)

"""

    level_1 = [f for f in health_results['pharmgkb_findings'] if f['level'] in ['1A', '1B']]
    if level_1:
        report += "| Gene | Level | Drugs | Your Genotype |\n"
        report += "|------|-------|-------|---------------|\n"
        for f in level_1:
            drugs = f['drugs'][:50] + '...' if len(f['drugs']) > 50 else f['drugs']
            report += f"| {f['gene']} | {f['level']} | {drugs} | `{f['genotype']}` |\n"
    else:
        report += "None detected.\n"

    report += "\n### PharmGKB Level 2 (Moderate Evidence)\n\n"

    level_2 = [f for f in health_results['pharmgkb_findings'] if f['level'] in ['2A', '2B']]
    if level_2:
        report += "| Gene | Level | Drugs | Your Genotype |\n"
        report += "|------|-------|-------|---------------|\n"
        for f in level_2[:15]:
            drugs = f['drugs'][:50] + '...' if len(f['drugs']) > 50 else f['drugs']
            report += f"| {f['gene']} | {f['level']} | {drugs} | `{f['genotype']}` |\n"
        if len(level_2) > 15:
            report += f"\n*...and {len(level_2) - 15} more Level 2 interactions*\n"
    else:
        report += "None detected.\n"

    report += "\n### ClinVar Drug Response Variants\n\n"

    if disease_findings and disease_findings.get('drug_response'):
        drug_resp = disease_findings['drug_response']
        report += "| Gene | RSID | Genotype | Drug/Response |\n"
        report += "|------|------|----------|---------------|\n"
        for f in drug_resp[:20]:
            traits = f['traits'][:60] + '...' if len(f['traits']) > 60 else f['traits']
            gene = f['gene'] if f['gene'] else '\u2014'
            report += f"| {gene} | {f['rsid']} | `{f['user_genotype']}` | {traits} |\n"
        if len(drug_resp) > 20:
            report += f"\n*...and {len(drug_resp) - 20} more drug response variants*\n"
    else:
        report += "None detected.\n"

    report += """

---

## Carrier Status Notes

"""

    carrier_notes = {
        'CFTR': """**Cystic Fibrosis Carrier (CFTR)**:
- CF carriers may have ~10% reduced lung function (FEV1)
- Increased risk of pancreatitis (2-3x general population)
- Higher prevalence of chronic sinusitis
- Possible male fertility effects (CBAVD spectrum)
- **Recommendation**: Baseline pulmonary function test, avoid smoking, genetic counseling if planning pregnancy
""",
        'HBB': """**Sickle Cell Trait Carrier (HBB)**:
- Generally asymptomatic under normal conditions
- Possible complications at extreme altitude or severe dehydration
- Malaria resistance (evolutionary advantage)
- **Recommendation**: Stay hydrated during intense exercise; inform physicians before surgery
""",
        'GBA': """**Gaucher Disease Carrier (GBA)**:
- Carriers have increased Parkinson's disease risk (5-8x)
- No Gaucher disease symptoms
- **Recommendation**: Awareness of early Parkinson's symptoms; inform neurologist
""",
        'SERPINA1': """**Alpha-1 Antitrypsin Carrier (SERPINA1)**:
- Carriers (MZ) have ~60% normal AAT levels
- Mildly increased risk of COPD, especially if smoking
- **Recommendation**: Absolutely avoid smoking; baseline liver function; consider AAT level testing
"""
    }

    found_carriers = False
    all_carrier_genes = [f['gene'].upper() for f in carriers + het_unknown if f.get('gene')]

    for gene, note in carrier_notes.items():
        if gene in all_carrier_genes:
            report += note + "\n"
            found_carriers = True

    cftr_finding = next((f for f in het_unknown if f.get('gene', '').upper() == 'CFTR'), None)
    if cftr_finding and 'CFTR' not in all_carrier_genes:
        report += carrier_notes['CFTR'] + "\n"
        found_carriers = True

    if not found_carriers:
        if carriers or het_unknown:
            report += "Carrier status detected but no specific phenotype notes available for these genes. General recommendation: genetic counseling if planning pregnancy.\n"
        else:
            report += "No carrier status detected.\n"

    report += """

---

## Risk Factor Summary

*These variants indicate increased susceptibility, not certainty of disease.*

"""

    if disease_findings and disease_findings.get('risk_factor'):
        conditions = defaultdict(list)
        for f in disease_findings['risk_factor']:
            traits = f.get('traits', '').lower()
            gene = f.get('gene', 'Unknown')

            if 'hypertension' in traits:
                conditions['Hypertension'].append(gene)
            elif 'diabetes' in traits:
                conditions['Diabetes'].append(gene)
            elif 'macular degeneration' in traits:
                conditions['Macular Degeneration'].append(gene)
            elif 'thrombosis' in traits or 'thromboembolism' in traits:
                conditions['Thrombosis/Clotting'].append(gene)
            elif 'obesity' in traits:
                conditions['Obesity'].append(gene)
            elif 'cancer' in traits or 'carcinoma' in traits:
                conditions['Cancer Risk'].append(gene)
            elif 'inflammatory bowel' in traits or 'crohn' in traits:
                conditions['Inflammatory Bowel Disease'].append(gene)

        if conditions:
            report += "| Condition | Genes Involved |\n"
            report += "|-----------|----------------|\n"
            for condition, genes in sorted(conditions.items()):
                unique_genes = list(set(g for g in genes if g))[:5]
                report += f"| {condition} | {', '.join(unique_genes)} |\n"
        else:
            report += "Risk factors detected but not categorizable. See full disease risk report for details.\n"
    else:
        report += "No significant risk factors detected.\n"

    report += """

---

## Disclaimer

This protocol synthesizes genetic findings from multiple sources for informational purposes.
It is NOT a clinical diagnosis or medical advice.

- Genetic associations are probabilistic, not deterministic
- Environmental factors, lifestyle, and other genes also influence outcomes
- Classifications evolve as research progresses
- Consult healthcare providers before making medical decisions

---

*Generated by Genetic Health Analysis Pipeline - combining lifestyle genetics, PharmGKB, and ClinVar*
"""

    with open(output_path, 'w') as f:
        f.write(report)

    print(f"    Markdown: {output_path}")
    _write_html(report, output_path)
