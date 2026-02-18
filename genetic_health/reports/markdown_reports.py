"""Three markdown report generators + HTML sidecar output.

- generate_exhaustive_genetic_report()  -> EXHAUSTIVE_GENETIC_REPORT.{md,html}
- generate_disease_risk_report()        -> EXHAUSTIVE_DISEASE_RISK_REPORT.{md,html}
- generate_actionable_protocol()        -> ACTIONABLE_HEALTH_PROTOCOL_V3.{md,html}
"""

from pathlib import Path
from datetime import datetime
from collections import defaultdict

from ..ancestry import get_population_warnings

from .html_converter import _write_html
from .section_builders import (
    generate_executive_summary,
    generate_priority_findings,
    generate_pathway_analysis,
    generate_full_findings,
    generate_pharmgkb_report,
    generate_action_summary,
    generate_disclaimer,
)

# Display limits for report sections
MAX_RISK_FACTORS = 30
MAX_DRUG_RESPONSES = 30
MAX_LEVEL2_INTERACTIONS = 15
MAX_CLINVAR_DRUG_RESPONSES = 20


def _print_step(text):
    print(f"\n>>> {text}")


def _first_trait(traits, fallback='Unknown'):
    """Extract first trait from semicolon-separated list, stripped."""
    if not traits:
        return fallback
    return traits.split(';')[0].strip() or fallback


def _truncate(text, limit=80):
    """Truncate text to limit, adding ellipsis if needed."""
    text = text.strip()
    if len(text) > limit:
        return text[:limit] + '...'
    return text


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
            report += f"""### {f['gene']} - {_first_trait(f['traits'])}

- **Position:** chr{f['chromosome']}:{f['position']}
- **RSID:** {f['rsid']}
- **Genotype:** `{f['user_genotype']}` ({'Homozygous' if f['is_homozygous'] else 'Heterozygous'})
- **Variant:** {f['ref']} -> {f['alt']}
- **Confidence:** {stars} ({f['gold_stars']}/4)
- **Condition:** {(f['traits'].strip() if f['traits'] else 'Not specified')}

"""
            pw = get_population_warnings(f['gene'], 'pathogenic')
            if not pw and f['is_homozygous']:
                pw = get_population_warnings(f['gene'], 'homozygous')
            if pw:
                for w in pw:
                    report += f"> **Population Note:** {w}\n\n"

            report += "---\n\n"

    if carriers:
        report += "## Carrier Status - Recessive Conditions\n\n"
        report += "**You are a carrier - no personal symptoms expected, but reproductive implications.**\n\n"
        for f in carriers:
            stars = '*' * f['gold_stars'] + '.' * (4 - f['gold_stars'])
            report += f"""### {f['gene']} - {_first_trait(f['traits'])}

- **RSID:** {f['rsid']}
- **Genotype:** `{f['user_genotype']}` (Carrier)
- **Confidence:** {stars} ({f['gold_stars']}/4)
- **Condition:** {(f['traits'].strip() if f['traits'] else 'Not specified')}

"""
            pw = get_population_warnings(f['gene'], 'carrier')
            if pw:
                for w in pw:
                    report += f"> **Population Note:** {w}\n\n"

            report += "---\n\n"

    if het_unknown:
        report += "## Pathogenic/Likely Pathogenic - Inheritance Unclear\n\n"
        for f in het_unknown:
            stars = '*' * f['gold_stars'] + '.' * (4 - f['gold_stars'])
            report += f"""### {f['gene']} - {_first_trait(f['traits'])}

- **RSID:** {f['rsid']}
- **Genotype:** `{f['user_genotype']}`
- **Confidence:** {stars} ({f['gold_stars']}/4)
- **Condition:** {(f['traits'].strip() if f['traits'] else 'Not specified')}

---

"""

    if findings['risk_factor']:
        report += "## Risk Factor Variants\n\n"
        for f in findings['risk_factor'][:MAX_RISK_FACTORS]:
            report += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` - {_truncate(f['traits']) if f['traits'] else 'Risk factor'}...\n"
        report += "\n---\n\n"

    if findings['drug_response']:
        report += "## Drug Response Variants\n\n"
        for f in findings['drug_response'][:MAX_DRUG_RESPONSES]:
            report += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` - {_truncate(f['traits']) if f['traits'] else 'Drug response'}...\n"
        report += "\n---\n\n"

    if findings['protective']:
        report += "## Protective Variants\n\n"
        for f in findings['protective']:
            report += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` - {_truncate(f['traits']) if f['traits'] else 'Protective'}...\n"
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


def _render_recommendations_section(recommendations):
    """Render the personalized recommendations section as Markdown."""
    if not recommendations:
        return ""

    section = "\n\n---\n\n## Personalized Recommendations\n\n"
    section += ("*Prioritized synthesis of all genetic signals — convergent risks "
                "are ranked highest.*\n\n")

    priorities = recommendations.get("priorities", [])
    if not priorities:
        section += "No convergent risk patterns detected.\n"
    else:
        for p in priorities:
            icon = {"high": "!!!", "moderate": "!!", "low": ""}.get(p["priority"], "")
            section += f"### {p['title']} ({p['priority'].title()}) {icon}\n\n"
            section += f"**Why:** {p['why']}\n\n"
            section += "**Actions:**\n"
            for i, action in enumerate(p["actions"], 1):
                section += f"{i}. {action}\n"
            if p.get("doctor_note"):
                section += f"\n> **Doctor note:** {p['doctor_note']}\n"
            if p.get("monitoring"):
                section += "\n**Monitoring:**\n"
                for m in p["monitoring"]:
                    section += f"- {m['test']} ({m['frequency']}): {m['reason']}\n"
            section += "\n---\n\n"

    # Drug card
    drug_card = recommendations.get("drug_card", [])
    if drug_card:
        section += "### Drug-Gene Card\n\n"
        section += "| Gene | rsID | Genotype | Status | Source |\n"
        section += "|------|------|----------|--------|--------|\n"
        for entry in drug_card:
            for e in entry["entries"]:
                status = e.get("status", "").replace("_", " ").title() or "\u2014"
                section += (f"| {entry['gene']} | {e['rsid']} | "
                            f"`{e['genotype']}` | {status} | {e['source']} |\n")
        section += "\n"

    # Monitoring schedule
    schedule = recommendations.get("monitoring_schedule", [])
    if schedule:
        section += "### Monitoring Schedule\n\n"
        section += "| Test | Frequency | Reason |\n"
        section += "|------|-----------|--------|\n"
        for m in schedule:
            section += f"| {m['test']} | {m['frequency']} | {m['reason']} |\n"
        section += "\n"

    # Good news
    good_news = recommendations.get("good_news", [])
    if good_news:
        section += "### Good News\n\n"
        for g in good_news:
            section += f"- **{g['gene']}**: {g['description']}\n"
        section += "\n"

    return section


def generate_actionable_protocol(health_results: dict, disease_findings: dict,
                                  output_path: Path, subject_name: str = None,
                                  ancestry_results: dict = None, prs_results: dict = None,
                                  epistasis_results: list = None,
                                  recommendations: dict = None,
                                  quality_metrics: dict = None,
                                  blood_type: dict = None,
                                  mt_haplogroup: dict = None,
                                  star_alleles: dict = None):
    """Generate comprehensive actionable health protocol combining ALL sources."""
    _print_step("Generating actionable health protocol (comprehensive)")

    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    subject_line = f"\n**Subject:** {subject_name}" if subject_name else ""

    # Build gene->finding dict keeping highest-magnitude finding per gene
    # (some genes like MTHFR, CYP2C9 have multiple rsIDs)
    findings_dict = {}
    for f in health_results['findings']:
        gene = f['gene']
        if gene not in findings_dict or f['magnitude'] > findings_dict[gene]['magnitude']:
            findings_dict[gene] = f

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

"""

    # --- Data Quality Section ---
    if quality_metrics:
        report += "## Data Quality\n\n"
        report += f"- **Total SNPs loaded:** {quality_metrics['total_snps']:,}\n"
        report += f"- **Call rate:** {quality_metrics['call_rate']:.1%}\n"
        report += f"- **No-call positions:** {quality_metrics['no_call_count']:,}\n"
        report += f"- **Autosomal SNPs:** {quality_metrics['autosomal_count']:,}\n"
        report += f"- **Heterozygosity rate:** {quality_metrics['het_rate']:.3f}\n"
        report += f"- **Mitochondrial SNPs:** {quality_metrics['mt_snp_count']}\n"
        if quality_metrics['has_y']:
            report += "- **Inferred sex:** Male (Y chromosome detected)\n"
        else:
            report += "- **Inferred sex:** Female (no Y chromosome)\n"

        # Chromosome coverage summary
        chroms = quality_metrics.get('chromosomes', {})
        autosomal = sorted(
            [(ch, cnt) for ch, cnt in chroms.items() if ch.isdigit()],
            key=lambda x: int(x[0])
        )
        if autosomal:
            report += "\n**Chromosome coverage (top 5):**\n\n"
            report += "| Chromosome | SNPs |\n|------------|------|\n"
            for ch, cnt in sorted(autosomal, key=lambda x: -x[1])[:5]:
                report += f"| chr{ch} | {cnt:,} |\n"

        report += "\n---\n\n"

    report += """## Executive Summary

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
            condition = _first_trait(f['traits'], 'Unknown condition')
            stars = f['gold_stars']
            confidence = f"({stars}/4 stars)" if stars > 0 else "(low confidence)"
            report += f"- **{f['gene']}**: {condition} {confidence}\n"
        report += "\n"

    if carriers:
        report += "**Carrier Status (Recessive):**\n"
        for f in carriers:
            condition = _first_trait(f['traits'], 'Unknown condition')
            stars = f['gold_stars']
            confidence = f"({stars}/4 stars)" if stars > 0 else "(low confidence)"
            report += f"- **{f['gene']}**: {condition} {confidence}\n"
        report += "\n"

    if het_unknown:
        report += "**Heterozygous (Inheritance Unclear):**\n"
        for f in het_unknown:
            condition = _first_trait(f['traits'], 'Unknown condition')
            stars = f['gold_stars']
            confidence = f"({stars}/4 stars)" if stars > 0 else "(low confidence)"
            report += f"- **{f['gene']}**: {condition} {confidence}\n"
        report += "\n"

    if not affected and not carriers and not het_unknown:
        report += "None detected.\n\n"

    report += "### Protective Variants\n\n"
    if disease_findings and disease_findings.get('protective'):
        for f in disease_findings['protective']:
            condition = _first_trait(f['traits'], 'Protective effect')
            report += f"- **{f['gene']}**: {condition}\n"
    else:
        report += "None detected.\n"

    # --- Personalized Recommendations ---
    report += _render_recommendations_section(recommendations)

    # --- Ancestry Section ---
    report += "\n\n---\n\n## Ancestry Estimation\n\n"

    if ancestry_results and ancestry_results.get('markers_found', 0) > 0:
        report += f"**Top Ancestry:** {ancestry_results['top_ancestry']}\n"
        report += f"**Confidence:** {ancestry_results['confidence'].title()} "
        report += f"({ancestry_results['markers_found']} markers analyzed)\n\n"

        report += "| Population | Proportion |\n"
        report += "|------------|------------|\n"
        for pop in sorted(ancestry_results['proportions'], key=lambda p: -ancestry_results['proportions'][p]):
            prop = ancestry_results['proportions'][pop]
            from ..ancestry import POPULATION_LABELS
            label = POPULATION_LABELS.get(pop, pop)
            report += f"| {label} ({pop}) | {prop:.1%} |\n"

        report += "\n*Ancestry estimation uses ~55 ancestry-informative markers (AIMs). "
        report += "This is a rough superpopulation estimate, not a detailed ethnicity breakdown. "
        report += "Commercial ancestry services use hundreds of thousands of markers for finer resolution.*\n"
    else:
        report += "Insufficient markers for ancestry estimation.\n"

    # --- Blood Type Section ---
    report += "\n\n---\n\n## Blood Type\n\n"

    if blood_type and blood_type.get("blood_type") != "Unknown":
        report += f"**Predicted Blood Type:** {blood_type['blood_type']}\n\n"
        report += f"- **ABO Group:** {blood_type['abo']}\n"
        report += f"- **Rh Factor:** {blood_type['rh']}\n"
        report += f"- **Confidence:** {blood_type['confidence'].title()}\n\n"
        if blood_type.get("details"):
            report += "**SNPs used:**\n"
            for d in blood_type["details"]:
                report += f"- {d}\n"
        report += "\n*Blood type prediction from SNP data has limitations. "
        report += "23andMe cannot detect all ABO variants (rs8176719 is an indel). "
        report += "Confirm with a clinical blood typing test.*\n"
    else:
        report += "Insufficient SNP data for blood type prediction.\n"

    # --- Mitochondrial Haplogroup Section ---
    report += "\n\n---\n\n## Mitochondrial Haplogroup (Maternal Lineage)\n\n"

    if mt_haplogroup and mt_haplogroup.get("haplogroup") != "Unknown":
        report += f"**Haplogroup:** {mt_haplogroup['haplogroup']}\n\n"
        report += f"- **Description:** {mt_haplogroup['description']}\n"
        report += f"- **Lineage:** {mt_haplogroup['lineage']}\n"
        report += f"- **Confidence:** {mt_haplogroup['confidence'].title()}\n"
        report += f"- **Markers found:** {mt_haplogroup['markers_found']}/{mt_haplogroup['markers_tested']}\n\n"
        report += "*Mitochondrial DNA is inherited exclusively from the mother and traces "
        report += "the direct maternal line. Your haplogroup reflects ancient migration patterns.*\n"
    else:
        report += "Insufficient mitochondrial SNP data for haplogroup estimation.\n"

    # --- PRS Section ---
    report += "\n\n---\n\n## Polygenic Risk Scores\n\n"

    if prs_results:
        report += "*PRS models estimate relative genetic risk based on common variants. "
        report += "They do not account for lifestyle, environment, or rare variants.*\n\n"

        # Check for ancestry warning
        any_non_applicable = any(not r['ancestry_applicable'] for r in prs_results.values())
        if any_non_applicable:
            report += "> **Note:** Your ancestry profile is substantially non-European. "
            report += "These PRS models were calibrated on European-ancestry populations "
            report += "and may not accurately reflect your risk.\n\n"

        report += "| Condition | Percentile | Risk Category | SNPs Found | Reference |\n"
        report += "|-----------|-----------|---------------|------------|----------|\n"
        for cid, r in prs_results.items():
            flag = "" if r['ancestry_applicable'] else " *"
            report += (f"| {r['name']} | {r['percentile']:.0f}th | "
                       f"{r['risk_category'].title()}{flag} | "
                       f"{r['snps_found']}/{r['snps_total']} | "
                       f"{r['reference']} |\n")

        # Highlight elevated/high
        elevated = [r for r in prs_results.values() if r['risk_category'] in ('elevated', 'high')]
        if elevated:
            report += "\n### Elevated Risk Conditions\n\n"
            for r in elevated:
                report += f"**{r['name']}** — {r['percentile']:.0f}th percentile ({r['risk_category']})\n\n"
                if r['contributing_snps']:
                    report += "Top contributing variants:\n"
                    for s in r['contributing_snps'][:5]:
                        report += f"- **{s['gene']}** ({s['rsid']}): {s['copies']} copies of risk allele\n"
                report += "\n"
    else:
        report += "PRS calculation not available.\n"

    # --- Epistasis Section ---
    report += "\n\n---\n\n## Gene-Gene Interactions (Epistasis)\n\n"

    if epistasis_results:
        report += ("*These interactions occur when the combined effect of multiple "
                   "gene variants differs from each individual effect.*\n\n")

        for interaction in epistasis_results:
            risk_icon = {"high": "!!!", "moderate": "!!", "low": "!"}.get(
                interaction['risk_level'], "")
            report += f"### {interaction['name']} {risk_icon}\n\n"
            report += f"**Risk Level:** {interaction['risk_level'].title()}\n\n"
            report += f"**Genes:** {', '.join(interaction['genes_involved'].keys())}\n\n"
            report += f"**Effect:** {interaction['effect']}\n\n"
            report += f"**Mechanism:** {interaction['mechanism']}\n\n"
            report += "**Recommended Actions:**\n"
            for action in interaction['actions']:
                report += f"- {action}\n"
            report += "\n---\n\n"
    else:
        report += "No significant gene-gene interactions detected.\n"

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

    # --- Star Alleles Section ---
    report += "\n\n---\n\n## Pharmacogenomic Star Alleles\n\n"

    if star_alleles:
        report += ("*CPIC-style star allele calling for key drug-metabolizing enzymes. "
                   "Share with prescribing physicians.*\n\n")

        report += "| Gene | Diplotype | Phenotype | SNPs Found | Note |\n"
        report += "|------|-----------|-----------|------------|------|\n"
        for gene, r in star_alleles.items():
            phenotype = r['phenotype'].replace('_', ' ').title()
            note = r['clinical_note'][:80] + '...' if len(r['clinical_note']) > 80 else r['clinical_note']
            report += (f"| {gene} | {r['diplotype']} | {phenotype} | "
                       f"{r['snps_found']}/{r['snps_total']} | {note} |\n")

        # Highlight actionable results
        actionable = [
            (g, r) for g, r in star_alleles.items()
            if r['phenotype'] in ('poor', 'intermediate', 'rapid', 'ultrarapid')
               and r['phenotype'] != 'Unknown'
        ]
        if actionable:
            report += "\n### Actionable Metabolizer Phenotypes\n\n"
            for gene, r in actionable:
                phenotype = r['phenotype'].replace('_', ' ').title()
                report += f"**{gene}** — {r['diplotype']} ({phenotype} Metabolizer)\n\n"
                report += f"{r['clinical_note']}\n\n"
    else:
        report += "Star allele data not available.\n"

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
            drugs = _truncate(f['drugs'], 50)
            report += f"| {f['gene']} | {f['level']} | {drugs} | `{f['genotype']}` |\n"
    else:
        report += "None detected.\n"

    report += "\n### PharmGKB Level 2 (Moderate Evidence)\n\n"

    level_2 = [f for f in health_results['pharmgkb_findings'] if f['level'] in ['2A', '2B']]
    if level_2:
        report += "| Gene | Level | Drugs | Your Genotype |\n"
        report += "|------|-------|-------|---------------|\n"
        for f in level_2[:MAX_LEVEL2_INTERACTIONS]:
            drugs = _truncate(f['drugs'], 50)
            report += f"| {f['gene']} | {f['level']} | {drugs} | `{f['genotype']}` |\n"
        if len(level_2) > MAX_LEVEL2_INTERACTIONS:
            report += f"\n*...and {len(level_2) - MAX_LEVEL2_INTERACTIONS} more Level 2 interactions*\n"
    else:
        report += "None detected.\n"

    report += "\n### ClinVar Drug Response Variants\n\n"

    if disease_findings and disease_findings.get('drug_response'):
        drug_resp = disease_findings['drug_response']
        report += "| Gene | RSID | Genotype | Drug/Response |\n"
        report += "|------|------|----------|---------------|\n"
        for f in drug_resp[:MAX_CLINVAR_DRUG_RESPONSES]:
            traits = _truncate(f['traits'], 60)
            gene = f['gene'] if f['gene'] else '\u2014'
            report += f"| {gene} | {f['rsid']} | `{f['user_genotype']}` | {traits} |\n"
        if len(drug_resp) > MAX_CLINVAR_DRUG_RESPONSES:
            report += f"\n*...and {len(drug_resp) - MAX_CLINVAR_DRUG_RESPONSES} more drug response variants*\n"
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
