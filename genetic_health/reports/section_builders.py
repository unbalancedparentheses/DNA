"""Section generators for the exhaustive genetic report.

Extracted from the original generate_exhaustive_report.py — data-building
functions only (CLINICAL_CONTEXT/PATHWAYS live in clinical_context.py).
"""

from datetime import datetime
from collections import defaultdict

from ..clinical_context import get_clinical_context, get_related_pathways, PATHWAYS
from ..ancestry import get_population_warnings


def format_magnitude(mag):
    """Format magnitude with color indicator."""
    if mag >= 3:
        return f"\U0001f534 HIGH ({mag}/6)"
    elif mag == 2:
        return f"\U0001f7e1 MODERATE ({mag}/6)"
    elif mag == 1:
        return f"\U0001f7e2 LOW ({mag}/6)"
    else:
        return f"\u26aa NEUTRAL ({mag}/6)"


def format_evidence_level(level):
    """Format PharmGKB evidence level."""
    if level == "1A":
        return "\U0001f535 1A - Clinical guideline annotation"
    elif level == "1B":
        return "\U0001f535 1B - Clinical guideline annotation"
    elif level == "2A":
        return "\U0001f7e3 2A - Variant has moderate evidence"
    elif level == "2B":
        return "\U0001f7e3 2B - Variant has moderate evidence"
    else:
        return f"\u26ab {level}"


def generate_finding_section(finding, index):
    """Generate a comprehensive section for a single finding."""
    gene = finding.get('gene', 'Unknown')
    rsid = finding.get('rsid', '')
    category = finding.get('category', 'Uncategorized')
    genotype = finding.get('genotype', '')
    status = finding.get('status', '')
    description = finding.get('description', '')
    magnitude = finding.get('magnitude', 0)
    note = finding.get('note', '')

    section = []
    section.append(f"### {index}. {gene} ({rsid})")
    section.append("")
    section.append(f"**Category:** {category}  ")
    section.append(f"**Your Genotype:** `{genotype}`  ")
    section.append(f"**Status:** {status.replace('_', ' ').title()}  ")
    section.append(f"**Impact:** {format_magnitude(magnitude)}")
    section.append("")
    section.append(f"**Description:** {description}")

    if note:
        section.append(f"")
        section.append(f"**Note:** {note}")

    pathways = get_related_pathways(gene)
    if pathways:
        section.append("")
        section.append(f"**Related Pathways:** {', '.join(pathways)}")

    context = get_clinical_context(gene, status)
    if context:
        section.append("")
        section.append("#### Mechanism")
        section.append(context['mechanism'])

        if context.get('implications'):
            section.append("")
            section.append("#### Implications")
            for imp in context['implications']:
                section.append(f"- {imp}")

        if context.get('actions'):
            section.append("")
            section.append("#### Recommended Actions")
            for action in context['actions']:
                section.append(f"- {action}")

        if context.get('interactions'):
            section.append("")
            section.append("#### Gene Interactions")
            for interaction in context['interactions']:
                section.append(f"- {interaction}")

    pop_warnings = get_population_warnings(gene, status)
    if pop_warnings:
        section.append("")
        section.append("#### Population Note")
        for w in pop_warnings:
            section.append(f"- {w}")

    section.append("")
    section.append("---")
    section.append("")

    return "\n".join(section)


def generate_pharmgkb_section(finding, index):
    """Generate a comprehensive section for a PharmGKB drug interaction."""
    gene = finding.get('gene', 'Unknown')
    rsid = finding.get('rsid', '')
    drugs = finding.get('drugs', '')
    genotype = finding.get('genotype', '')
    annotation = finding.get('annotation', '')
    level = finding.get('level', '')
    category = finding.get('category', 'Other')

    section = []
    section.append(f"### {index}. {gene} - {rsid}")
    section.append("")
    section.append(f"**Evidence Level:** {format_evidence_level(level)}  ")
    section.append(f"**Category:** {category}  ")
    section.append(f"**Your Genotype:** `{genotype}`  ")
    section.append(f"**Affected Drugs:** {drugs}")
    section.append("")
    section.append("#### Clinical Annotation")
    section.append(annotation)
    section.append("")

    if level in ["1A", "1B"]:
        section.append("#### Clinical Significance")
        section.append("This is a high-evidence drug-gene interaction with clinical guideline support. Discuss with prescribing physicians before starting these medications.")

    section.append("---")
    section.append("")

    return "\n".join(section)


def generate_executive_summary(data):
    """Generate executive summary."""
    findings = data.get('findings', [])
    pharmgkb = data.get('pharmgkb_findings', [])
    summary = data.get('summary', {})

    high_impact = [f for f in findings if f.get('magnitude', 0) >= 3]
    mod_impact = [f for f in findings if f.get('magnitude', 0) == 2]
    low_impact = [f for f in findings if f.get('magnitude', 0) == 1]

    level_1 = [f for f in pharmgkb if f.get('level', '').startswith('1')]
    level_2 = [f for f in pharmgkb if f.get('level', '').startswith('2')]

    categories = set(f.get('category') for f in findings)

    lines = []
    lines.append("# Exhaustive Genetic Health Report")
    lines.append("")
    lines.append(f"**Generated:** {datetime.now().strftime('%Y-%m-%d %H:%M')}")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## Executive Summary")
    lines.append("")
    lines.append("### Genome Overview")
    lines.append(f"- **Total SNPs in Raw Data:** {summary.get('total_snps', 'N/A'):,}")
    lines.append(f"- **Clinically Relevant SNPs Analyzed:** {len(findings)}")
    lines.append(f"- **PharmGKB Drug Interactions:** {len(pharmgkb)}")
    lines.append("")
    lines.append("### Impact Distribution")
    lines.append(f"- \U0001f534 **High Impact (magnitude \u22653):** {len(high_impact)}")
    lines.append(f"- \U0001f7e1 **Moderate Impact (magnitude 2):** {len(mod_impact)}")
    lines.append(f"- \U0001f7e2 **Low Impact (magnitude 1):** {len(low_impact)}")
    lines.append(f"- \u26aa **Informational (magnitude 0):** {len(findings) - len(high_impact) - len(mod_impact) - len(low_impact)}")
    lines.append("")
    lines.append("### Pharmacogenomics")
    lines.append(f"- \U0001f535 **Level 1 (Clinical Guidelines):** {len(level_1)}")
    lines.append(f"- \U0001f7e3 **Level 2 (Moderate Evidence):** {len(level_2)}")
    lines.append("")
    lines.append("### Categories Covered")
    for cat in sorted(categories):
        count = len([f for f in findings if f.get('category') == cat])
        lines.append(f"- {cat}: {count} findings")
    lines.append("")
    lines.append("---")
    lines.append("")

    return "\n".join(lines)


def generate_priority_findings(findings):
    """Generate the priority findings section."""
    high_impact = [f for f in findings if f.get('magnitude', 0) >= 3]
    mod_impact = [f for f in findings if f.get('magnitude', 0) == 2]

    lines = []
    lines.append("## \U0001f534 Priority Findings (High Impact)")
    lines.append("")
    lines.append("These findings have the most significant implications for your health decisions.")
    lines.append("")

    for i, finding in enumerate(high_impact, 1):
        lines.append(generate_finding_section(finding, i))

    if mod_impact:
        lines.append("## \U0001f7e1 Moderate Impact Findings")
        lines.append("")
        lines.append("These findings warrant attention and may influence health decisions.")
        lines.append("")

        for i, finding in enumerate(mod_impact, 1):
            lines.append(generate_finding_section(finding, i))

    return "\n".join(lines)


def generate_full_findings(findings):
    """Generate all findings organized by category."""
    categories = sorted(set(f.get('category', 'Other') for f in findings))

    lines = []
    lines.append("## Complete Findings by Category")
    lines.append("")
    lines.append("Every genetic finding analyzed, organized by category.")
    lines.append("")

    for category in categories:
        cat_findings = [f for f in findings if f.get('category') == category]
        if not cat_findings:
            continue

        lines.append(f"### \U0001f4c2 {category}")
        lines.append("")

        cat_findings.sort(key=lambda x: x.get('magnitude', 0), reverse=True)

        for i, finding in enumerate(cat_findings, 1):
            gene = finding.get('gene', 'Unknown')
            rsid = finding.get('rsid', '')
            genotype = finding.get('genotype', '')
            status = finding.get('status', '').replace('_', ' ').title()
            description = finding.get('description', '')
            magnitude = finding.get('magnitude', 0)

            mag_icon = "\U0001f534" if magnitude >= 3 else "\U0001f7e1" if magnitude == 2 else "\U0001f7e2" if magnitude == 1 else "\u26aa"

            lines.append(f"#### {i}. {gene} ({rsid}) {mag_icon}")
            lines.append(f"- **Genotype:** `{genotype}` | **Status:** {status} | **Impact:** {magnitude}/6")
            lines.append(f"- {description}")

            context = get_clinical_context(gene, finding.get('status', ''))
            if context:
                lines.append(f"- **Mechanism:** {context['mechanism'][:200]}...")
                if context.get('actions'):
                    lines.append(f"- **Key Action:** {context['actions'][0]}")

            lines.append("")

        lines.append("---")
        lines.append("")

    return "\n".join(lines)


def generate_pharmgkb_report(pharmgkb):
    """Generate comprehensive pharmacogenomics section."""
    lines = []
    lines.append("## \U0001f48a Pharmacogenomics - Complete Drug-Gene Interactions")
    lines.append("")
    lines.append("This section contains all drug-gene interactions from PharmGKB with clinical annotations.")
    lines.append("Share this information with prescribing physicians before starting new medications.")
    lines.append("")

    level_1a = [f for f in pharmgkb if f.get('level') == '1A']
    level_1b = [f for f in pharmgkb if f.get('level') == '1B']
    level_2a = [f for f in pharmgkb if f.get('level') == '2A']
    level_2b = [f for f in pharmgkb if f.get('level') == '2B']

    if level_1a:
        lines.append("### Level 1A - Highest Evidence (Clinical Guideline Annotations)")
        lines.append("")
        for i, finding in enumerate(level_1a, 1):
            lines.append(generate_pharmgkb_section(finding, i))

    if level_1b:
        lines.append("### Level 1B - High Evidence (Clinical Guideline Annotations)")
        lines.append("")
        for i, finding in enumerate(level_1b, 1):
            lines.append(generate_pharmgkb_section(finding, i))

    if level_2a:
        lines.append("### Level 2A - Moderate Evidence")
        lines.append("")
        for i, finding in enumerate(level_2a, 1):
            lines.append(generate_pharmgkb_section(finding, i))

    if level_2b:
        lines.append("### Level 2B - Moderate Evidence")
        lines.append("")
        for i, finding in enumerate(level_2b, 1):
            lines.append(generate_pharmgkb_section(finding, i))

    return "\n".join(lines)


def generate_pathway_analysis(findings):
    """Generate pathway analysis section."""
    lines = []
    lines.append("## \U0001f517 Pathway Analysis")
    lines.append("")
    lines.append("Your genes grouped by biological pathway, showing how multiple variants may interact.")
    lines.append("")

    for pathway_name, pathway_genes in PATHWAYS.items():
        pathway_findings = [f for f in findings if f.get('gene') in pathway_genes]
        if not pathway_findings:
            continue

        lines.append(f"### {pathway_name}")
        lines.append("")

        for finding in pathway_findings:
            gene = finding.get('gene', '')
            status = finding.get('status', '').replace('_', ' ').title()
            magnitude = finding.get('magnitude', 0)
            mag_icon = "\U0001f534" if magnitude >= 3 else "\U0001f7e1" if magnitude == 2 else "\U0001f7e2" if magnitude == 1 else "\u26aa"

            lines.append(f"- {mag_icon} **{gene}:** {status}")

        lines.append("")

        if pathway_name == "Methylation Cycle":
            mthfr = next((f for f in pathway_findings if f.get('gene') == 'MTHFR' and f.get('magnitude', 0) >= 2), None)
            mtrr = next((f for f in pathway_findings if f.get('gene') == 'MTRR' and f.get('magnitude', 0) >= 2), None)
            if mthfr and mtrr:
                lines.append("\u26a0\ufe0f **Pathway Impact:** Multiple methylation cycle variants detected. Consider comprehensive methylation support (methylfolate + methylcobalamin + B2).")
                lines.append("")

        elif pathway_name == "Blood Pressure":
            bp_findings = [f for f in pathway_findings if f.get('magnitude', 0) >= 1]
            if len(bp_findings) >= 2:
                lines.append("\u26a0\ufe0f **Pathway Impact:** Multiple blood pressure-related variants. Recommend regular monitoring and lifestyle optimization.")
                lines.append("")

        lines.append("---")
        lines.append("")

    return "\n".join(lines)


def generate_action_summary(findings):
    """Generate comprehensive action summary."""
    lines = []
    lines.append("## \U0001f4cb Comprehensive Action Summary")
    lines.append("")

    supplement_actions = []
    diet_actions = []
    lifestyle_actions = []
    monitoring_actions = []
    medical_actions = []

    for finding in findings:
        context = get_clinical_context(finding.get('gene'), finding.get('status'))
        if context and context.get('actions'):
            for action in context['actions']:
                action_lower = action.lower()
                if any(word in action_lower for word in ['supplement', 'vitamin', 'mg', 'mcg', 'iu', 'dose']):
                    supplement_actions.append(f"- {action} *(from {finding.get('gene')})*")
                elif any(word in action_lower for word in ['diet', 'eat', 'food', 'limit', 'avoid', 'meal']):
                    diet_actions.append(f"- {action} *(from {finding.get('gene')})*")
                elif any(word in action_lower for word in ['exercise', 'sleep', 'stress', 'meditation']):
                    lifestyle_actions.append(f"- {action} *(from {finding.get('gene')})*")
                elif any(word in action_lower for word in ['test', 'monitor', 'check', 'measure']):
                    monitoring_actions.append(f"- {action} *(from {finding.get('gene')})*")
                elif any(word in action_lower for word in ['doctor', 'physician', 'medical', 'prescrib']):
                    medical_actions.append(f"- {action} *(from {finding.get('gene')})*")

    if supplement_actions:
        lines.append("### \U0001f48a Supplement Considerations")
        lines.append("*Discuss with healthcare provider before starting*")
        lines.append("")
        lines.extend(list(set(supplement_actions))[:15])
        lines.append("")

    if diet_actions:
        lines.append("### \U0001f957 Dietary Recommendations")
        lines.append("")
        lines.extend(list(set(diet_actions))[:10])
        lines.append("")

    if lifestyle_actions:
        lines.append("### \U0001f3c3 Lifestyle Actions")
        lines.append("")
        lines.extend(list(set(lifestyle_actions))[:10])
        lines.append("")

    if monitoring_actions:
        lines.append("### \U0001f4ca Monitoring Recommendations")
        lines.append("")
        lines.extend(list(set(monitoring_actions))[:10])
        lines.append("")

    if medical_actions:
        lines.append("### \U0001f3e5 Medical Considerations")
        lines.append("")
        lines.extend(list(set(medical_actions))[:10])
        lines.append("")

    lines.append("---")
    lines.append("")

    return "\n".join(lines)


def generate_disclaimer():
    """Generate disclaimer section."""
    return """\
## \u26a0\ufe0f Important Disclaimer

This report is for **informational and educational purposes only**. It is NOT medical advice.

### Key Points:
- Genetic associations are probabilistic, not deterministic
- Your genes are just one factor - environment, lifestyle, and other genes matter
- "Risk" variants don't guarantee outcomes; "protective" variants don't guarantee safety
- Consult healthcare providers before making medical decisions
- Some findings may have different implications in different populations
- Genetic science evolves - recommendations may change as research advances

### How to Use This Report:
1. **Share with providers** - Especially the pharmacogenomics section before new medications
2. **Focus on actionable items** - Prioritize evidence-based interventions
3. **Don't over-interpret** - One gene doesn't define your health destiny
4. **Combine with testing** - Many recommendations include follow-up lab tests

---

*Report generated using comprehensive genetic analysis. Data source: Personal genome analysis with PharmGKB clinical annotations.*
"""
