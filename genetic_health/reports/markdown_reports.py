"""Unified markdown report generator + HTML sidecar output.

Produces a single comprehensive GENETIC_HEALTH_REPORT.md combining all
analysis modules into one report.
"""

from pathlib import Path
from datetime import datetime
from collections import defaultdict

from ..ancestry import get_population_warnings, POPULATION_LABELS
from ..clinical_context import get_clinical_context, get_related_pathways, PATHWAYS

from .html_converter import _write_html

# Display limits
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


def generate_unified_report(
    health_results, disease_findings, disease_stats,
    output_path, subject_name=None,
    ancestry_results=None, prs_results=None,
    epistasis_results=None, recommendations=None,
    quality_metrics=None, blood_type=None,
    mt_haplogroup=None, star_alleles=None,
    apoe=None, acmg=None, carrier_screen=None, traits=None,
):
    """Generate the unified comprehensive genetic health report."""
    _print_step("Generating unified genetic health report")

    now = datetime.now().strftime("%Y-%m-%d %H:%M")
    subject_line = f"\n**Subject:** {subject_name}" if subject_name else ""

    # Classify disease findings
    affected, carriers, het_unknown = [], [], []
    if disease_findings:
        for f in disease_findings.get('pathogenic', []) + disease_findings.get('likely_pathogenic', []):
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

    # Build gene->finding dict for supplement/diet logic
    findings_dict = {}
    for f in health_results['findings']:
        gene = f['gene']
        if gene not in findings_dict or f['magnitude'] > findings_dict[gene]['magnitude']:
            findings_dict[gene] = f

    report = f"""# Genetic Health Report
{subject_line}
**Generated:** {now}

---

"""
    # === 1. Data Quality ===
    report += _section_quality(quality_metrics)

    # === 2. Executive Summary ===
    report += _section_executive_summary(
        health_results, affected, carriers, het_unknown, disease_findings, apoe, acmg)

    # === 3. APOE Haplotype ===
    report += _section_apoe(apoe)

    # === 4. Blood Type + Traits ===
    report += _section_blood_type_and_traits(blood_type, traits)

    # === 5. Mitochondrial Haplogroup ===
    report += _section_mt_haplogroup(mt_haplogroup)

    # === 6. Ancestry Estimation ===
    report += _section_ancestry(ancestry_results)

    # === 7. Polygenic Risk Scores ===
    report += _section_prs(prs_results)

    # === 8. Pharmacogenomic Star Alleles ===
    report += _section_star_alleles(star_alleles)

    # === 9. ACMG Secondary Findings ===
    report += _section_acmg(acmg)

    # === 10. Epistasis ===
    report += _section_epistasis(epistasis_results)

    # === 11. Personalized Recommendations ===
    report += _section_recommendations(recommendations, findings_dict)

    # === 12. Complete Lifestyle/Health Findings ===
    report += _section_lifestyle_findings(health_results['findings'])

    # === 13. Pathway Analysis ===
    report += _section_pathway_analysis(health_results['findings'])

    # === 14. Disease Risk Analysis ===
    report += _section_disease_risk(
        affected, carriers, het_unknown, disease_findings, disease_stats,
        health_results['summary'].get('total_snps', 0) if isinstance(health_results['summary'], dict) else 0)

    # === 15. Carrier Screening ===
    report += _section_carrier_screening(carrier_screen)

    # === 16. Drug-Gene Interactions ===
    report += _section_drug_interactions(health_results, disease_findings)

    # === 17. Doctor Card ===
    report += _section_doctor_card(
        quality_metrics, blood_type, apoe, star_alleles, acmg,
        affected, carriers, prs_results)

    # === 18. References & Disclaimer ===
    report += _section_disclaimer()

    with open(output_path, 'w') as f:
        f.write(report)

    print(f"    Markdown: {output_path}")
    _write_html(report, output_path)


# =============================================================================
# Individual section generators
# =============================================================================

def _section_quality(quality_metrics):
    if not quality_metrics:
        return ""
    s = "## 1. Data Quality\n\n"
    s += f"- **Total SNPs loaded:** {quality_metrics['total_snps']:,}\n"
    s += f"- **Call rate:** {quality_metrics['call_rate']:.1%}\n"
    s += f"- **No-call positions:** {quality_metrics['no_call_count']:,}\n"
    s += f"- **Autosomal SNPs:** {quality_metrics['autosomal_count']:,}\n"
    s += f"- **Heterozygosity rate:** {quality_metrics['het_rate']:.3f}\n"
    s += f"- **Mitochondrial SNPs:** {quality_metrics['mt_snp_count']}\n"
    if quality_metrics['has_y']:
        s += "- **Inferred sex:** Male (Y chromosome detected)\n"
    else:
        s += "- **Inferred sex:** Female (no Y chromosome)\n"

    chroms = quality_metrics.get('chromosomes', {})
    autosomal = sorted(
        [(ch, cnt) for ch, cnt in chroms.items() if ch.isdigit()],
        key=lambda x: int(x[0])
    )
    if autosomal:
        s += "\n**Chromosome coverage (top 5):**\n\n"
        s += "| Chromosome | SNPs |\n|------------|------|\n"
        for ch, cnt in sorted(autosomal, key=lambda x: -x[1])[:5]:
            s += f"| chr{ch} | {cnt:,} |\n"

    s += "\n---\n\n"
    return s


def _section_executive_summary(health_results, affected, carriers, het_unknown,
                                disease_findings, apoe, acmg):
    findings = health_results['findings']
    pharmgkb = health_results['pharmgkb_findings']
    high_impact = [f for f in findings if f.get('magnitude', 0) >= 3]
    total_lifestyle = len(findings)
    total_pharmgkb = len(pharmgkb)

    s = "## 2. Executive Summary\n\n"

    s += f"- **Lifestyle/health findings:** {total_lifestyle}\n"
    s += f"- **Drug-gene interactions (PharmGKB):** {total_pharmgkb}\n"
    s += f"- **High-impact findings (magnitude >= 3):** {len(high_impact)}\n"
    s += f"- **Pathogenic variants (affected):** {len(affected)}\n"
    s += f"- **Carrier status:** {len(carriers)}\n"

    if apoe and apoe.get('apoe_type') != 'Unknown':
        s += f"- **APOE haplotype:** {apoe['apoe_type']} ({apoe['risk_level']} Alzheimer's risk)\n"

    if acmg and acmg.get('genes_with_variants', 0) > 0:
        s += f"- **ACMG actionable genes:** {acmg['genes_with_variants']} variant(s) found\n"

    if high_impact:
        s += "\n### High-Impact Lifestyle Findings\n\n"
        for f in high_impact:
            s += f"- **{f['gene']}** ({f['category']}): {f['description']}\n"

    if affected:
        s += "\n### Critical Disease Variants\n\n"
        for f in affected:
            condition = _first_trait(f['traits'], 'Unknown condition')
            stars = f['gold_stars']
            s += f"- **{f['gene']}**: {condition} ({stars}/4 stars)\n"

    s += "\n---\n\n"
    return s


def _section_apoe(apoe):
    s = "## 3. APOE Haplotype\n\n"
    if not apoe or apoe.get('apoe_type') == 'Unknown':
        s += "Insufficient SNP data for APOE haplotype determination.\n"
    else:
        s += f"**Haplotype:** {apoe['apoe_type']}\n\n"
        s += f"- **Risk level:** {apoe['risk_level'].title()}\n"
        if apoe.get('alzheimer_or') is not None:
            s += f"- **Alzheimer's odds ratio:** {apoe['alzheimer_or']}x\n"
        s += f"- **Confidence:** {apoe['confidence'].title()}\n\n"
        s += f"{apoe['description']}\n\n"
        s += ("*APOE is the strongest common genetic risk factor for late-onset Alzheimer's disease. "
              "Risk is modifiable through cardiovascular health, exercise, sleep, and cognitive engagement.*\n")
    s += "\n---\n\n"
    return s


def _section_blood_type_and_traits(blood_type, traits):
    s = "## 4. Blood Type & Traits\n\n"

    # Blood type
    s += "### Blood Type\n\n"
    if blood_type and blood_type.get("blood_type") != "Unknown":
        s += f"**Predicted Blood Type:** {blood_type['blood_type']}\n\n"
        s += f"- **ABO Group:** {blood_type['abo']}\n"
        s += f"- **Rh Factor:** {blood_type['rh']}\n"
        s += f"- **Confidence:** {blood_type['confidence'].title()}\n\n"
        if blood_type.get("details"):
            for d in blood_type["details"]:
                s += f"- {d}\n"
        s += "\n*Confirm with clinical blood typing.*\n\n"
    else:
        s += "Insufficient SNP data for blood type prediction.\n\n"

    # Traits
    if traits:
        s += "### Predicted Traits\n\n"
        s += "| Trait | Prediction | Confidence | Description |\n"
        s += "|-------|-----------|------------|-------------|\n"
        trait_labels = {
            "eye_color": "Eye Color",
            "hair_color": "Hair Color",
            "earwax_type": "Earwax Type",
            "freckling": "Freckling/Sun Sensitivity",
        }
        for trait_id, trait in traits.items():
            label = trait_labels.get(trait_id, trait_id.replace('_', ' ').title())
            desc = _truncate(trait['description'], 60)
            s += f"| {label} | {trait['prediction']} | {trait['confidence'].title()} | {desc} |\n"
        s += "\n"

    s += "\n---\n\n"
    return s


def _section_mt_haplogroup(mt_haplogroup):
    s = "## 5. Mitochondrial Haplogroup (Maternal Lineage)\n\n"
    if mt_haplogroup and mt_haplogroup.get("haplogroup") != "Unknown":
        s += f"**Haplogroup:** {mt_haplogroup['haplogroup']}\n\n"
        s += f"- **Description:** {mt_haplogroup['description']}\n"
        s += f"- **Lineage:** {mt_haplogroup['lineage']}\n"
        s += f"- **Confidence:** {mt_haplogroup['confidence'].title()}\n"
        s += f"- **Markers found:** {mt_haplogroup['markers_found']}/{mt_haplogroup['markers_tested']}\n\n"
        s += "*Mitochondrial DNA traces the direct maternal line.*\n"
    else:
        s += "Insufficient mitochondrial SNP data for haplogroup estimation.\n"
    s += "\n---\n\n"
    return s


def _section_ancestry(ancestry_results):
    s = "## 6. Ancestry Estimation\n\n"
    if ancestry_results and ancestry_results.get('markers_found', 0) > 0:
        s += f"**Top Ancestry:** {ancestry_results['top_ancestry']}\n"
        s += f"**Confidence:** {ancestry_results['confidence'].title()} "
        s += f"({ancestry_results['markers_found']} markers analyzed)\n\n"

        s += "| Population | Proportion |\n"
        s += "|------------|------------|\n"
        for pop in sorted(ancestry_results['proportions'],
                          key=lambda p: -ancestry_results['proportions'][p]):
            prop = ancestry_results['proportions'][pop]
            label = POPULATION_LABELS.get(pop, pop)
            s += f"| {label} ({pop}) | {prop:.1%} |\n"

        s += ("\n*Ancestry estimation uses ~55 ancestry-informative markers. "
              "This is a rough superpopulation estimate.*\n")
    else:
        s += "Insufficient markers for ancestry estimation.\n"
    s += "\n---\n\n"
    return s


def _section_prs(prs_results):
    s = "## 7. Polygenic Risk Scores\n\n"
    if not prs_results:
        s += "PRS calculation not available.\n\n---\n\n"
        return s

    s += ("*PRS models estimate relative genetic risk based on common variants. "
          "They do not account for lifestyle, environment, or rare variants.*\n\n")

    any_non_applicable = any(not r['ancestry_applicable'] for r in prs_results.values())
    if any_non_applicable:
        s += ("> **Note:** Your ancestry profile is substantially non-European. "
              "These PRS models may not accurately reflect your risk.\n\n")

    s += "| Condition | Percentile | Risk Category | SNPs Found | Reference |\n"
    s += "|-----------|-----------|---------------|------------|----------|\n"
    for cid, r in prs_results.items():
        flag = "" if r['ancestry_applicable'] else " *"
        s += (f"| {r['name']} | {r['percentile']:.0f}th | "
              f"{r['risk_category'].title()}{flag} | "
              f"{r['snps_found']}/{r['snps_total']} | "
              f"{r['reference']} |\n")

    elevated = [r for r in prs_results.values() if r['risk_category'] in ('elevated', 'high')]
    if elevated:
        s += "\n### Elevated Risk Conditions\n\n"
        for r in elevated:
            s += f"**{r['name']}** — {r['percentile']:.0f}th percentile ({r['risk_category']})\n\n"
            if r['contributing_snps']:
                s += "Top contributing variants:\n"
                for snp in r['contributing_snps'][:5]:
                    s += f"- **{snp['gene']}** ({snp['rsid']}): {snp['copies']} copies of risk allele\n"
                s += "\n"
    s += "\n---\n\n"
    return s


def _section_star_alleles(star_alleles):
    s = "## 8. Pharmacogenomic Star Alleles\n\n"
    if not star_alleles:
        s += "Star allele data not available.\n\n---\n\n"
        return s

    s += ("*CPIC-style star allele calling for drug-metabolizing enzymes. "
          "Share with prescribing physicians.*\n\n")

    s += "| Gene | Diplotype | Phenotype | SNPs Found | Note |\n"
    s += "|------|-----------|-----------|------------|------|\n"
    for gene, r in star_alleles.items():
        phenotype = r['phenotype'].replace('_', ' ').title()
        note = r['clinical_note'][:80] + '...' if len(r['clinical_note']) > 80 else r['clinical_note']
        s += (f"| {gene} | {r['diplotype']} | {phenotype} | "
              f"{r['snps_found']}/{r['snps_total']} | {note} |\n")

    actionable = [
        (g, r) for g, r in star_alleles.items()
        if r['phenotype'] in ('poor', 'intermediate', 'rapid', 'ultrarapid')
    ]
    if actionable:
        s += "\n### Actionable Metabolizer Phenotypes\n\n"
        for gene, r in actionable:
            phenotype = r['phenotype'].replace('_', ' ').title()
            s += f"**{gene}** — {r['diplotype']} ({phenotype} Metabolizer)\n\n"
            s += f"{r['clinical_note']}\n\n"
    s += "\n---\n\n"
    return s


def _section_acmg(acmg):
    s = "## 9. ACMG Secondary Findings\n\n"
    if not acmg:
        s += "ACMG screening not available.\n\n---\n\n"
        return s

    s += (f"*Screened {acmg['genes_screened']} medically actionable genes from the "
          f"ACMG SF v3.2 recommendation list.*\n\n")

    if not acmg['acmg_findings']:
        s += "**No pathogenic/likely pathogenic variants found in ACMG genes.**\n"
    else:
        s += f"**{len(acmg['acmg_findings'])} variant(s) found in {acmg['genes_with_variants']} gene(s):**\n\n"
        for f in acmg['acmg_findings']:
            gene = f.get('gene', 'Unknown')
            condition = _first_trait(f.get('traits'), 'Unknown condition')
            stars = f.get('gold_stars', 0)
            s += f"### {gene} — {condition}\n\n"
            s += f"- **RSID:** {f.get('rsid', '')}\n"
            s += f"- **Genotype:** `{f.get('user_genotype', '')}`\n"
            s += f"- **Confidence:** {stars}/4 stars\n"
            s += f"- **Actionability:** {f.get('acmg_actionability', '')}\n\n"
            s += "---\n\n"

    s += f"\n{acmg['summary']}\n"
    s += "\n---\n\n"
    return s


def _section_epistasis(epistasis_results):
    s = "## 10. Gene-Gene Interactions (Epistasis)\n\n"
    if not epistasis_results:
        s += "No significant gene-gene interactions detected.\n\n---\n\n"
        return s

    s += ("*Combined effects of multiple gene variants that differ from individual effects.*\n\n")
    for interaction in epistasis_results:
        risk_icon = {"high": "!!!", "moderate": "!!", "low": "!"}.get(
            interaction['risk_level'], "")
        s += f"### {interaction['name']} {risk_icon}\n\n"
        s += f"**Risk Level:** {interaction['risk_level'].title()}\n\n"
        s += f"**Genes:** {', '.join(interaction['genes_involved'].keys())}\n\n"
        s += f"**Effect:** {interaction['effect']}\n\n"
        s += f"**Mechanism:** {interaction['mechanism']}\n\n"
        s += "**Recommended Actions:**\n"
        for action in interaction['actions']:
            s += f"- {action}\n"
        s += "\n---\n\n"
    return s


def _section_recommendations(recommendations, findings_dict):
    s = "## 11. Personalized Recommendations\n\n"

    if not recommendations:
        s += "No personalized recommendations generated (module not run).\n\n"
    else:
        # Priority actions
        s += _subsection_recommendation_details(recommendations)

    # Supplement recommendations (always based on lifestyle findings)
    s += _subsection_supplements(findings_dict)
    s += _subsection_diet(findings_dict)
    s += _subsection_lifestyle(findings_dict)

    s += "\n---\n\n"
    return s


def _subsection_recommendation_details(recommendations):
    s = ""
    # Priority actions
    priorities = recommendations.get("priorities", [])
    if priorities:
        s += "*Prioritized synthesis of all genetic signals.*\n\n"
        for p in priorities:
            icon = {"high": "!!!", "moderate": "!!", "low": ""}.get(p["priority"], "")
            s += f"### {p['title']} ({p['priority'].title()}) {icon}\n\n"
            s += f"**Why:** {p['why']}\n\n"
            s += "**Actions:**\n"
            for i, action in enumerate(p["actions"], 1):
                s += f"{i}. {action}\n"
            if p.get("doctor_note"):
                s += f"\n> **Doctor note:** {p['doctor_note']}\n"
            if p.get("monitoring"):
                s += "\n**Monitoring:**\n"
                for m in p["monitoring"]:
                    s += f"- {m['test']} ({m['frequency']}): {m['reason']}\n"
            s += "\n---\n\n"

    # Drug card
    drug_card = recommendations.get("drug_card", [])
    if drug_card:
        s += "### Drug-Gene Card\n\n"
        s += "| Gene | rsID | Genotype | Status | Source |\n"
        s += "|------|------|----------|--------|--------|\n"
        for entry in drug_card:
            for e in entry["entries"]:
                status = e.get("status", "").replace("_", " ").title() or "\u2014"
                s += (f"| {entry['gene']} | {e['rsid']} | "
                      f"`{e['genotype']}` | {status} | {e['source']} |\n")
        s += "\n"

    # Monitoring schedule
    schedule = recommendations.get("monitoring_schedule", [])
    if schedule:
        s += "### Monitoring Schedule\n\n"
        s += "| Test | Frequency | Reason |\n"
        s += "|------|-----------|--------|\n"
        for m in schedule:
            s += f"| {m['test']} | {m['frequency']} | {m['reason']} |\n"
        s += "\n"

    # Good news
    good_news = recommendations.get("good_news", [])
    if good_news:
        s += "### Good News\n\n"
        for g in good_news:
            s += f"- **{g['gene']}**: {g['description']}\n"
        s += "\n"

    return s


def _subsection_supplements(findings_dict):
    s = "### Supplement Recommendations\n\n"
    s += "*Discuss with healthcare provider before starting any supplements*\n\n"

    supplements = []
    if 'MTHFR' in findings_dict and findings_dict['MTHFR']['magnitude'] >= 2:
        supplements.append(("Methylfolate (L-5-MTHF)", "400-800mcg daily",
                            "MTHFR variant reduces folic acid conversion",
                            "Avoid synthetic folic acid. Start low if slow COMT."))
        supplements.append(("Methylcobalamin (B12)", "1000mcg sublingual",
                            "Supports methylation cycle",
                            "Prefer methylcobalamin over cyanocobalamin"))
    if 'MTRR' in findings_dict and findings_dict['MTRR']['magnitude'] >= 2:
        if not any('B12' in s[0] for s in supplements):
            supplements.append(("Methylcobalamin (B12)", "1000-5000mcg sublingual",
                                "MTRR variant impairs B12 recycling",
                                "May need higher doses than typical"))
    if 'GC' in findings_dict and findings_dict['GC'].get('status') == 'low':
        supplements.append(("Vitamin D3", "2000-5000 IU daily",
                            "Genetically low vitamin D binding protein",
                            "Take with fat. Test 25-OH-D after 2-3 months."))
    if 'FADS1' in findings_dict and findings_dict['FADS1'].get('status') == 'low_conversion':
        supplements.append(("Fish Oil (EPA/DHA)", "1-2g daily",
                            "Poor conversion from plant omega-3s",
                            "Direct marine source required."))
    if 'COMT' in findings_dict and findings_dict['COMT'].get('status') == 'slow':
        supplements.append(("Magnesium Glycinate", "300-400mg evening",
                            "Supports COMT function", "Glycinate form preferred."))

    if supplements:
        s += "| Supplement | Dose | Reason | Notes |\n"
        s += "|------------|------|--------|-------|\n"
        for name, dose, reason, notes in supplements:
            s += f"| {name} | {dose} | {reason} | {notes} |\n"
    else:
        s += "No specific supplements indicated by genetic profile.\n"
    s += "\n"
    return s


def _subsection_diet(findings_dict):
    s = "### Dietary Recommendations\n\n"
    recs = []
    if 'APOA2' in findings_dict and findings_dict['APOA2'].get('status') == 'sensitive':
        recs.append("**Limit saturated fat (<7% calories)**: APOA2 variant links sat fat intake to weight gain.")
    if 'MTHFR' in findings_dict and findings_dict['MTHFR']['magnitude'] >= 2:
        recs.append("**Emphasize folate-rich foods**: Leafy greens, legumes, liver.")
    if 'MCM6/LCT' in findings_dict and 'intolerant' in findings_dict['MCM6/LCT'].get('status', ''):
        recs.append("**Lactose intolerance**: May tolerate fermented dairy. Ensure calcium from other sources.")

    caffeine_issues = []
    if 'CYP1A2' in findings_dict and findings_dict['CYP1A2'].get('status') in ['slow', 'intermediate']:
        caffeine_issues.append("slow metabolizer")
    if 'ADORA2A' in findings_dict and findings_dict['ADORA2A'].get('status') == 'anxiety_prone':
        caffeine_issues.append("anxiety-prone")
    if caffeine_issues:
        recs.append(f"**Caffeine caution** ({', '.join(caffeine_issues)}): Limit to morning only.")

    if recs:
        for rec in recs:
            s += f"- {rec}\n\n"
    else:
        s += "No specific dietary modifications indicated.\n"
    s += "\n"
    return s


def _subsection_lifestyle(findings_dict):
    s = "### Lifestyle Recommendations\n\n"
    recs = []
    if 'COMT' in findings_dict and findings_dict['COMT'].get('status') == 'slow':
        recs.append("**Stress management is critical**: Slow COMT means catecholamines build up under stress.")
    if 'BDNF' in findings_dict and findings_dict['BDNF']['magnitude'] >= 2:
        recs.append("**Exercise is essential**: BDNF variant reduces activity-dependent brain growth factor.")
    if 'ACTN3' in findings_dict:
        status = findings_dict['ACTN3'].get('status', '')
        if status == 'endurance':
            recs.append("**Training style (ACTN3 endurance)**: Genetics favor endurance/aerobic training.")
        elif status == 'power':
            recs.append("**Training style (ACTN3 power)**: Genetics favor explosive/strength training.")
    if 'MC1R' in findings_dict:
        recs.append("**Sun protection (MC1R)**: Daily SPF 30+, avoid excessive sun exposure.")

    if recs:
        for rec in recs:
            s += f"- {rec}\n\n"
    else:
        s += "Standard healthy lifestyle recommendations apply.\n"
    s += "\n"
    return s


def _section_lifestyle_findings(findings):
    """Section 12: Complete findings by category."""
    categories = sorted(set(f.get('category', 'Other') for f in findings))

    s = "## 12. Complete Lifestyle/Health Findings by Category\n\n"
    s += f"*{len(findings)} total findings across {len(categories)} categories.*\n\n"

    for category in categories:
        cat_findings = [f for f in findings if f.get('category') == category]
        if not cat_findings:
            continue

        s += f"### {category}\n\n"
        cat_findings.sort(key=lambda x: x.get('magnitude', 0), reverse=True)

        for f in cat_findings:
            gene = f.get('gene', 'Unknown')
            rsid = f.get('rsid', '')
            genotype = f.get('genotype', '')
            status = f.get('status', '').replace('_', ' ').title()
            description = f.get('description', '')
            magnitude = f.get('magnitude', 0)

            s += f"- **{gene}** ({rsid}): `{genotype}` | {status} | Impact {magnitude}/6 — {description}\n"

            context = get_clinical_context(gene, f.get('status', ''))
            if context and context.get('actions'):
                s += f"  - Key action: {context['actions'][0]}\n"

        s += "\n---\n\n"

    return s


def _section_pathway_analysis(findings):
    s = "## 13. Pathway Analysis\n\n"
    s += "*Genes grouped by biological pathway.*\n\n"

    found_any = False
    for pathway_name, pathway_genes in PATHWAYS.items():
        pathway_findings = [f for f in findings if f.get('gene') in pathway_genes]
        if not pathway_findings:
            continue
        found_any = True
        s += f"### {pathway_name}\n\n"
        for f in pathway_findings:
            gene = f.get('gene', '')
            status = f.get('status', '').replace('_', ' ').title()
            magnitude = f.get('magnitude', 0)
            s += f"- **{gene}**: {status} (impact {magnitude}/6)\n"
        s += "\n---\n\n"

    if not found_any:
        s += "No pathway interactions detected.\n\n---\n\n"
    return s


def _section_disease_risk(affected, carriers, het_unknown, disease_findings,
                           disease_stats, genome_count):
    s = "## 14. Disease Risk Analysis\n\n"

    if disease_stats:
        s += f"*ClinVar variants scanned: {disease_stats.get('total_clinvar', 0):,}*\n\n"

    # Pathogenic - affected
    if affected:
        s += "### Pathogenic Variants — Affected Status\n\n"
        for f in affected:
            stars = '*' * f['gold_stars'] + '.' * (4 - f['gold_stars'])
            s += f"#### {f['gene']} — {_first_trait(f['traits'])}\n\n"
            s += f"- **Position:** chr{f['chromosome']}:{f['position']}\n"
            s += f"- **RSID:** {f['rsid']}\n"
            s += f"- **Genotype:** `{f['user_genotype']}` ({'Homozygous' if f['is_homozygous'] else 'Heterozygous'})\n"
            s += f"- **Confidence:** {stars} ({f['gold_stars']}/4)\n\n"

            pw = get_population_warnings(f['gene'], 'pathogenic')
            if pw:
                for w in pw:
                    s += f"> **Population Note:** {w}\n\n"
            s += "---\n\n"

    # Carrier
    if carriers:
        s += "### Carrier Status — Recessive Conditions\n\n"
        s += "**Carrier — no personal symptoms expected, but reproductive implications.**\n\n"
        for f in carriers:
            s += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` — {_first_trait(f['traits'])} ({f['gold_stars']}/4 stars)\n"
        s += "\n---\n\n"

    # Het unknown
    if het_unknown:
        s += "### Pathogenic/Likely Pathogenic — Inheritance Unclear\n\n"
        for f in het_unknown:
            s += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` — {_first_trait(f['traits'])} ({f['gold_stars']}/4 stars)\n"
        s += "\n---\n\n"

    # Risk factors
    if disease_findings and disease_findings.get('risk_factor'):
        s += "### Risk Factor Variants\n\n"
        for f in sorted(disease_findings['risk_factor'], key=lambda x: -x['gold_stars'])[:MAX_RISK_FACTORS]:
            s += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` — {_truncate(f['traits']) if f['traits'] else 'Risk factor'}\n"
        s += "\n---\n\n"

    # Drug response
    if disease_findings and disease_findings.get('drug_response'):
        s += "### Drug Response Variants\n\n"
        for f in disease_findings['drug_response'][:MAX_DRUG_RESPONSES]:
            s += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` — {_truncate(f['traits']) if f['traits'] else 'Drug response'}\n"
        s += "\n---\n\n"

    # Protective
    if disease_findings and disease_findings.get('protective'):
        s += "### Protective Variants\n\n"
        for f in disease_findings['protective']:
            s += f"- **{f['gene']}** ({f['rsid']}): `{f['user_genotype']}` — {_truncate(f['traits']) if f['traits'] else 'Protective'}\n"
        s += "\n---\n\n"

    if not affected and not carriers and not het_unknown:
        if not disease_findings:
            s += "ClinVar data not available.\n\n---\n\n"
        else:
            s += "No pathogenic or likely pathogenic variants detected.\n\n---\n\n"
    return s


def _section_carrier_screening(carrier_screen):
    s = "## 15. Carrier Screening\n\n"
    if not carrier_screen or carrier_screen['total_carriers'] == 0:
        s += "No carrier findings to report.\n\n---\n\n"
        return s

    s += (f"*{carrier_screen['total_carriers']} carrier finding(s) organized by disease system. "
          f"Relevant for reproductive planning.*\n\n")

    for system, carriers in sorted(carrier_screen['by_system'].items()):
        s += f"### {system}\n\n"
        for c in carriers:
            s += f"- **{c['gene']}** ({c['rsid']}): {c['condition']}\n"
            s += f"  - Inheritance: {c['inheritance']}\n"
            if c['reproductive_note']:
                s += f"  - {c['reproductive_note']}\n"
        s += "\n"

    if carrier_screen['couples_relevant']:
        s += "### Couples-Relevant Conditions\n\n"
        s += "These conditions are commonly included in expanded carrier screening panels:\n\n"
        for c in carrier_screen['couples_relevant']:
            s += f"- **{c['gene']}**: {c['condition']}\n"
        s += "\n"

    s += "\n---\n\n"
    return s


def _section_drug_interactions(health_results, disease_findings):
    s = "## 16. Drug-Gene Interactions\n\n"
    s += "**Share this section with prescribing physicians.**\n\n"

    # PharmGKB Level 1
    s += "### PharmGKB Level 1 (Clinical Guidelines)\n\n"
    level_1 = [f for f in health_results['pharmgkb_findings'] if f['level'] in ['1A', '1B']]
    if level_1:
        s += "| Gene | Level | Drugs | Your Genotype |\n"
        s += "|------|-------|-------|---------------|\n"
        for f in level_1:
            drugs = _truncate(f['drugs'], 50)
            s += f"| {f['gene']} | {f['level']} | {drugs} | `{f['genotype']}` |\n"
    else:
        s += "None detected.\n"

    # PharmGKB Level 2
    s += "\n### PharmGKB Level 2 (Moderate Evidence)\n\n"
    level_2 = [f for f in health_results['pharmgkb_findings'] if f['level'] in ['2A', '2B']]
    if level_2:
        s += "| Gene | Level | Drugs | Your Genotype |\n"
        s += "|------|-------|-------|---------------|\n"
        for f in level_2[:MAX_LEVEL2_INTERACTIONS]:
            drugs = _truncate(f['drugs'], 50)
            s += f"| {f['gene']} | {f['level']} | {drugs} | `{f['genotype']}` |\n"
        if len(level_2) > MAX_LEVEL2_INTERACTIONS:
            s += f"\n*...and {len(level_2) - MAX_LEVEL2_INTERACTIONS} more Level 2 interactions*\n"
    else:
        s += "None detected.\n"

    # ClinVar drug response
    s += "\n### ClinVar Drug Response Variants\n\n"
    if disease_findings and disease_findings.get('drug_response'):
        drug_resp = disease_findings['drug_response']
        s += "| Gene | RSID | Genotype | Drug/Response |\n"
        s += "|------|------|----------|---------------|\n"
        for f in drug_resp[:MAX_CLINVAR_DRUG_RESPONSES]:
            traits = _truncate(f['traits'], 60)
            gene = f['gene'] if f['gene'] else '\u2014'
            s += f"| {gene} | {f['rsid']} | `{f['user_genotype']}` | {traits} |\n"
    else:
        s += "None detected.\n"

    s += "\n---\n\n"
    return s


def _section_doctor_card(quality_metrics, blood_type, apoe, star_alleles,
                          acmg, affected, carriers, prs_results):
    s = "## 17. Doctor Card (Print-Friendly Summary)\n\n"
    s += "*Bring this section to medical appointments.*\n\n"

    if blood_type and blood_type.get('blood_type') != 'Unknown':
        s += f"- **Blood Type:** {blood_type['blood_type']}\n"
    if apoe and apoe.get('apoe_type') != 'Unknown':
        s += f"- **APOE:** {apoe['apoe_type']} ({apoe['risk_level']})\n"

    if star_alleles:
        for gene, r in star_alleles.items():
            if r['phenotype'] not in ('normal', 'Unknown'):
                phenotype = r['phenotype'].replace('_', ' ').title()
                s += f"- **{gene}:** {r['diplotype']} ({phenotype} Metabolizer)\n"

    if affected:
        s += "\n**Pathogenic Variants:**\n"
        for f in affected[:5]:
            s += f"- {f['gene']}: {_first_trait(f['traits'])}\n"

    if carriers:
        s += "\n**Carrier Status:**\n"
        for f in carriers[:5]:
            s += f"- {f['gene']}: {_first_trait(f['traits'])}\n"

    if acmg and acmg.get('acmg_findings'):
        s += "\n**ACMG Actionable Findings:**\n"
        for f in acmg['acmg_findings']:
            s += f"- {f.get('gene', 'Unknown')}: {f.get('acmg_actionability', '')}\n"

    if prs_results:
        elevated = [r for r in prs_results.values() if r['risk_category'] in ('elevated', 'high')]
        if elevated:
            s += "\n**Elevated PRS:**\n"
            for r in elevated:
                s += f"- {r['name']}: {r['percentile']:.0f}th percentile ({r['risk_category']})\n"

    s += "\n---\n\n"
    return s


def _section_disclaimer():
    return """## 18. References & Disclaimer

This report is for **informational and educational purposes only**. It is NOT medical advice.

- Genetic associations are probabilistic, not deterministic
- Environmental factors, lifestyle, and other genes also influence outcomes
- Classifications evolve as research progresses
- Consult healthcare providers before making medical decisions
- Some findings may have different implications in different populations

### How to Use This Report

1. **Share with providers** — Especially pharmacogenomics and ACMG sections
2. **Focus on actionable items** — Prioritize evidence-based interventions
3. **Don't over-interpret** — One gene doesn't define your health destiny
4. **Combine with testing** — Many recommendations include follow-up lab tests

### Key References

- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
- PharmGKB: https://www.pharmgkb.org/
- CPIC Guidelines: https://cpicpgx.org/
- ACMG SF v3.2: https://www.acmg.net/

---

*Generated by Genetic Health Analysis Pipeline*
"""
