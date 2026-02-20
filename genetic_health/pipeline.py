"""Pipeline orchestrator and CLI entry point."""

import sys
import json
from pathlib import Path
from datetime import datetime

from .config import DATA_DIR, REPORTS_DIR
from .loading import load_genome, load_pharmgkb
from .analysis import analyze_lifestyle_health, load_clinvar_and_analyze
from .ancestry import estimate_ancestry
from .prs import calculate_prs
from .epistasis import evaluate_epistasis
from .recommendations import generate_recommendations
from .insights import generate_insights
from .quality_metrics import compute_quality_metrics
from .blood_type import predict_blood_type
from .mt_haplogroup import estimate_mt_haplogroup
from .star_alleles import call_star_alleles
from .apoe import call_apoe_haplotype
from .acmg import flag_acmg_findings
from .carrier_screen import organize_carrier_findings
from .traits import predict_traits
from .polypharmacy import assess_polypharmacy
from .longevity import profile_longevity
from .sleep_profile import profile_sleep
from .nutrigenomics import profile_nutrigenomics
from .mental_health import profile_mental_health
from .drug_dosing import generate_drug_dosing
from .preventive_care import generate_preventive_timeline
from .reports import generate_html_report


def print_header(text):
    """Print a formatted header."""
    print()
    print("=" * 70)
    print(text)
    print("=" * 70)


def print_step(text):
    """Print a step indicator."""
    print(f"\n>>> {text}")


def run_full_analysis(genome_path: Path = None, subject_name: str = None,
                      skip_ancestry: bool = False, skip_prs: bool = False,
                      export_pdf: bool = False):
    """Run the complete genetic analysis pipeline."""

    print_header("FULL GENETIC HEALTH ANALYSIS")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    # Default genome path
    if genome_path is None:
        genome_path = DATA_DIR / "genome.txt"

    if not genome_path.exists():
        print(f"\nERROR: Genome file not found: {genome_path}")
        print("Please provide a valid 23andMe genome file.")
        sys.exit(1)

    # Create reports directory
    REPORTS_DIR.mkdir(exist_ok=True)

    # Load genome
    genome_by_rsid, genome_by_position = load_genome(genome_path)

    # Compute data quality metrics
    print_step("Computing data quality metrics")
    quality_metrics = compute_quality_metrics(genome_by_rsid, genome_path)
    print(f"    Total SNPs: {quality_metrics['total_snps']:,}")
    print(f"    Call rate: {quality_metrics['call_rate']:.1%}")
    print(f"    Autosomal: {quality_metrics['autosomal_count']:,}, "
          f"MT: {quality_metrics['mt_snp_count']}, "
          f"Het rate: {quality_metrics['het_rate']:.3f}")
    if quality_metrics['has_y']:
        print("    Sex chromosomes: X + Y (biological male)")
    else:
        print("    Sex chromosomes: X only (biological female)")

    # Predict blood type
    print_step("Predicting blood type")
    blood_type = predict_blood_type(genome_by_rsid)
    print(f"    Blood type: {blood_type['blood_type']} "
          f"({blood_type['confidence']} confidence)")

    # Estimate mitochondrial haplogroup
    print_step("Estimating mitochondrial haplogroup")
    mt_haplogroup = estimate_mt_haplogroup(genome_by_rsid)
    print(f"    Haplogroup: {mt_haplogroup['haplogroup']} "
          f"({mt_haplogroup['confidence']} confidence, "
          f"{mt_haplogroup['markers_found']}/{mt_haplogroup['markers_tested']} markers)")
    print(f"    Lineage: {mt_haplogroup['lineage']}")

    # Call pharmacogenomic star alleles (CYP2C19, CYP2C9, CYP2D6, DPYD, TPMT, UGT1A1)
    print_step("Calling pharmacogenomic star alleles")
    star_alleles = call_star_alleles(genome_by_rsid)
    for gene, result in star_alleles.items():
        print(f"    {gene}: {result['diplotype']} -> {result['phenotype']} metabolizer "
              f"[{result['snps_found']}/{result['snps_total']} SNPs]")

    # Call APOE haplotype
    print_step("Calling APOE haplotype")
    apoe = call_apoe_haplotype(genome_by_rsid)
    print(f"    APOE: {apoe['apoe_type']} ({apoe['risk_level']} risk, "
          f"{apoe['confidence']} confidence)")

    # Predict traits
    print_step("Predicting visible traits")
    traits = predict_traits(genome_by_rsid)
    for trait_name, trait in traits.items():
        print(f"    {trait_name.replace('_', ' ').title()}: {trait['prediction']} "
              f"({trait['confidence']} confidence)")

    # Load PharmGKB
    pharmgkb = load_pharmgkb()

    # Run lifestyle/health analysis
    health_results = analyze_lifestyle_health(genome_by_rsid, pharmgkb)

    # Run ancestry estimation
    if skip_ancestry:
        print_step("Skipping ancestry estimation (--no-ancestry)")
        ancestry_results = None
    else:
        print_step("Estimating ancestry from AIMs")
        ancestry_results = estimate_ancestry(genome_by_rsid)
        print(f"    Top ancestry: {ancestry_results['top_ancestry']} "
              f"({ancestry_results['markers_found']} markers, "
              f"{ancestry_results['confidence']} confidence)")

    # Run polygenic risk scores
    if skip_prs:
        print_step("Skipping PRS calculation (--no-prs)")
        prs_results = None
    else:
        print_step("Calculating polygenic risk scores")
        ancestry_props = ancestry_results['proportions'] if ancestry_results else None
        prs_results = calculate_prs(genome_by_rsid, ancestry_props)
        for cid, r in prs_results.items():
            print(f"    {r['name']}: {r['percentile']:.0f}th percentile "
                  f"({r['risk_category']}) [{r['snps_found']}/{r['snps_total']} SNPs]")

    # Evaluate gene-gene interactions
    print_step("Evaluating gene-gene interactions (epistasis)")
    epistasis_results = evaluate_epistasis(health_results['findings'])
    if epistasis_results:
        print(f"    {len(epistasis_results)} interactions detected")
        for e in epistasis_results:
            print(f"    - {e['name']} ({e['risk_level']})")
    else:
        print("    No significant gene-gene interactions detected")

    # Run disease risk analysis (needed for recommendations)
    disease_findings, disease_stats = load_clinvar_and_analyze(genome_by_position)

    # ACMG secondary findings
    print_step("Screening ACMG SF v3.2 medically actionable genes")
    acmg = flag_acmg_findings(disease_findings)
    print(f"    {acmg['genes_screened']} genes screened, "
          f"{acmg['genes_with_variants']} with variants")

    # Carrier screening
    print_step("Organizing carrier findings")
    carrier_screen = organize_carrier_findings(disease_findings)
    print(f"    {carrier_screen['total_carriers']} carrier findings, "
          f"{len(carrier_screen['couples_relevant'])} couples-relevant")

    # Generate personalized recommendations synthesis
    print_step("Generating personalized recommendations")
    recommendations = generate_recommendations(
        health_results['findings'],
        disease_findings=disease_findings,
        ancestry_results=ancestry_results,
        prs_results=prs_results,
        epistasis_results=epistasis_results,
        star_alleles=star_alleles,
        acmg=acmg,
    )
    high_count = sum(1 for p in recommendations['priorities'] if p['priority'] == 'high')
    mod_count = sum(1 for p in recommendations['priorities'] if p['priority'] == 'moderate')
    print(f"    {len(recommendations['priorities'])} priority areas "
          f"({high_count} high, {mod_count} moderate)")
    ref_count = len(recommendations.get('specialist_referrals', []))
    if ref_count:
        print(f"    {ref_count} specialist referral(s)")

    # Generate research-backed insights
    print_step("Generating research-backed insights")
    insights = generate_insights(
        health_results['findings'],
        apoe=apoe,
        star_alleles=star_alleles,
        ancestry_results=ancestry_results,
        epistasis_results=epistasis_results,
        disease_findings=disease_findings,
    )
    print(f"    {len(insights['single_gene'])} gene insights, "
          f"{len(insights['narratives'])} narratives, "
          f"{len(insights['genome_highlights'])} highlights")

    # Polypharmacy risk assessment
    print_step("Assessing polypharmacy risks")
    polypharmacy = assess_polypharmacy(
        genome_by_rsid, star_alleles=star_alleles,
        lifestyle_findings=health_results['findings'])
    if polypharmacy['total_warnings']:
        print(f"    {polypharmacy['total_warnings']} polypharmacy warnings")
        for w in polypharmacy['warnings']:
            print(f"    - [{w['severity'].upper()}] {w['name']}")
    else:
        print("    No polypharmacy risks detected")

    # Longevity profiling
    print_step("Profiling longevity and healthspan")
    longevity = profile_longevity(
        genome_by_rsid, lifestyle_findings=health_results['findings'],
        apoe=apoe, prs_results=prs_results)
    print(f"    Longevity score: {longevity['longevity_score']}/100 "
          f"({longevity['alleles_checked']} longevity alleles checked)")

    # Sleep/circadian profiling
    print_step("Profiling sleep and circadian rhythm")
    sleep = profile_sleep(genome_by_rsid, lifestyle_findings=health_results['findings'])
    print(f"    Chronotype: {sleep['chronotype']} "
          f"({sleep['confidence']} confidence, {sleep['markers_found']} markers)")
    print(f"    Optimal sleep: {sleep['optimal_sleep_window']}")

    # Nutrigenomics profiling
    print_step("Profiling nutrigenomics")
    nutrigenomics = profile_nutrigenomics(
        genome_by_rsid, lifestyle_findings=health_results['findings'])
    high_needs = [n['name'] for n in nutrigenomics['nutrient_needs'] if n['need_level'] == 'high']
    if high_needs:
        print(f"    High nutrient needs: {', '.join(high_needs)}")
    else:
        print("    No high-priority nutrient deficiency risks")

    # Mental health profiling
    print_step("Profiling mental health genetics")
    mental_health = profile_mental_health(
        genome_by_rsid, lifestyle_findings=health_results['findings'],
        star_alleles=star_alleles)
    elevated = [d for d, info in mental_health['domains'].items()
                if info['risk_level'] == 'elevated']
    if elevated:
        print(f"    Elevated domains: {', '.join(elevated)}")
    else:
        print("    No elevated mental health genetic signals")

    # Drug-specific dosing recommendations
    print_step("Generating drug-specific dosing recommendations")
    drug_dosing = generate_drug_dosing(
        star_alleles, lifestyle_findings=health_results['findings'])
    n_recs = len(drug_dosing['recommendations'])
    n_warn = len(drug_dosing['warnings'])
    if n_warn:
        print(f"    {n_recs} dosing recommendations, {n_warn} critical warning(s)")
    elif n_recs:
        print(f"    {n_recs} dosing recommendations")
    else:
        print("    No drug dosing adjustments needed")

    # Preventive care timeline
    print_step("Generating preventive care timeline")
    preventive_care = generate_preventive_timeline(
        prs_results=prs_results, apoe=apoe, acmg=acmg,
        star_alleles=star_alleles, carrier_screen=carrier_screen)
    print(f"    {len(preventive_care['timeline'])} screening recommendations "
          f"({preventive_care['early_screenings']} genetically modified)")

    # Save intermediate results for HTML report generator
    results_json = {
        'findings': health_results['findings'],
        'pharmgkb_findings': health_results['pharmgkb_findings'],
        'summary': health_results['summary'],
        'ancestry': ancestry_results,
        'prs': prs_results,
        'epistasis': [
            {k: v for k, v in e.items() if k != 'genes_involved'}
            for e in epistasis_results
        ],
        'recommendations': recommendations,
        'quality_metrics': quality_metrics,
        'blood_type': blood_type,
        'mt_haplogroup': mt_haplogroup,
        'star_alleles': star_alleles,
        'apoe': apoe,
        'acmg': acmg,
        'carrier_screen': carrier_screen,
        'traits': traits,
        'insights': insights,
        'disease_findings': disease_findings,
        'polypharmacy': polypharmacy,
        'longevity': longevity,
        'sleep_profile': sleep,
        'nutrigenomics': nutrigenomics,
        'mental_health': mental_health,
        'drug_dosing': drug_dosing,
        'preventive_care': preventive_care,
    }
    intermediate_path = REPORTS_DIR / "comprehensive_results.json"
    with open(intermediate_path, 'w') as f:
        json.dump(results_json, f, indent=2)

    # Generate interactive HTML report
    print_step("Generating interactive HTML report")
    generate_html_report()

    # PDF export
    if export_pdf:
        from .reports.pdf_export import export_pdf as _export_pdf
        print_step("Exporting reports to PDF")
        html_path = REPORTS_DIR / "GENETIC_HEALTH_REPORT.html"
        if html_path.exists():
            pdf_path = _export_pdf(html_path)
            if pdf_path:
                print(f"    PDF: {pdf_path}")

    # Summary
    print_header("ANALYSIS COMPLETE")
    print(f"\nReports generated in: {REPORTS_DIR}")
    print(f"\n  GENETIC_HEALTH_REPORT.html")
    print(f"     - Interactive HTML report (SVG charts, collapsible sections)")
    print(f"     - {len(health_results['findings'])} lifestyle/health findings")
    print(f"     - {len(health_results['pharmgkb_findings'])} drug-gene interactions")

    if disease_findings:
        print(f"     - {len(disease_findings['pathogenic'])} pathogenic variants")
        print(f"     - {len(disease_findings['likely_pathogenic'])} likely pathogenic")
        print(f"     - {len(disease_findings['risk_factor'])} risk factors")

    if ancestry_results:
        print(f"\n  Ancestry: {ancestry_results['top_ancestry']} "
              f"({ancestry_results['confidence']} confidence, "
              f"{ancestry_results['markers_found']} markers)")

    if prs_results:
        elevated_prs = [r for r in prs_results.values() if r['risk_category'] in ('elevated', 'high')]
        if elevated_prs:
            print(f"\n  PRS Alerts:")
            for r in elevated_prs:
                print(f"     - {r['name']}: {r['percentile']:.0f}th percentile ({r['risk_category']})")
        else:
            print(f"\n  PRS: All conditions within average range")

    print(f"\n  APOE: {apoe['apoe_type']} ({apoe['risk_level']})")
    if acmg['genes_with_variants'] > 0:
        print(f"\n  ACMG: {acmg['genes_with_variants']} medically actionable gene(s) found")

    print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    return {
        'health_results': health_results,
        'disease_findings': disease_findings,
        'disease_stats': disease_stats,
        'ancestry_results': ancestry_results,
        'prs_results': prs_results,
        'epistasis_results': epistasis_results,
        'recommendations': recommendations,
        'quality_metrics': quality_metrics,
        'blood_type': blood_type,
        'mt_haplogroup': mt_haplogroup,
        'star_alleles': star_alleles,
        'apoe': apoe,
        'acmg': acmg,
        'carrier_screen': carrier_screen,
        'traits': traits,
        'insights': insights,
        'polypharmacy': polypharmacy,
        'longevity': longevity,
        'sleep_profile': sleep,
        'nutrigenomics': nutrigenomics,
        'mental_health': mental_health,
        'drug_dosing': drug_dosing,
        'preventive_care': preventive_care,
    }


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Run full genetic health analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python -m genetic_health                        # Use default genome.txt
  python -m genetic_health /path/to/genome.txt   # Custom genome file
  python -m genetic_health --name "John Doe"     # Add name to reports
        """
    )
    parser.add_argument('genome', nargs='?', type=Path, default=None,
                       help='Path to 23andMe genome file (default: data/genome.txt)')
    parser.add_argument('--name', '-n', type=str, default=None,
                       help='Subject name to include in reports')
    parser.add_argument('--no-ancestry', action='store_true',
                       help='Skip ancestry estimation')
    parser.add_argument('--no-prs', action='store_true',
                       help='Skip polygenic risk score calculation')
    parser.add_argument('--pdf', action='store_true',
                       help='Export reports to PDF (requires Chrome or weasyprint)')

    args = parser.parse_args()

    run_full_analysis(args.genome, args.name,
                      skip_ancestry=args.no_ancestry, skip_prs=args.no_prs,
                      export_pdf=args.pdf)
