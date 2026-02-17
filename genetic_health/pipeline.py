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
from .reports import (
    generate_exhaustive_genetic_report,
    generate_disease_risk_report,
    generate_actionable_protocol,
)


def print_header(text):
    """Print a formatted header."""
    print()
    print("=" * 70)
    print(text)
    print("=" * 70)


def print_step(text):
    """Print a step indicator."""
    print(f"\n>>> {text}")


def run_full_analysis(genome_path: Path = None, subject_name: str = None):
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

    # Load PharmGKB
    pharmgkb = load_pharmgkb()

    # Run lifestyle/health analysis
    health_results = analyze_lifestyle_health(genome_by_rsid, pharmgkb)

    # Run ancestry estimation
    print_step("Estimating ancestry from AIMs")
    ancestry_results = estimate_ancestry(genome_by_rsid)
    print(f"    Top ancestry: {ancestry_results['top_ancestry']} "
          f"({ancestry_results['markers_found']} markers, "
          f"{ancestry_results['confidence']} confidence)")

    # Run polygenic risk scores
    print_step("Calculating polygenic risk scores")
    prs_results = calculate_prs(genome_by_rsid, ancestry_results['proportions'])
    for cid, r in prs_results.items():
        print(f"    {r['name']}: {r['percentile']:.0f}th percentile "
              f"({r['risk_category']}) [{r['snps_found']}/{r['snps_total']} SNPs]")

    # Save intermediate results for exhaustive report generator
    results_json = {
        'findings': health_results['findings'],
        'pharmgkb_findings': health_results['pharmgkb_findings'],
        'summary': health_results['summary'],
        'ancestry': ancestry_results,
        'prs': prs_results,
    }
    intermediate_path = REPORTS_DIR / "comprehensive_results.json"
    with open(intermediate_path, 'w') as f:
        json.dump(results_json, f, indent=2)

    # Generate exhaustive genetic report
    genetic_report_path = REPORTS_DIR / "EXHAUSTIVE_GENETIC_REPORT.md"
    generate_exhaustive_genetic_report(health_results, genetic_report_path, subject_name)

    # Run disease risk analysis
    disease_findings, disease_stats = load_clinvar_and_analyze(genome_by_position)

    # Generate disease risk report
    if disease_findings:
        disease_report_path = REPORTS_DIR / "EXHAUSTIVE_DISEASE_RISK_REPORT.md"
        generate_disease_risk_report(disease_findings, disease_stats, len(genome_by_rsid),
                                      disease_report_path, subject_name)

    # Generate actionable protocol - use versioned filename
    protocol_path = REPORTS_DIR / "ACTIONABLE_HEALTH_PROTOCOL_V3.md"
    generate_actionable_protocol(health_results, disease_findings, protocol_path, subject_name,
                                 ancestry_results=ancestry_results, prs_results=prs_results)

    # Generate enhanced all-in-one HTML report
    from .reports.enhanced_html import main as generate_enhanced_report
    print_step("Generating enhanced all-in-one HTML report")
    generate_enhanced_report()

    # Summary
    print_header("ANALYSIS COMPLETE")
    print(f"\nReports generated in: {REPORTS_DIR}  (Markdown + HTML)")
    print(f"\n  1. EXHAUSTIVE_GENETIC_REPORT.{{md,html}}")
    print(f"     - {len(health_results['findings'])} lifestyle/health findings")
    print(f"     - {len(health_results['pharmgkb_findings'])} drug-gene interactions")

    if disease_findings:
        print(f"\n  2. EXHAUSTIVE_DISEASE_RISK_REPORT.{{md,html}}")
        print(f"     - {len(disease_findings['pathogenic'])} pathogenic variants")
        print(f"     - {len(disease_findings['likely_pathogenic'])} likely pathogenic")
        print(f"     - {len(disease_findings['risk_factor'])} risk factors")

    print(f"\n  3. ACTIONABLE_HEALTH_PROTOCOL_V3.{{md,html}}")
    print(f"     - Comprehensive protocol (lifestyle + disease risk + carrier status)")

    print(f"\n  4. ENHANCED_HEALTH_REPORT.html")
    print(f"     - All-in-one interactive report (SVG charts, collapsible sections)")

    print(f"\n  Ancestry: {ancestry_results['top_ancestry']} "
          f"({ancestry_results['confidence']} confidence, "
          f"{ancestry_results['markers_found']} markers)")

    elevated_prs = [r for r in prs_results.values() if r['risk_category'] in ('elevated', 'high')]
    if elevated_prs:
        print(f"\n  PRS Alerts:")
        for r in elevated_prs:
            print(f"     - {r['name']}: {r['percentile']:.0f}th percentile ({r['risk_category']})")
    else:
        print(f"\n  PRS: All conditions within average range")

    print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    return {
        'health_results': health_results,
        'disease_findings': disease_findings,
        'disease_stats': disease_stats,
        'ancestry_results': ancestry_results,
        'prs_results': prs_results,
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

    args = parser.parse_args()

    run_full_analysis(args.genome, args.name)
