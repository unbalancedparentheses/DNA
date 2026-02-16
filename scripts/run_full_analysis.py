#!/usr/bin/env python3
"""
Full Genetic Health Analysis Pipeline

This master script runs the complete genetic analysis workflow:
1. Lifestyle/Health Analysis -> EXHAUSTIVE_GENETIC_REPORT.md
2. Disease Risk Analysis -> EXHAUSTIVE_DISEASE_RISK_REPORT.md
3. Actionable Health Protocol -> ACTIONABLE_HEALTH_PROTOCOL.md

Usage:
    python run_full_analysis.py                     # Uses default genome.txt
    python run_full_analysis.py path/to/genome.txt  # Custom genome file
    python run_full_analysis.py --name "John Doe"   # Add name to reports
"""

import sys
import json
import csv
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# Add scripts directory to path for imports
SCRIPT_DIR = Path(__file__).parent
sys.path.insert(0, str(SCRIPT_DIR))

from comprehensive_snp_database import COMPREHENSIVE_SNPS
from report_generators import (
    classify_zygosity,
    generate_exhaustive_genetic_report,
    generate_disease_risk_report,
    generate_actionable_protocol,
)

# Directory configuration
BASE_DIR = SCRIPT_DIR.parent
DATA_DIR = BASE_DIR / "data"
REPORTS_DIR = BASE_DIR / "reports"


def print_header(text):
    """Print a formatted header."""
    print()
    print("=" * 70)
    print(text)
    print("=" * 70)


def print_step(text):
    """Print a step indicator."""
    print(f"\n>>> {text}")


# =============================================================================
# GENOME LOADING
# =============================================================================

def load_genome(genome_path: Path) -> tuple:
    """Load 23andMe genome file into dictionaries."""
    print_step(f"Loading genome from {genome_path}")

    genome_by_rsid = {}
    genome_by_position = {}

    with open(genome_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                rsid, chrom, pos, genotype = parts[0], parts[1], parts[2], parts[3]
                if genotype != '--':
                    genome_by_rsid[rsid] = {
                        'chromosome': chrom,
                        'position': pos,
                        'genotype': genotype
                    }
                    pos_key = f"{chrom}:{pos}"
                    genome_by_position[pos_key] = {
                        'rsid': rsid,
                        'genotype': genotype
                    }

    print(f"    Loaded {len(genome_by_rsid):,} SNPs")
    return genome_by_rsid, genome_by_position


# =============================================================================
# PHARMGKB LOADING
# =============================================================================

def load_pharmgkb() -> dict:
    """Load PharmGKB drug-gene annotations."""
    annotations_path = DATA_DIR / "clinical_annotations.tsv"
    alleles_path = DATA_DIR / "clinical_ann_alleles.tsv"

    if not annotations_path.exists() or not alleles_path.exists():
        print("    PharmGKB files not found, skipping drug interactions")
        return {}

    print_step("Loading PharmGKB data")

    pharmgkb = {}
    annotations = {}

    with open(annotations_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ann_id = row.get('Clinical Annotation ID', '')
            variant = row.get('Variant/Haplotypes', '')
            if variant.startswith('rs'):
                annotations[ann_id] = {
                    'rsid': variant,
                    'gene': row.get('Gene', ''),
                    'drugs': row.get('Drug(s)', ''),
                    'phenotype': row.get('Phenotype(s)', ''),
                    'level': row.get('Level of Evidence', ''),
                    'category': row.get('Phenotype Category', ''),
                }

    with open(alleles_path, 'r') as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            ann_id = row.get('Clinical Annotation ID', '')
            if ann_id in annotations:
                rsid = annotations[ann_id]['rsid']
                genotype = row.get('Genotype/Allele', '')
                if rsid not in pharmgkb:
                    pharmgkb[rsid] = {
                        'gene': annotations[ann_id]['gene'],
                        'drugs': annotations[ann_id]['drugs'],
                        'phenotype': annotations[ann_id]['phenotype'],
                        'level': annotations[ann_id]['level'],
                        'category': annotations[ann_id]['category'],
                        'genotypes': {}
                    }
                pharmgkb[rsid]['genotypes'][genotype] = row.get('Annotation Text', '')

    print(f"    Loaded {len(pharmgkb):,} drug-gene interactions")
    return pharmgkb


# =============================================================================
# LIFESTYLE/HEALTH ANALYSIS
# =============================================================================

def analyze_lifestyle_health(genome_by_rsid: dict, pharmgkb: dict) -> dict:
    """Analyze genome against lifestyle/health SNP database."""
    print_step("Running lifestyle/health analysis")

    results = {
        'findings': [],
        'pharmgkb_findings': [],
        'by_category': defaultdict(list),
        'summary': {
            'total_snps': len(genome_by_rsid),
            'analyzed_snps': 0,
            'high_impact': 0,
            'moderate_impact': 0,
            'low_impact': 0,
        }
    }

    # Check against comprehensive database
    for rsid, info in COMPREHENSIVE_SNPS.items():
        if rsid in genome_by_rsid:
            genotype = genome_by_rsid[rsid]['genotype']
            genotype_rev = genotype[::-1] if len(genotype) == 2 else genotype

            variant_info = info['variants'].get(genotype) or info['variants'].get(genotype_rev)

            if variant_info:
                finding = {
                    'rsid': rsid,
                    'gene': info['gene'],
                    'category': info['category'],
                    'genotype': genotype,
                    'status': variant_info['status'],
                    'description': variant_info['desc'],
                    'magnitude': variant_info['magnitude'],
                    'note': info.get('note', ''),
                }
                results['findings'].append(finding)
                results['by_category'][info['category']].append(finding)
                results['summary']['analyzed_snps'] += 1

                if variant_info['magnitude'] >= 3:
                    results['summary']['high_impact'] += 1
                elif variant_info['magnitude'] >= 2:
                    results['summary']['moderate_impact'] += 1
                elif variant_info['magnitude'] >= 1:
                    results['summary']['low_impact'] += 1

    # Check PharmGKB
    for rsid, info in pharmgkb.items():
        if rsid in genome_by_rsid:
            genotype = genome_by_rsid[rsid]['genotype']
            genotype_rev = genotype[::-1] if len(genotype) == 2 else genotype
            annotation = info['genotypes'].get(genotype) or info['genotypes'].get(genotype_rev)
            if annotation and info['level'] in ['1A', '1B', '2A', '2B']:
                finding = {
                    'rsid': rsid,
                    'gene': info['gene'],
                    'drugs': info['drugs'],
                    'genotype': genotype,
                    'annotation': annotation,
                    'level': info['level'],
                    'category': info['category'],
                }
                results['pharmgkb_findings'].append(finding)

    # Sort findings
    results['findings'].sort(key=lambda x: -x['magnitude'])
    results['pharmgkb_findings'].sort(key=lambda x: x['level'])

    print(f"    Found {len(results['findings'])} lifestyle/health findings")
    print(f"    Found {len(results['pharmgkb_findings'])} drug-gene interactions")

    return results


# =============================================================================
# DISEASE RISK ANALYSIS
# =============================================================================

def load_clinvar_and_analyze(genome_by_position: dict) -> tuple:
    """Load ClinVar and analyze for disease variants."""
    clinvar_path = DATA_DIR / "clinvar_alleles.tsv"

    if not clinvar_path.exists():
        print("    ClinVar file not found, skipping disease risk analysis")
        return None, None

    print_step("Loading ClinVar and analyzing disease risk")

    findings = {
        'pathogenic': [],
        'likely_pathogenic': [],
        'risk_factor': [],
        'drug_response': [],
        'protective': [],
        'other_significant': []
    }

    stats = {
        'total_clinvar': 0,
        'matched': 0,
        'pathogenic_matched': 0,
        'likely_pathogenic_matched': 0
    }

    with open(clinvar_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f, delimiter='\t')

        for row in reader:
            stats['total_clinvar'] += 1

            chrom = row['chrom']
            pos = row['pos']
            pos_key = f"{chrom}:{pos}"

            if pos_key not in genome_by_position:
                continue

            stats['matched'] += 1

            user_data = genome_by_position[pos_key]
            user_genotype = user_data['genotype']
            ref_allele = row['ref']
            alt_allele = row['alt']
            clinical_sig = row['clinical_significance'].lower()

            # Only process true SNPs
            if len(ref_allele) != 1 or len(alt_allele) != 1:
                continue

            has_variant = alt_allele in user_genotype
            is_homozygous = user_genotype == alt_allele + alt_allele
            is_heterozygous = has_variant and not is_homozygous
            has_ref_only = user_genotype == ref_allele + ref_allele

            if has_ref_only or not has_variant:
                continue

            finding = {
                'chromosome': chrom,
                'position': pos,
                'rsid': user_data['rsid'],
                'gene': row['symbol'],
                'ref': ref_allele,
                'alt': alt_allele,
                'user_genotype': user_genotype,
                'is_homozygous': is_homozygous,
                'is_heterozygous': is_heterozygous,
                'clinical_significance': row['clinical_significance'],
                'review_status': row['review_status'],
                'gold_stars': int(row['gold_stars']) if row['gold_stars'] else 0,
                'traits': row['all_traits'],
                'inheritance': row.get('inheritance_modes', ''),
                'hgvs_p': row.get('hgvs_p', ''),
                'hgvs_c': row.get('hgvs_c', ''),
                'molecular_consequence': row.get('molecular_consequence', ''),
                'xrefs': row.get('xrefs', '')
            }

            if 'pathogenic' in clinical_sig and 'likely' not in clinical_sig and 'conflict' not in clinical_sig:
                findings['pathogenic'].append(finding)
                stats['pathogenic_matched'] += 1
            elif 'likely pathogenic' in clinical_sig or 'likely_pathogenic' in clinical_sig:
                findings['likely_pathogenic'].append(finding)
                stats['likely_pathogenic_matched'] += 1
            elif 'risk factor' in clinical_sig or 'risk_factor' in clinical_sig:
                findings['risk_factor'].append(finding)
            elif 'drug response' in clinical_sig or 'drug_response' in clinical_sig:
                findings['drug_response'].append(finding)
            elif 'protective' in clinical_sig:
                findings['protective'].append(finding)
            elif 'association' in clinical_sig or 'affects' in clinical_sig:
                findings['other_significant'].append(finding)

    print(f"    ClinVar entries scanned: {stats['total_clinvar']:,}")
    print(f"    Pathogenic variants: {stats['pathogenic_matched']}")
    print(f"    Likely pathogenic: {stats['likely_pathogenic_matched']}")
    print(f"    Risk factors: {len(findings['risk_factor'])}")

    return findings, stats




# =============================================================================
# MAIN PIPELINE
# =============================================================================

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

    # Save intermediate results for exhaustive report generator
    results_json = {
        'findings': health_results['findings'],
        'pharmgkb_findings': health_results['pharmgkb_findings'],
        'summary': health_results['summary'],
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
    generate_actionable_protocol(health_results, disease_findings, protocol_path, subject_name)

    # Summary
    print_header("ANALYSIS COMPLETE")
    print(f"\nReports generated in: {REPORTS_DIR}  (Markdown + HTML)")
    print(f"\n  1. EXHAUSTIVE_GENETIC_REPORT.{{md,html}}")
    print(f"     - {len(health_results['findings'])} lifestyle/health findings")
    print(f"     - {len(health_results['pharmgkb_findings'])} drug-gene interactions")

    if disease_findings:
        total_disease = (len(disease_findings['pathogenic']) +
                        len(disease_findings['likely_pathogenic']) +
                        len(disease_findings['risk_factor']))
        print(f"\n  2. EXHAUSTIVE_DISEASE_RISK_REPORT.{{md,html}}")
        print(f"     - {len(disease_findings['pathogenic'])} pathogenic variants")
        print(f"     - {len(disease_findings['likely_pathogenic'])} likely pathogenic")
        print(f"     - {len(disease_findings['risk_factor'])} risk factors")

    print(f"\n  3. ACTIONABLE_HEALTH_PROTOCOL_V3.{{md,html}}")
    print(f"     - Comprehensive protocol (lifestyle + disease risk + carrier status)")

    print(f"\nFinished: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

    return {
        'health_results': health_results,
        'disease_findings': disease_findings,
        'disease_stats': disease_stats
    }


def main():
    """Main entry point."""
    import argparse

    parser = argparse.ArgumentParser(
        description="Run full genetic health analysis pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_full_analysis.py                        # Use default genome.txt
  python run_full_analysis.py /path/to/genome.txt   # Custom genome file
  python run_full_analysis.py --name "John Doe"     # Add name to reports
        """
    )
    parser.add_argument('genome', nargs='?', type=Path, default=None,
                       help='Path to 23andMe genome file (default: data/genome.txt)')
    parser.add_argument('--name', '-n', type=str, default=None,
                       help='Subject name to include in reports')

    args = parser.parse_args()

    run_full_analysis(args.genome, args.name)


if __name__ == "__main__":
    main()
