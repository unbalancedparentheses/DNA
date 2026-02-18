"""Lifestyle/health and disease risk analysis logic."""

import csv
from pathlib import Path
from collections import defaultdict

from .config import DATA_DIR
from .snp_database import COMPREHENSIVE_SNPS


def _lookup_genotype(variants_dict, genotype):
    """Look up genotype in a variants dict, trying reverse complement."""
    genotype_rev = genotype[::-1] if len(genotype) == 2 else genotype
    return variants_dict.get(genotype) or variants_dict.get(genotype_rev)


def _safe_int(value, default=0):
    """Convert to int, returning default on failure."""
    try:
        return int(value) if value else default
    except (ValueError, TypeError):
        return default


def analyze_lifestyle_health(genome_by_rsid: dict, pharmgkb: dict) -> dict:
    """Analyze genome against lifestyle/health SNP database."""
    print("\n>>> Running lifestyle/health analysis")

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

    for rsid, info in COMPREHENSIVE_SNPS.items():
        if rsid in genome_by_rsid:
            genotype = genome_by_rsid[rsid]['genotype']
            variant_info = _lookup_genotype(info['variants'], genotype)

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
                    'freq': info.get('freq'),
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

    for rsid, info in pharmgkb.items():
        if rsid in genome_by_rsid:
            genotype = genome_by_rsid[rsid]['genotype']
            annotation = _lookup_genotype(info['genotypes'], genotype)
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

    results['findings'].sort(key=lambda x: -x['magnitude'])
    results['pharmgkb_findings'].sort(key=lambda x: x['level'])

    print(f"    Found {len(results['findings'])} lifestyle/health findings")
    print(f"    Found {len(results['pharmgkb_findings'])} drug-gene interactions")

    return results


def load_clinvar_and_analyze(genome_by_position: dict, data_dir: Path = None) -> tuple:
    """Load ClinVar and analyze for disease variants."""
    if data_dir is None:
        data_dir = DATA_DIR

    clinvar_path = data_dir / "clinvar_alleles.tsv"

    if not clinvar_path.exists():
        print("    ClinVar file not found, skipping disease risk analysis")
        return None, None

    print("\n>>> Loading ClinVar and analyzing disease risk")

    findings = {
        'pathogenic': [],
        'likely_pathogenic': [],
        'risk_factor': [],
        'drug_response': [],
        'protective': [],
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

            chrom = row.get('chrom', '')
            pos = row.get('pos', '')
            if not chrom or not pos:
                continue
            pos_key = f"{chrom}:{pos}"

            if pos_key not in genome_by_position:
                continue

            stats['matched'] += 1

            user_data = genome_by_position[pos_key]
            user_genotype = user_data['genotype']
            ref_allele = row.get('ref', '')
            alt_allele = row.get('alt', '')
            clinical_sig = row.get('clinical_significance', '').lower().replace('_', ' ')

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
                'gene': row.get('symbol', ''),
                'ref': ref_allele,
                'alt': alt_allele,
                'user_genotype': user_genotype,
                'is_homozygous': is_homozygous,
                'is_heterozygous': is_heterozygous,
                'clinical_significance': row.get('clinical_significance', ''),
                'review_status': row.get('review_status', ''),
                'gold_stars': _safe_int(row.get('gold_stars')),
                'traits': row.get('all_traits', ''),
                'inheritance': row.get('inheritance_modes', ''),
                'hgvs_p': row.get('hgvs_p', ''),
                'hgvs_c': row.get('hgvs_c', ''),
                'molecular_consequence': row.get('molecular_consequence', ''),
                'xrefs': row.get('xrefs', '')
            }

            if 'pathogenic' in clinical_sig and 'likely' not in clinical_sig and 'conflict' not in clinical_sig:
                findings['pathogenic'].append(finding)
                stats['pathogenic_matched'] += 1
            elif 'likely pathogenic' in clinical_sig:
                findings['likely_pathogenic'].append(finding)
                stats['likely_pathogenic_matched'] += 1
            elif 'risk factor' in clinical_sig:
                findings['risk_factor'].append(finding)
            elif 'drug response' in clinical_sig:
                findings['drug_response'].append(finding)
            elif 'protective' in clinical_sig:
                findings['protective'].append(finding)

    print(f"    ClinVar entries scanned: {stats['total_clinvar']:,}")
    print(f"    Pathogenic variants: {stats['pathogenic_matched']}")
    print(f"    Likely pathogenic: {stats['likely_pathogenic_matched']}")
    print(f"    Risk factors: {len(findings['risk_factor'])}")

    return findings, stats
