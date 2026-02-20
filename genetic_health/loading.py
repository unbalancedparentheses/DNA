"""Genome and PharmGKB data loading."""

import csv
import json
from pathlib import Path

from .config import DATA_DIR

_VALID_BASES = frozenset("ACGT")


def _load_position_to_rsid_map():
    """Build a chrom:pos -> rsID reverse lookup from rsid_positions_grch37.json.

    This allows WGS-derived genome files (which use positional IDs like
    chr1_871334 instead of rsIDs) to be matched to known rsIDs.
    """
    lookup_path = DATA_DIR / "rsid_positions_grch37.json"
    if not lookup_path.exists():
        return {}
    with open(lookup_path, 'r') as f:
        rsid_to_pos = json.load(f)
    pos_to_rsid = {}
    for rsid, info in rsid_to_pos.items():
        key = f"{info['chrom']}:{info['pos']}"
        pos_to_rsid[key] = rsid
    return pos_to_rsid


def load_genome(genome_path: Path) -> tuple:
    """Load 23andMe or WGS-converted genome file into dictionaries.

    Handles both rsID-based entries (rs12345) and positional entries
    (chr1_871334). Positional entries are matched to known rsIDs using
    the rsid_positions_grch37.json lookup file.
    """
    print(f"\n>>> Loading genome from {genome_path}")

    genome_by_rsid = {}
    genome_by_position = {}
    skipped = 0
    resolved = 0

    # Load position -> rsID map for WGS data
    pos_to_rsid = _load_position_to_rsid_map()

    with open(genome_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                rsid, chrom, pos, genotype = parts[0], parts[1], parts[2], parts[3]
                if genotype == '--':
                    continue
                # Validate genotype: must be 1-2 valid bases
                if not (1 <= len(genotype) <= 2 and all(b in _VALID_BASES for b in genotype)):
                    skipped += 1
                    continue

                # Resolve positional IDs to rsIDs if possible
                pos_key = f"{chrom}:{pos}"
                if not rsid.startswith("rs"):
                    mapped = pos_to_rsid.get(pos_key)
                    if mapped:
                        rsid = mapped
                        resolved += 1

                genome_by_rsid[rsid] = {
                    'chromosome': chrom,
                    'position': pos,
                    'genotype': genotype
                }
                genome_by_position[pos_key] = {
                    'rsid': rsid,
                    'genotype': genotype
                }

    print(f"    Loaded {len(genome_by_rsid):,} SNPs")
    if resolved:
        print(f"    Resolved {resolved:,} positional IDs to rsIDs via lookup")
    if skipped:
        print(f"    Skipped {skipped:,} entries with invalid genotypes")
    return genome_by_rsid, genome_by_position


def load_pharmgkb(data_dir: Path = None) -> dict:
    """Load PharmGKB drug-gene annotations."""
    if data_dir is None:
        data_dir = DATA_DIR

    annotations_path = data_dir / "clinical_annotations.tsv"
    alleles_path = data_dir / "clinical_ann_alleles.tsv"

    missing = []
    if not annotations_path.exists():
        missing.append(str(annotations_path))
    if not alleles_path.exists():
        missing.append(str(alleles_path))
    if missing:
        print(f"    PharmGKB files not found, skipping drug interactions: {', '.join(missing)}")
        return {}

    print("\n>>> Loading PharmGKB data")

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
