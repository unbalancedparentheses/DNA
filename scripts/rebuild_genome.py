#!/usr/bin/env python3
"""Rebuild genome.txt with proper rsID annotations.

Collects all rsIDs from every analysis module, queries Ensembl GRCh37 for
positions, re-converts the VCF with proper rsID mapping, and extracts any
missing positions directly from the BAM.

Usage (run from project root inside nix shell):
    nix develop --command python scripts/rebuild_genome.py
"""

import json
import gzip
import subprocess
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

# Project paths
BASE_DIR = Path(__file__).parent.parent
DATA_DIR = BASE_DIR / "data"
WORK_DIR = BASE_DIR / "wgs_work"
RSID_LOOKUP = DATA_DIR / "rsid_positions_grch37.json"
GENOME_OUT = DATA_DIR / "genome.txt"
VCF_PATH = WORK_DIR / "variants.vcf.gz"
BAM_PATH = WORK_DIR / "aligned.bam"
REF_GENOME = BASE_DIR / "reference" / "human_g1k_v37.fasta"

sys.path.insert(0, str(BASE_DIR))


def log(msg):
    print(f"  {msg}", flush=True)


def collect_all_rsids():
    """Collect every rsID used across all analysis modules."""
    from genetic_health.snp_database import COMPREHENSIVE_SNPS
    from genetic_health.ancestry import ANCESTRY_MARKERS
    from genetic_health.prs import PRS_MODELS
    from genetic_health.mt_haplogroup import MT_HAPLOGROUP_TREE
    from genetic_health.star_alleles import STAR_ALLELE_DEFINITIONS

    all_rsids = set()

    all_rsids.update(COMPREHENSIVE_SNPS.keys())
    all_rsids.update(ANCESTRY_MARKERS.keys())

    for model in PRS_MODELS.values():
        for snp in model['snps']:
            all_rsids.add(snp['rsid'])

    # Blood type
    all_rsids.update(['rs505922', 'rs8176746', 'rs590787'])

    # MT haplogroup
    def _collect_mt(tree):
        rsids = set()
        for node in tree:
            if 'rsid' in node:
                rsids.add(node['rsid'])
            if 'children' in node:
                rsids.update(_collect_mt(node['children']))
        return rsids
    all_rsids.update(_collect_mt(MT_HAPLOGROUP_TREE))

    # Star alleles
    for gene, defn in STAR_ALLELE_DEFINITIONS.items():
        all_rsids.update(defn.get('snps', []))

    # APOE
    all_rsids.update(['rs429358', 'rs7412'])

    # Traits
    all_rsids.update(['rs12913832', 'rs1800407', 'rs1805007', 'rs1805008', 'rs17822931'])

    return sorted(all_rsids)


def query_ensembl(rsids):
    """Query Ensembl GRCh37 REST API for rsID positions in batches."""
    lookup = {}
    batch_size = 200
    valid_chroms = {str(c) for c in range(1, 23)} | {'X', 'Y', 'MT'}

    for i in range(0, len(rsids), batch_size):
        batch = rsids[i:i + batch_size]
        payload = json.dumps({"ids": batch}).encode()

        req = urllib.request.Request(
            "https://grch37.rest.ensembl.org/variation/homo_sapiens",
            data=payload,
            headers={
                "Content-Type": "application/json",
                "Accept": "application/json"
            }
        )

        retries = 3
        for attempt in range(retries):
            try:
                with urllib.request.urlopen(req, timeout=60) as resp:
                    results = json.loads(resp.read())

                for rsid, info in results.items():
                    if isinstance(info, dict) and 'mappings' in info:
                        for m in info['mappings']:
                            chrom = m.get('seq_region_name', '').replace('chr', '')
                            if chrom in valid_chroms:
                                lookup[rsid] = {
                                    'chrom': chrom,
                                    'pos': str(m['start']),
                                }
                                break
                break  # success
            except (urllib.error.URLError, json.JSONDecodeError, KeyError) as e:
                if attempt < retries - 1:
                    wait = 2 ** (attempt + 1)
                    log(f"Batch {i // batch_size + 1} attempt {attempt + 1} failed: {e}, retrying in {wait}s...")
                    time.sleep(wait)
                else:
                    log(f"Batch {i // batch_size + 1} failed after {retries} attempts: {e}")

        log(f"Batch {i // batch_size + 1}/{(len(rsids) + batch_size - 1) // batch_size}: "
            f"{len(lookup)} rsIDs resolved so far")

        # Ensembl rate limit: max 15 requests/sec
        time.sleep(0.2)

    return lookup


def convert_vcf(vcf_path, pos_to_rsid):
    """Re-convert VCF to 23andMe format using the complete rsID lookup."""
    entries = {}  # pos_key -> (rsid, chrom, pos, genotype)

    with gzip.open(str(vcf_path), 'rt') as fin:
        for line in fin:
            if line.startswith('#'):
                continue

            fields = line.split('\t', 10)
            if len(fields) < 10:
                continue

            chrom, pos, vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]

            # Only SNPs
            alt_alleles = alt.split(',')
            if len(ref) != 1 or any(len(a) != 1 for a in alt_alleles):
                continue

            # Parse genotype
            fmt_fields = fields[8].split(':')
            sample_fields = fields[9].split(':')
            try:
                gt_idx = fmt_fields.index('GT')
                gt = sample_fields[gt_idx]
            except (ValueError, IndexError):
                continue

            sep = '/' if '/' in gt else '|'
            gt_parts = gt.split(sep)
            if len(gt_parts) != 2 or '.' in gt_parts:
                continue

            alleles = [ref] + alt_alleles
            try:
                a1 = alleles[int(gt_parts[0])]
                a2 = alleles[int(gt_parts[1])]
            except (ValueError, IndexError):
                continue

            genotype = a1 + a2

            # Assign rsID
            pos_key = f"{chrom}:{pos}"
            if vid != '.' and vid.startswith('rs'):
                rsid = vid
            elif pos_key in pos_to_rsid:
                rsid = pos_to_rsid[pos_key]
            else:
                rsid = f"chr{chrom}_{pos}"

            entries[pos_key] = (rsid, chrom, pos, genotype)

    return entries


def extract_from_bam(bam_path, ref_genome, missing_positions):
    """Extract genotypes at specific positions from BAM using bcftools.

    missing_positions: list of (chrom, pos, rsid) not found in VCF.
    """
    if not missing_positions or not bam_path.exists():
        return {}

    # Write positions to a temporary regions file
    regions_file = WORK_DIR / "missing_regions.txt"
    with open(regions_file, 'w') as f:
        for chrom, pos, rsid in missing_positions:
            f.write(f"{chrom}\t{int(pos) - 1}\t{pos}\n")

    log(f"Extracting {len(missing_positions)} positions from BAM...")

    # Use bcftools mpileup + call for targeted positions
    cmd = (
        f"bcftools mpileup -f {ref_genome} -q 20 -Q 20 "
        f"-T {regions_file} -a FORMAT/DP {bam_path} "
        f"| bcftools call -m --ploidy GRCh37 -Ov"
    )

    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        log(f"bcftools failed: {result.stderr[:500]}")
        return {}

    # Build position→rsid lookup
    pos_to_rsid = {}
    for chrom, pos, rsid in missing_positions:
        pos_to_rsid[f"{chrom}:{pos}"] = rsid

    entries = {}
    for line in result.stdout.split('\n'):
        if line.startswith('#') or not line.strip():
            continue

        fields = line.split('\t', 10)
        if len(fields) < 10:
            continue

        chrom, pos, vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]

        # Parse genotype
        fmt_fields = fields[8].split(':')
        sample_fields = fields[9].strip().split(':')
        try:
            gt_idx = fmt_fields.index('GT')
            gt = sample_fields[gt_idx]
        except (ValueError, IndexError):
            continue

        sep = '/' if '/' in gt else '|'
        gt_parts = gt.split(sep)
        if len(gt_parts) != 2 or '.' in gt_parts:
            continue

        alleles = [ref] + (alt.split(',') if alt != '.' else [])
        try:
            a1 = alleles[int(gt_parts[0])]
            a2 = alleles[int(gt_parts[1])]
        except (ValueError, IndexError):
            continue

        # Only single-nucleotide
        if len(a1) != 1 or len(a2) != 1:
            continue

        genotype = a1 + a2
        pos_key = f"{chrom}:{pos}"
        rsid = pos_to_rsid.get(pos_key, f"chr{chrom}_{pos}")

        entries[pos_key] = (rsid, chrom, pos, genotype)

    return entries


def write_genome(entries, output_path):
    """Write genome.txt in 23andMe format, sorted by chromosome and position."""
    chrom_order = {str(i): i for i in range(1, 23)}
    chrom_order.update({'X': 23, 'Y': 24, 'MT': 25})

    sorted_entries = sorted(
        entries.values(),
        key=lambda x: (chrom_order.get(x[1], 99), int(x[2]))
    )

    with open(output_path, 'w') as f:
        f.write("# rsid\tchromosome\tposition\tgenotype\n")
        for rsid, chrom, pos, genotype in sorted_entries:
            f.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")

    return len(sorted_entries)


def main():
    print("=" * 60)
    print("Rebuild genome.txt with complete rsID annotations")
    print("=" * 60)

    # Step 1: Collect all rsIDs
    print("\n>>> Collecting rsIDs from all analysis modules")
    all_rsids = collect_all_rsids()
    log(f"{len(all_rsids)} unique rsIDs across all modules")

    # Step 2: Load existing lookup + query Ensembl for missing
    print("\n>>> Building rsID position lookup")
    existing = {}
    if RSID_LOOKUP.exists():
        with open(RSID_LOOKUP) as f:
            existing = json.load(f)
        log(f"Existing lookup: {len(existing)} rsIDs")

    missing = [r for r in all_rsids if r not in existing]
    if missing:
        log(f"Querying Ensembl GRCh37 for {len(missing)} missing rsIDs...")
        new_lookup = query_ensembl(missing)
        existing.update(new_lookup)
        log(f"Resolved {len(new_lookup)} new rsIDs")
    else:
        log("All rsIDs already in lookup")

    # Save updated lookup
    with open(RSID_LOOKUP, 'w') as f:
        json.dump(existing, f, indent=2)
    log(f"Saved {len(existing)} rsIDs to {RSID_LOOKUP.name}")

    # Build position→rsid reverse lookup
    pos_to_rsid = {}
    for rsid, info in existing.items():
        pos_to_rsid[f"{info['chrom']}:{info['pos']}"] = rsid

    # Step 3: Re-convert VCF with complete rsID mapping
    print("\n>>> Re-converting VCF with complete rsID mapping")
    if not VCF_PATH.exists():
        print(f"ERROR: VCF not found: {VCF_PATH}")
        sys.exit(1)

    entries = convert_vcf(VCF_PATH, pos_to_rsid)
    rsid_count = sum(1 for r, _, _, _ in entries.values() if r.startswith('rs'))
    log(f"VCF: {len(entries)} SNPs ({rsid_count} with rsIDs)")

    # Step 4: Find positions missing from VCF and extract from BAM
    print("\n>>> Checking for positions missing from VCF")
    vcf_positions = set(entries.keys())
    missing_from_vcf = []
    for rsid, info in existing.items():
        pos_key = f"{info['chrom']}:{info['pos']}"
        if pos_key not in vcf_positions:
            missing_from_vcf.append((info['chrom'], info['pos'], rsid))

    log(f"{len(missing_from_vcf)} positions not in VCF")

    if missing_from_vcf and BAM_PATH.exists() and REF_GENOME.exists():
        print("\n>>> Extracting missing positions from BAM")
        bam_entries = extract_from_bam(BAM_PATH, REF_GENOME, missing_from_vcf)
        log(f"Extracted {len(bam_entries)} genotypes from BAM")
        entries.update(bam_entries)
    elif missing_from_vcf:
        log("BAM or reference not available, skipping BAM extraction")

    # Step 5: Write genome.txt
    print("\n>>> Writing genome.txt")
    total = write_genome(entries, GENOME_OUT)
    rsid_count = sum(1 for r, _, _, _ in entries.values() if r.startswith('rs'))
    log(f"Wrote {total:,} SNPs ({rsid_count} with rsIDs, {total - rsid_count} position-based)")
    log(f"Output: {GENOME_OUT}")

    print("\n" + "=" * 60)
    print(f"Done! {total:,} variants in genome.txt ({rsid_count} annotated with rsIDs)")
    print("=" * 60)


if __name__ == "__main__":
    main()
