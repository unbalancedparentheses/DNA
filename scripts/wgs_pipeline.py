#!/usr/bin/env python3
"""
WGS Pipeline: FASTQ -> Genetic Health Reports

Processes whole genome sequencing data through:
  1. Quality control (fastp)
  2. Alignment to GRCh37 (minimap2)
  3. BAM sorting and indexing (samtools)
  4. Variant calling (bcftools)
  5. rsID annotation (Ensembl GRCh37 API)
  6. Conversion to 23andMe format
  7. Full genetic health analysis

Usage:
    python scripts/wgs_pipeline.py /path/to/reads.fastq
    python scripts/wgs_pipeline.py /path/to/reads.fastq --name "Subject Name"
    python scripts/wgs_pipeline.py /path/to/reads.fastq --keep-intermediates
"""

import subprocess
import sys
import os
import gzip
import json
import shutil
import argparse
from pathlib import Path
from datetime import datetime

# Directory layout
SCRIPT_DIR = Path(__file__).parent.resolve()
PROJECT_DIR = SCRIPT_DIR.parent.resolve()
DATA_DIR = PROJECT_DIR / "data"
REPORTS_DIR = PROJECT_DIR / "reports"
REF_DIR = PROJECT_DIR / "reference"
WORK_DIR = PROJECT_DIR / "wgs_work"

# Reference files
REF_GENOME = REF_DIR / "human_g1k_v37.fasta"
REF_MMI = REF_DIR / "human_g1k_v37.mmi"
RSID_LOOKUP = DATA_DIR / "rsid_positions_grch37.json"


def log(msg):
    print(f"[{datetime.now().strftime('%H:%M:%S')}] {msg}", flush=True)


def run(cmd, desc=None):
    """Run a shell command, exit on failure."""
    if desc:
        log(desc)
    result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"Command failed: {cmd}")
        print(f"STDERR: {result.stderr[:3000]}")
        sys.exit(1)
    return result


def check_prereqs():
    """Verify tools and reference genome are available."""
    required = ['minimap2', 'samtools', 'bcftools']
    missing = [t for t in required if not shutil.which(t)]
    if missing:
        print(f"Missing tools: {', '.join(missing)}")
        print("Enter the Nix shell first: nix develop")
        sys.exit(1)

    if not shutil.which('fastp'):
        log("Warning: fastp not found. QC step will be skipped.")
        log("  Install via: brew install fastp")

    if not REF_GENOME.exists():
        print(f"Reference genome not found: {REF_GENOME}")
        print("Run setup first: bash scripts/setup_reference.sh")
        sys.exit(1)


def step_qc(fastq_in, threads):
    """Quality control and adapter trimming with fastp."""
    out = WORK_DIR / "trimmed.fastq.gz"
    json_out = WORK_DIR / "fastp.json"

    run(
        f"fastp -i {fastq_in} -o {out} "
        f"-j {json_out} -h {WORK_DIR}/fastp.html "
        f"-w {min(threads, 16)} --length_required 50",
        "Step 1/6: Quality control (fastp)"
    )

    with open(json_out) as f:
        qc = json.load(f)
    before = qc['summary']['before_filtering']['total_reads']
    after = qc['summary']['after_filtering']['total_reads']
    q30 = qc['summary']['after_filtering']['q30_rate']
    log(f"  Reads: {before:,} -> {after:,} ({after/before*100:.1f}% passed, Q30: {q30*100:.1f}%)")

    return out


def step_align(fastq_in, threads):
    """Align to GRCh37 with minimap2, pipe into samtools sort."""
    bam = WORK_DIR / "aligned.bam"

    ref = str(REF_MMI) if REF_MMI.exists() else str(REF_GENOME)
    align_threads = max(1, threads - 2)
    sort_threads = min(4, threads)

    run(
        f"minimap2 -a -x sr -t {align_threads} --secondary=no {ref} {fastq_in} "
        f"| samtools sort -@ {sort_threads} -m 2G -o {bam}",
        "Step 2/6: Alignment + sorting (minimap2 -> samtools sort)"
    )

    run(f"samtools index -@ {sort_threads} {bam}", "  Indexing BAM")

    result = run(f"samtools flagstat {bam}")
    for line in result.stdout.strip().split('\n')[:5]:
        log(f"  {line.strip()}")

    return bam


def step_call(bam_in, threads):
    """Call variants with bcftools mpileup + call."""
    vcf = WORK_DIR / "variants.vcf.gz"

    run(
        f"bcftools mpileup -f {REF_GENOME} -q 20 -Q 20 "
        f"--threads {max(1, threads//2)} -a FORMAT/DP {bam_in} "
        f"| bcftools call -m -v --ploidy GRCh37 "
        f"--threads {max(1, threads//4)} -Oz -o {vcf}",
        "Step 3/6: Variant calling (bcftools mpileup + call)"
    )

    run(f"bcftools index {vcf}")

    result = run(f"bcftools stats {vcf}")
    for line in result.stdout.split('\n'):
        if line.startswith('SN') and 'number of records' in line:
            log(f"  Variants called: {line.strip().split(chr(9))[-1]}")
        if line.startswith('SN') and 'number of SNPs' in line:
            log(f"  SNPs: {line.strip().split(chr(9))[-1]}")

    return vcf


def step_rsid_lookup():
    """Build rsID -> GRCh37 position lookup by querying Ensembl API."""
    log("Step 4/6: Building rsID position lookup")

    if RSID_LOOKUP.exists():
        with open(RSID_LOOKUP) as f:
            lookup = json.load(f)
        log(f"  Lookup already exists ({len(lookup)} rsIDs)")
        return

    sys.path.insert(0, str(SCRIPT_DIR))
    from comprehensive_snp_database import COMPREHENSIVE_SNPS

    rsids = list(COMPREHENSIVE_SNPS.keys())
    log(f"  Querying Ensembl GRCh37 API for {len(rsids)} rsIDs...")

    import urllib.request

    lookup = {}
    batch_size = 50

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

        try:
            with urllib.request.urlopen(req, timeout=30) as resp:
                results = json.loads(resp.read())

            for rsid, info in results.items():
                for m in info.get('mappings', []):
                    if m.get('seq_region_name', '').replace('chr', '') in \
                       [str(c) for c in range(1, 23)] + ['X', 'Y', 'MT']:
                        lookup[rsid] = {
                            'chrom': m['seq_region_name'],
                            'pos': str(m['start']),
                        }
                        break
        except Exception as e:
            log(f"  Warning: API batch {i//batch_size + 1} failed: {e}")

    with open(RSID_LOOKUP, 'w') as f:
        json.dump(lookup, f, indent=2)

    log(f"  Saved {len(lookup)} rsID positions")


def step_convert(vcf_in, genome_out):
    """Convert VCF to 23andMe tab-separated format."""
    log("Step 5/6: Converting VCF to 23andMe format")

    # Load rsID position lookup (position -> rsid)
    pos_to_rsid = {}
    if RSID_LOOKUP.exists():
        with open(RSID_LOOKUP) as f:
            lookup = json.load(f)
        for rsid, info in lookup.items():
            pos_to_rsid[f"{info['chrom']}:{info['pos']}"] = rsid

    count = 0
    with gzip.open(str(vcf_in), 'rt') as fin, open(str(genome_out), 'w') as fout:
        fout.write("# rsid\tchromosome\tposition\tgenotype\n")

        for line in fin:
            if line.startswith('#'):
                continue

            fields = line.split('\t', 10)
            if len(fields) < 10:
                continue

            chrom, pos, vid, ref, alt = fields[0], fields[1], fields[2], fields[3], fields[4]

            # Only single-nucleotide variants
            alt_alleles = alt.split(',')
            if len(ref) != 1 or any(len(a) != 1 for a in alt_alleles):
                continue

            # Parse genotype (first field in FORMAT)
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

            # Assign rsID: from VCF ID field, from lookup, or position-based
            pos_key = f"{chrom}:{pos}"
            if vid != '.' and vid.startswith('rs'):
                rsid = vid
            elif pos_key in pos_to_rsid:
                rsid = pos_to_rsid[pos_key]
            else:
                rsid = f"chr{chrom}_{pos}"

            fout.write(f"{rsid}\t{chrom}\t{pos}\t{genotype}\n")
            count += 1

    log(f"  Converted {count:,} SNPs to 23andMe format -> {genome_out}")
    return count


def step_analyze(genome_path, subject_name):
    """Run the existing genetic health analysis pipeline."""
    log("Step 6/6: Running genetic health analysis")

    sys.path.insert(0, str(SCRIPT_DIR))
    from run_full_analysis import run_full_analysis
    run_full_analysis(Path(genome_path), subject_name)


def main():
    parser = argparse.ArgumentParser(
        description="WGS Pipeline: FASTQ -> Genetic Health Reports",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python scripts/wgs_pipeline.py ../22374_1.fastq
  python scripts/wgs_pipeline.py ../22374_1.fastq --name "John Doe"
  python scripts/wgs_pipeline.py ../22374_1.fastq.gz --threads 8
  python scripts/wgs_pipeline.py ../22374_1.fastq --skip-analysis
        """
    )
    parser.add_argument('fastq', type=Path, help='Input FASTQ file (.fastq or .fastq.gz)')
    parser.add_argument('--name', '-n', type=str, default=None, help='Subject name for reports')
    parser.add_argument('--threads', '-t', type=int, default=12, help='CPU threads (default: 12)')
    parser.add_argument('--keep-intermediates', action='store_true',
                        help='Keep BAM and intermediate files')
    parser.add_argument('--skip-qc', action='store_true', help='Skip fastp QC step')
    parser.add_argument('--skip-analysis', action='store_true',
                        help='Stop after VCF conversion (no health analysis)')
    parser.add_argument('--output', '-o', type=Path, default=None,
                        help='Output genome.txt path (default: data/genome.txt)')

    args = parser.parse_args()

    fastq = args.fastq.resolve()
    if not fastq.exists():
        print(f"FASTQ file not found: {fastq}")
        sys.exit(1)

    size_gb = fastq.stat().st_size / 1e9

    print("=" * 60)
    print("WGS Pipeline: FASTQ -> Genetic Health Reports")
    print(f"Started: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Input:   {fastq.name} ({size_gb:.1f} GB)")
    print(f"Threads: {args.threads}")
    print("=" * 60)

    check_prereqs()
    WORK_DIR.mkdir(exist_ok=True)
    REPORTS_DIR.mkdir(exist_ok=True)

    genome_out = args.output or (DATA_DIR / "genome.txt")

    # Back up existing genome.txt
    if genome_out.exists():
        backup = genome_out.with_suffix('.txt.bak')
        shutil.copy2(genome_out, backup)
        log(f"Backed up existing genome.txt -> {backup.name}")

    # Step 1: QC
    if args.skip_qc or not shutil.which('fastp'):
        if not shutil.which('fastp'):
            log("Skipping QC (fastp not installed)")
        trimmed = fastq
    else:
        trimmed = step_qc(fastq, args.threads)

    # Step 2: Align + sort
    bam = step_align(trimmed, args.threads)

    # Step 3: Call variants
    vcf = step_call(bam, args.threads)

    # Step 4: Build rsID lookup
    step_rsid_lookup()

    # Step 5: Convert to 23andMe format
    variant_count = step_convert(vcf, genome_out)

    # Step 6: Run analysis
    if not args.skip_analysis:
        step_analyze(genome_out, args.name)

    # Cleanup
    if not args.keep_intermediates:
        log("Cleaning up intermediate files...")
        for f in WORK_DIR.glob("trimmed.*"):
            f.unlink()
        for f in WORK_DIR.glob("fastp.*"):
            f.unlink()

    print()
    print("=" * 60)
    print(f"Pipeline complete! {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    print(f"Variants: {variant_count:,} SNPs extracted")
    print(f"Genome:   {genome_out}")
    if not args.skip_analysis:
        print(f"Reports:  {REPORTS_DIR}/")
    print()
    print("Note: At ~1.7x coverage, many genomic positions lack reads.")
    print("Disease analysis (ClinVar) uses position matching and works well.")
    print("Lifestyle analysis may have fewer matches than 23andMe data.")
    print("=" * 60)


if __name__ == "__main__":
    main()
