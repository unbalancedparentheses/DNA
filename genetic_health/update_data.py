"""Automated data updates for ClinVar and PharmGKB.

Downloads ClinVar variant_summary.txt.gz from NCBI FTP, filters to
GRCh37, and converts to the clinvar_alleles.tsv format used by the
analysis pipeline. Also validates manually-downloaded PharmGKB files.

Usage:
    python -m genetic_health.update_data clinvar      # Download + process ClinVar
    python -m genetic_health.update_data pharmgkb     # Validate PharmGKB files
    python -m genetic_health.update_data --status     # Show data versions
"""

import csv
import gzip
import json
import sys
import urllib.request
from datetime import datetime
from pathlib import Path

from .config import DATA_DIR

CLINVAR_URL = "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz"
DATA_VERSIONS_FILE = "data_versions.json"

# Review status -> gold stars mapping
REVIEW_STATUS_STARS = {
    "practice guideline": 4,
    "reviewed by expert panel": 4,
    "criteria provided, multiple submitters, no conflicts": 3,
    "criteria provided, multiple submitters, conflicting interpretations": 2,
    "criteria provided, conflicting interpretations": 2,
    "criteria provided, single submitter": 1,
    "no assertion for the individual variant": 0,
    "no assertion criteria provided": 0,
    "no assertion provided": 0,
}


def _load_versions(data_dir: Path) -> dict:
    """Load data version metadata."""
    versions_path = data_dir / DATA_VERSIONS_FILE
    if versions_path.exists():
        with open(versions_path) as f:
            return json.load(f)
    return {}


def _save_versions(data_dir: Path, versions: dict):
    """Save data version metadata."""
    versions_path = data_dir / DATA_VERSIONS_FILE
    with open(versions_path, "w") as f:
        json.dump(versions, f, indent=2)


def _map_gold_stars(review_status: str) -> int:
    """Map ClinVar review status string to gold star count."""
    status_lower = review_status.strip().lower()
    for key, stars in REVIEW_STATUS_STARS.items():
        if key in status_lower:
            return stars
    return 0


def update_clinvar(data_dir: Path = None):
    """Download and process ClinVar variant_summary.txt.gz.

    Steps:
        1. Download variant_summary.txt.gz from NCBI FTP
        2. Filter to GRCh37 assembly
        3. Filter to SNPs (single nucleotide changes)
        4. Map to clinvar_alleles.tsv format
        5. Update data_versions.json
    """
    if data_dir is None:
        data_dir = DATA_DIR
    data_dir.mkdir(parents=True, exist_ok=True)

    download_path = data_dir / "variant_summary.txt.gz"
    output_path = data_dir / "clinvar_alleles.tsv"

    # Step 1: Download
    print(f">>> Downloading ClinVar from {CLINVAR_URL}")
    print("    This may take a few minutes (~80MB)...")
    urllib.request.urlretrieve(CLINVAR_URL, str(download_path))
    print(f"    Downloaded to {download_path}")

    # Step 2-4: Process
    print(">>> Processing ClinVar data (filtering to GRCh37 SNPs)...")
    output_fields = [
        "chrom", "pos", "ref", "alt", "clinical_significance",
        "review_status", "gold_stars", "all_traits", "symbol",
        "inheritance_modes", "hgvs_p", "hgvs_c",
        "molecular_consequence", "xrefs",
    ]

    row_count = 0
    written = 0

    with gzip.open(str(download_path), "rt", encoding="utf-8") as gz_in, \
         open(output_path, "w", newline="") as tsv_out:

        reader = csv.DictReader(gz_in, delimiter="\t")
        writer = csv.DictWriter(tsv_out, fieldnames=output_fields, delimiter="\t")
        writer.writeheader()

        for row in reader:
            row_count += 1

            # Filter to GRCh37
            if row.get("Assembly") != "GRCh37":
                continue

            chrom = row.get("Chromosome", "")
            pos = row.get("Start", "")
            if not chrom or not pos:
                continue

            # We keep all variants (the analysis module filters indels)
            ref = row.get("ReferenceAlleleVCF", "")
            alt = row.get("AlternateAlleleVCF", "")
            if not ref or not alt:
                continue

            review_status = row.get("ReviewStatus", "")
            gold_stars = _map_gold_stars(review_status)

            out_row = {
                "chrom": chrom,
                "pos": pos,
                "ref": ref,
                "alt": alt,
                "clinical_significance": row.get("ClinicalSignificance", ""),
                "review_status": review_status,
                "gold_stars": str(gold_stars),
                "all_traits": row.get("PhenotypeList", ""),
                "symbol": row.get("GeneSymbol", ""),
                "inheritance_modes": row.get("OriginSimple", ""),
                "hgvs_p": row.get("ProteinChange", "") if "ProteinChange" in row else "",
                "hgvs_c": row.get("HGVS(c)", "") if "HGVS(c)" in row else "",
                "molecular_consequence": row.get("MolecularConsequence", ""),
                "xrefs": row.get("OtherIDs", ""),
            }
            writer.writerow(out_row)
            written += 1

    # Step 5: Clean up and update metadata
    download_path.unlink(missing_ok=True)

    versions = _load_versions(data_dir)
    versions["clinvar"] = {
        "updated": datetime.now().isoformat(),
        "source": CLINVAR_URL,
        "total_variants_processed": row_count,
        "grch37_variants_written": written,
        "file": str(output_path.name),
    }
    _save_versions(data_dir, versions)

    print(f"    Total ClinVar rows processed: {row_count:,}")
    print(f"    GRCh37 variants written: {written:,}")
    print(f"    Output: {output_path}")
    print(f"    Metadata updated: {data_dir / DATA_VERSIONS_FILE}")


def validate_pharmgkb(data_dir: Path = None):
    """Validate manually-downloaded PharmGKB files.

    PharmGKB requires account-based download, so we just validate
    the files exist and have expected columns.
    """
    if data_dir is None:
        data_dir = DATA_DIR

    annotations_path = data_dir / "clinical_annotations.tsv"
    alleles_path = data_dir / "clinical_ann_alleles.tsv"

    errors = []

    # Check annotations file
    if not annotations_path.exists():
        errors.append(f"Missing: {annotations_path}")
        print(f"    ERROR: {annotations_path} not found")
        print("    Download from: https://www.pharmgkb.org/downloads")
        print("    -> Clinical Annotation -> clinical_annotations.tsv")
    else:
        with open(annotations_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            headers = reader.fieldnames or []
            expected = ["Clinical Annotation ID", "Variant/Haplotypes", "Gene", "Drug(s)"]
            missing = [h for h in expected if h not in headers]
            if missing:
                errors.append(f"Missing columns in {annotations_path.name}: {missing}")
                print(f"    WARNING: Missing columns: {missing}")
            else:
                row_count = sum(1 for _ in reader)
                print(f"    {annotations_path.name}: OK ({row_count:,} rows)")

    # Check alleles file
    if not alleles_path.exists():
        errors.append(f"Missing: {alleles_path}")
        print(f"    ERROR: {alleles_path} not found")
        print("    Download from: https://www.pharmgkb.org/downloads")
        print("    -> Clinical Annotation -> clinical_ann_alleles.tsv")
    else:
        with open(alleles_path) as f:
            reader = csv.DictReader(f, delimiter="\t")
            headers = reader.fieldnames or []
            expected = ["Clinical Annotation ID", "Genotype/Allele", "Annotation Text"]
            missing = [h for h in expected if h not in headers]
            if missing:
                errors.append(f"Missing columns in {alleles_path.name}: {missing}")
                print(f"    WARNING: Missing columns: {missing}")
            else:
                row_count = sum(1 for _ in reader)
                print(f"    {alleles_path.name}: OK ({row_count:,} rows)")

    if not errors:
        versions = _load_versions(data_dir)
        versions["pharmgkb"] = {
            "validated": datetime.now().isoformat(),
            "files": [annotations_path.name, alleles_path.name],
        }
        _save_versions(data_dir, versions)
        print("    PharmGKB validation passed. Metadata updated.")
    else:
        print(f"\n    PharmGKB validation failed: {len(errors)} error(s)")

    return len(errors) == 0


def show_status(data_dir: Path = None):
    """Print data version information and file status."""
    if data_dir is None:
        data_dir = DATA_DIR

    print("=" * 60)
    print("Data Status")
    print("=" * 60)

    versions = _load_versions(data_dir)

    # ClinVar
    clinvar_path = data_dir / "clinvar_alleles.tsv"
    print(f"\nClinVar ({clinvar_path.name}):")
    if clinvar_path.exists():
        size_mb = clinvar_path.stat().st_size / (1024 * 1024)
        print(f"    File: EXISTS ({size_mb:.1f} MB)")
        if "clinvar" in versions:
            cv = versions["clinvar"]
            print(f"    Last updated: {cv.get('updated', 'unknown')}")
            print(f"    Source: {cv.get('source', 'unknown')}")
            print(f"    GRCh37 variants: {cv.get('grch37_variants_written', 'unknown'):,}" if isinstance(cv.get('grch37_variants_written'), int) else f"    GRCh37 variants: {cv.get('grch37_variants_written', 'unknown')}")
        else:
            print("    No version metadata (manually placed file)")
    else:
        print("    File: MISSING")
        print("    Run: python -m genetic_health.update_data clinvar")

    # PharmGKB
    print(f"\nPharmGKB:")
    for fname in ["clinical_annotations.tsv", "clinical_ann_alleles.tsv"]:
        fpath = data_dir / fname
        if fpath.exists():
            size_kb = fpath.stat().st_size / 1024
            print(f"    {fname}: EXISTS ({size_kb:.1f} KB)")
        else:
            print(f"    {fname}: MISSING")

    if "pharmgkb" in versions:
        pg = versions["pharmgkb"]
        print(f"    Last validated: {pg.get('validated', 'unknown')}")
    else:
        print("    Not validated yet")

    # Genome
    genome_path = data_dir / "genome.txt"
    print(f"\nGenome ({genome_path.name}):")
    if genome_path.exists():
        size_mb = genome_path.stat().st_size / (1024 * 1024)
        print(f"    File: EXISTS ({size_mb:.1f} MB)")
    else:
        print("    File: MISSING")

    print()


def main():
    """CLI entry point."""
    if len(sys.argv) < 2 or sys.argv[1] == "--status":
        show_status()
        return

    command = sys.argv[1]

    if command == "clinvar":
        update_clinvar()
    elif command == "pharmgkb":
        print(">>> Validating PharmGKB files")
        validate_pharmgkb()
    else:
        print(f"Unknown command: {command}")
        print("Usage: python -m genetic_health.update_data [clinvar|pharmgkb|--status]")
        sys.exit(1)


if __name__ == "__main__":
    main()
