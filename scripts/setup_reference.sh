#!/bin/bash
# One-time setup: download GRCh37 reference genome and create indexes.
#
# Usage:
#   nix develop          # enter shell with tools
#   bash scripts/setup_reference.sh
#
# Downloads ~900MB compressed, ~3GB uncompressed.
# Creates minimap2 index (~3.5GB) and samtools faidx.
# Total disk: ~10GB in reference/

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
REF_DIR="$PROJECT_DIR/reference"
FASTA="$REF_DIR/human_g1k_v37.fasta"
MMI="$REF_DIR/human_g1k_v37.mmi"
URL="https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz"

echo "=== Reference Genome Setup ==="
echo "Project: $PROJECT_DIR"
echo "Reference dir: $REF_DIR"
echo ""

# Check tools
for tool in minimap2 samtools curl; do
    if ! command -v "$tool" &>/dev/null; then
        echo "ERROR: $tool not found. Enter Nix shell first: nix develop"
        exit 1
    fi
done

mkdir -p "$REF_DIR"

# Download reference genome
if [ -f "$FASTA" ]; then
    echo "Reference genome already exists: $FASTA"
else
    echo "Downloading GRCh37 reference genome (~900MB)..."
    echo "Source: $URL"
    curl -L --progress-bar -o "$FASTA.gz" "$URL"

    echo "Decompressing..."
    gunzip "$FASTA.gz"
    echo "Reference genome: $FASTA ($(du -h "$FASTA" | cut -f1))"
fi

# Create samtools faidx index
if [ -f "$FASTA.fai" ]; then
    echo "samtools faidx index exists"
else
    echo "Creating samtools faidx index..."
    samtools faidx "$FASTA"
fi

# Create minimap2 index (speeds up repeated alignment)
if [ -f "$MMI" ]; then
    echo "minimap2 index exists"
else
    echo "Creating minimap2 index (this takes ~5 minutes)..."
    minimap2 -x sr -d "$MMI" "$FASTA"
    echo "minimap2 index: $MMI ($(du -h "$MMI" | cut -f1))"
fi

echo ""
echo "=== Setup Complete ==="
echo "Reference: $FASTA"
echo "Indexes:   $FASTA.fai, $MMI"
echo ""
echo "Disk usage:"
du -sh "$REF_DIR"
echo ""
echo "Next: python scripts/wgs_pipeline.py /path/to/reads.fastq"
