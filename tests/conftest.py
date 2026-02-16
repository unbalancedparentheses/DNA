"""Shared test fixtures for the genetic health analysis test suite."""

import csv
import gzip
import json
import sys
import tempfile
from pathlib import Path

import pytest

# Add scripts/ to path once for all tests
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

PROJECT_ROOT = Path(__file__).parent.parent


@pytest.fixture
def project_root():
    return PROJECT_ROOT


# --- Genome fixtures ---

@pytest.fixture
def genome_file(tmp_path):
    """Create a factory for temporary genome files."""
    def _make(lines: list[str]) -> Path:
        p = tmp_path / "genome.txt"
        p.write_text("\n".join(lines) + "\n")
        return p
    return _make


@pytest.fixture
def sample_genome(genome_file):
    """A small but realistic genome file."""
    return genome_file([
        "# rsid\tchromosome\tposition\tgenotype",
        "rs762551\t15\t75041917\tAC",
        "rs4244285\t10\t96541616\tGG",
        "rs123\t1\t12345\tAG",
        "rs456\t2\t67890\tCC",
        "rs999\tMT\t100\tA",
        "rs000\t1\t200\t--",
    ])


# --- ClinVar fixtures ---

CLINVAR_FIELDS = [
    "chrom", "pos", "ref", "alt", "clinical_significance",
    "review_status", "gold_stars", "all_traits", "symbol",
    "inheritance_modes", "hgvs_p", "hgvs_c",
    "molecular_consequence", "xrefs",
]


@pytest.fixture
def clinvar_file(tmp_path):
    """Create a factory for temporary ClinVar TSV files."""
    def _make(rows: list[dict]) -> Path:
        p = tmp_path / "clinvar_alleles.tsv"
        with open(p, "w", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=CLINVAR_FIELDS, delimiter="\t")
            writer.writeheader()
            for row in rows:
                full = {k: "" for k in CLINVAR_FIELDS}
                full.update(row)
                writer.writerow(full)
        return p
    return _make


# --- VCF fixtures ---

VCF_HEADER = [
    "##fileformat=VCFv4.2",
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
]


@pytest.fixture
def vcf_file(tmp_path):
    """Create a factory for gzipped VCF files."""
    def _make(lines: list[str]) -> Path:
        p = tmp_path / "variants.vcf.gz"
        with gzip.open(str(p), "wt") as gz:
            gz.write("\n".join(VCF_HEADER + lines) + "\n")
        return p
    return _make


@pytest.fixture
def rsid_lookup_file(tmp_path):
    """Create a factory for rsID lookup JSON files."""
    def _make(lookup: dict) -> Path:
        p = tmp_path / "rsid_positions_grch37.json"
        p.write_text(json.dumps(lookup))
        return p
    return _make


# --- Disease finding fixtures ---

@pytest.fixture
def make_finding():
    """Factory for disease finding dicts."""
    def _make(**overrides):
        base = {
            "chromosome": "1",
            "position": "100",
            "rsid": "rs000",
            "gene": "TEST",
            "ref": "A",
            "alt": "G",
            "user_genotype": "AG",
            "is_homozygous": False,
            "is_heterozygous": True,
            "clinical_significance": "Pathogenic",
            "review_status": "criteria provided",
            "gold_stars": 2,
            "traits": "Test condition",
            "inheritance": "",
            "hgvs_p": "",
            "hgvs_c": "",
            "molecular_consequence": "",
            "xrefs": "",
        }
        base.update(overrides)
        return base
    return _make
