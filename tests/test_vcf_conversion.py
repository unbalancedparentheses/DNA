"""Tests for VCF-to-23andMe conversion (step_convert in wgs_pipeline.py)."""

import gzip
import json
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))

# We can't import step_convert directly because it uses module-level globals.
# Instead, we inline the core logic into a testable helper.
import wgs_pipeline


def _make_vcf_gz(lines: list[str]) -> Path:
    """Write VCF lines to a gzipped temporary file."""
    header = [
        "##fileformat=VCFv4.2",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    f = tempfile.NamedTemporaryFile(suffix=".vcf.gz", delete=False)
    with gzip.open(f.name, "wt") as gz:
        gz.write("\n".join(header + lines) + "\n")
    return Path(f.name)


def _run_convert(vcf_lines: list[str], rsid_lookup: dict | None = None) -> list[str]:
    """Run step_convert and return output lines (non-comment)."""
    vcf_path = _make_vcf_gz(vcf_lines)
    out_path = Path(tempfile.mktemp(suffix=".txt"))

    # Set up rsID lookup
    lookup_path = Path(tempfile.mktemp(suffix=".json"))
    if rsid_lookup:
        with open(lookup_path, "w") as f:
            json.dump(rsid_lookup, f)
    else:
        with open(lookup_path, "w") as f:
            json.dump({}, f)

    # Monkey-patch the module-level constant
    orig = wgs_pipeline.RSID_LOOKUP
    wgs_pipeline.RSID_LOOKUP = lookup_path
    try:
        wgs_pipeline.step_convert(vcf_path, out_path)
    finally:
        wgs_pipeline.RSID_LOOKUP = orig

    with open(out_path) as f:
        return [l.strip() for l in f if not l.startswith("#")]


class TestVCFConversion:
    def test_basic_snp(self):
        lines = _run_convert([
            "1\t12345\trs100\tA\tG\t30\tPASS\t.\tGT\t0/1",
        ])
        assert len(lines) == 1
        parts = lines[0].split("\t")
        assert parts[0] == "rs100"
        assert parts[1] == "1"
        assert parts[2] == "12345"
        assert parts[3] == "AG"

    def test_homozygous_alt(self):
        lines = _run_convert([
            "1\t100\trs200\tA\tG\t30\tPASS\t.\tGT\t1/1",
        ])
        parts = lines[0].split("\t")
        assert parts[3] == "GG"

    def test_homozygous_ref(self):
        lines = _run_convert([
            "1\t100\trs300\tA\tG\t30\tPASS\t.\tGT\t0/0",
        ])
        parts = lines[0].split("\t")
        assert parts[3] == "AA"

    def test_indels_filtered(self):
        lines = _run_convert([
            "1\t100\t.\tA\tAG\t30\tPASS\t.\tGT\t0/1",    # insertion
            "1\t200\t.\tAG\tA\t30\tPASS\t.\tGT\t0/1",     # deletion
            "1\t300\t.\tA\tG\t30\tPASS\t.\tGT\t0/1",      # SNP — kept
        ])
        assert len(lines) == 1
        assert "300" in lines[0]

    def test_multi_allele_indel_filtered(self):
        lines = _run_convert([
            "1\t100\t.\tA\tG,AGT\t30\tPASS\t.\tGT\t0/1",  # mixed, one alt is indel
        ])
        assert len(lines) == 0

    def test_phased_genotype(self):
        """Phased genotypes use '|' separator."""
        lines = _run_convert([
            "1\t100\trs400\tA\tG\t30\tPASS\t.\tGT\t0|1",
        ])
        parts = lines[0].split("\t")
        assert parts[3] == "AG"

    def test_missing_genotype_skipped(self):
        lines = _run_convert([
            "1\t100\trs500\tA\tG\t30\tPASS\t.\tGT\t./.",
        ])
        assert len(lines) == 0

    def test_rsid_from_lookup(self):
        """When VCF ID is '.', use rsID lookup by position."""
        lookup = {"rs999": {"chrom": "1", "pos": "100"}}
        lines = _run_convert(
            ["1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1"],
            rsid_lookup=lookup,
        )
        parts = lines[0].split("\t")
        assert parts[0] == "rs999"

    def test_positional_rsid_fallback(self):
        """When no rsID and no lookup match, use chr_pos format."""
        lines = _run_convert([
            "5\t42000\t.\tA\tG\t30\tPASS\t.\tGT\t0/1",
        ])
        parts = lines[0].split("\t")
        assert parts[0] == "chr5_42000"

    def test_format_with_extra_fields(self):
        """GT may not be the only FORMAT field."""
        lines = _run_convert([
            "1\t100\trs600\tA\tG\t30\tPASS\t.\tGT:DP:GQ\t0/1:30:99",
        ])
        parts = lines[0].split("\t")
        assert parts[3] == "AG"
