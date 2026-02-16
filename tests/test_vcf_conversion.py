"""Tests for VCF-to-23andMe conversion (step_convert in wgs_pipeline.py)."""

import wgs_pipeline


def _run_convert(vcf_file, vcf_lines, rsid_lookup_file, rsid_lookup=None):
    """Run step_convert and return output lines (non-comment)."""
    vcf_path = vcf_file(vcf_lines)
    out_path = vcf_path.parent / "output.txt"
    lookup_path = rsid_lookup_file(rsid_lookup or {})

    orig = wgs_pipeline.RSID_LOOKUP
    wgs_pipeline.RSID_LOOKUP = lookup_path
    try:
        wgs_pipeline.step_convert(vcf_path, out_path)
    finally:
        wgs_pipeline.RSID_LOOKUP = orig

    return [l.strip() for l in out_path.read_text().splitlines() if not l.startswith("#")]


class TestVCFConversion:
    def test_basic_snp(self, vcf_file, rsid_lookup_file):
        lines = _run_convert(vcf_file, ["1\t12345\trs100\tA\tG\t30\tPASS\t.\tGT\t0/1"], rsid_lookup_file)
        assert len(lines) == 1
        parts = lines[0].split("\t")
        assert parts[0] == "rs100"
        assert parts[1] == "1"
        assert parts[2] == "12345"
        assert parts[3] == "AG"

    def test_homozygous_alt(self, vcf_file, rsid_lookup_file):
        lines = _run_convert(vcf_file, ["1\t100\trs200\tA\tG\t30\tPASS\t.\tGT\t1/1"], rsid_lookup_file)
        assert lines[0].split("\t")[3] == "GG"

    def test_homozygous_ref(self, vcf_file, rsid_lookup_file):
        lines = _run_convert(vcf_file, ["1\t100\trs300\tA\tG\t30\tPASS\t.\tGT\t0/0"], rsid_lookup_file)
        assert lines[0].split("\t")[3] == "AA"

    def test_indels_filtered(self, vcf_file, rsid_lookup_file):
        lines = _run_convert(vcf_file, [
            "1\t100\t.\tA\tAG\t30\tPASS\t.\tGT\t0/1",
            "1\t200\t.\tAG\tA\t30\tPASS\t.\tGT\t0/1",
            "1\t300\t.\tA\tG\t30\tPASS\t.\tGT\t0/1",
        ], rsid_lookup_file)
        assert len(lines) == 1
        assert "300" in lines[0]

    def test_multi_allele_indel_filtered(self, vcf_file, rsid_lookup_file):
        lines = _run_convert(vcf_file, ["1\t100\t.\tA\tG,AGT\t30\tPASS\t.\tGT\t0/1"], rsid_lookup_file)
        assert len(lines) == 0

    def test_phased_genotype(self, vcf_file, rsid_lookup_file):
        lines = _run_convert(vcf_file, ["1\t100\trs400\tA\tG\t30\tPASS\t.\tGT\t0|1"], rsid_lookup_file)
        assert lines[0].split("\t")[3] == "AG"

    def test_missing_genotype_skipped(self, vcf_file, rsid_lookup_file):
        lines = _run_convert(vcf_file, ["1\t100\trs500\tA\tG\t30\tPASS\t.\tGT\t./."], rsid_lookup_file)
        assert len(lines) == 0

    def test_rsid_from_lookup(self, vcf_file, rsid_lookup_file):
        lookup = {"rs999": {"chrom": "1", "pos": "100"}}
        lines = _run_convert(vcf_file, ["1\t100\t.\tA\tG\t30\tPASS\t.\tGT\t0/1"], rsid_lookup_file, lookup)
        assert lines[0].split("\t")[0] == "rs999"

    def test_positional_rsid_fallback(self, vcf_file, rsid_lookup_file):
        lines = _run_convert(vcf_file, ["5\t42000\t.\tA\tG\t30\tPASS\t.\tGT\t0/1"], rsid_lookup_file)
        assert lines[0].split("\t")[0] == "chr5_42000"

    def test_format_with_extra_fields(self, vcf_file, rsid_lookup_file):
        lines = _run_convert(vcf_file, ["1\t100\trs600\tA\tG\t30\tPASS\t.\tGT:DP:GQ\t0/1:30:99"], rsid_lookup_file)
        assert lines[0].split("\t")[3] == "AG"
