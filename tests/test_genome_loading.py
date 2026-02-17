"""Tests for 23andMe genome file parsing and PharmGKB loading."""

from genetic_health.loading import load_genome, load_pharmgkb


class TestGenomeLoading:
    def test_basic_parsing(self, genome_file):
        path = genome_file([
            "# rsid\tchromosome\tposition\tgenotype",
            "rs123\t1\t12345\tAG",
            "rs456\t2\t67890\tCC",
        ])
        by_rsid, by_pos = load_genome(path)
        assert "rs123" in by_rsid
        assert by_rsid["rs123"]["genotype"] == "AG"
        assert by_rsid["rs123"]["chromosome"] == "1"
        assert by_rsid["rs123"]["position"] == "12345"
        assert "rs456" in by_rsid
        assert len(by_rsid) == 2

    def test_position_index(self, genome_file):
        path = genome_file([
            "# comment",
            "rs123\t1\t12345\tAG",
            "rs456\tX\t99999\tTT",
        ])
        _, by_pos = load_genome(path)
        assert "1:12345" in by_pos
        assert by_pos["1:12345"]["rsid"] == "rs123"
        assert by_pos["1:12345"]["genotype"] == "AG"
        assert "X:99999" in by_pos

    def test_no_call_skipped(self, genome_file):
        path = genome_file([
            "# header",
            "rs111\t1\t100\tAA",
            "rs222\t1\t200\t--",
            "rs333\t1\t300\tGG",
        ])
        by_rsid, by_pos = load_genome(path)
        assert "rs111" in by_rsid
        assert "rs222" not in by_rsid
        assert "rs333" in by_rsid
        assert len(by_rsid) == 2

    def test_comment_lines_ignored(self, genome_file):
        path = genome_file([
            "# This is a header",
            "# Another comment",
            "rs100\t5\t500\tAT",
        ])
        by_rsid, _ = load_genome(path)
        assert len(by_rsid) == 1
        assert "rs100" in by_rsid

    def test_empty_file(self, genome_file):
        path = genome_file(["# just a comment"])
        by_rsid, by_pos = load_genome(path)
        assert len(by_rsid) == 0
        assert len(by_pos) == 0

    def test_short_lines_skipped(self, genome_file):
        path = genome_file([
            "# header",
            "rs100\t1",
            "rs200\t2\t300\tAG",
        ])
        by_rsid, _ = load_genome(path)
        assert len(by_rsid) == 1
        assert "rs200" in by_rsid

    def test_single_allele_genotype(self, genome_file):
        """MT and Y chromosomes may have single-allele genotypes."""
        path = genome_file([
            "# header",
            "rs999\tMT\t100\tA",
        ])
        by_rsid, _ = load_genome(path)
        assert by_rsid["rs999"]["genotype"] == "A"


class TestPharmGKBLoading:
    def test_missing_annotations_file(self, tmp_path, capsys):
        """Should report which file is missing."""
        # Create only the alleles file
        (tmp_path / "clinical_ann_alleles.tsv").write_text("")
        result = load_pharmgkb(data_dir=tmp_path)
        assert result == {}
        output = capsys.readouterr().out
        assert "clinical_annotations.tsv" in output

    def test_missing_alleles_file(self, tmp_path, capsys):
        """Should report which file is missing."""
        (tmp_path / "clinical_annotations.tsv").write_text("")
        result = load_pharmgkb(data_dir=tmp_path)
        assert result == {}
        output = capsys.readouterr().out
        assert "clinical_ann_alleles.tsv" in output

    def test_both_files_missing(self, tmp_path, capsys):
        """Should report both missing files."""
        result = load_pharmgkb(data_dir=tmp_path)
        assert result == {}
        output = capsys.readouterr().out
        assert "clinical_annotations.tsv" in output
        assert "clinical_ann_alleles.tsv" in output
