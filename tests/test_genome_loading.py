"""Tests for 23andMe genome file parsing (load_genome in run_full_analysis.py)."""

import sys
import tempfile
from pathlib import Path

# Add scripts/ to path so we can import
sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from run_full_analysis import load_genome


def _write_genome(lines: list[str]) -> Path:
    """Write lines to a temporary genome file and return its path."""
    f = tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False)
    f.write("\n".join(lines) + "\n")
    f.close()
    return Path(f.name)


class TestGenomeLoading:
    def test_basic_parsing(self):
        path = _write_genome([
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

    def test_position_index(self):
        path = _write_genome([
            "# comment",
            "rs123\t1\t12345\tAG",
            "rs456\tX\t99999\tTT",
        ])
        _, by_pos = load_genome(path)
        assert "1:12345" in by_pos
        assert by_pos["1:12345"]["rsid"] == "rs123"
        assert by_pos["1:12345"]["genotype"] == "AG"
        assert "X:99999" in by_pos

    def test_no_call_skipped(self):
        path = _write_genome([
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

    def test_comment_lines_ignored(self):
        path = _write_genome([
            "# This is a header",
            "# Another comment",
            "rs100\t5\t500\tAT",
        ])
        by_rsid, _ = load_genome(path)
        assert len(by_rsid) == 1
        assert "rs100" in by_rsid

    def test_empty_file(self):
        path = _write_genome(["# just a comment"])
        by_rsid, by_pos = load_genome(path)
        assert len(by_rsid) == 0
        assert len(by_pos) == 0

    def test_short_lines_skipped(self):
        path = _write_genome([
            "# header",
            "rs100\t1",
            "rs200\t2\t300\tAG",
        ])
        by_rsid, _ = load_genome(path)
        assert len(by_rsid) == 1
        assert "rs200" in by_rsid

    def test_single_allele_genotype(self):
        """MT and Y chromosomes may have single-allele genotypes."""
        path = _write_genome([
            "# header",
            "rs999\tMT\t100\tA",
        ])
        by_rsid, _ = load_genome(path)
        assert by_rsid["rs999"]["genotype"] == "A"
