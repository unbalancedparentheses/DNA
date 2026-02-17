"""Tests for data update module."""

import csv
import json
from pathlib import Path

from genetic_health.update_data import (
    _load_versions,
    _save_versions,
    _map_gold_stars,
    validate_pharmgkb,
    show_status,
    DATA_VERSIONS_FILE,
    REVIEW_STATUS_STARS,
)


class TestGoldStarsMapping:
    def test_practice_guideline(self):
        assert _map_gold_stars("practice guideline") == 4

    def test_expert_panel(self):
        assert _map_gold_stars("reviewed by expert panel") == 4

    def test_multiple_no_conflicts(self):
        assert _map_gold_stars("criteria provided, multiple submitters, no conflicts") == 3

    def test_conflicting(self):
        assert _map_gold_stars("criteria provided, conflicting interpretations") == 2

    def test_single_submitter(self):
        assert _map_gold_stars("criteria provided, single submitter") == 1

    def test_no_assertion(self):
        assert _map_gold_stars("no assertion criteria provided") == 0

    def test_unknown_returns_zero(self):
        assert _map_gold_stars("some unknown string") == 0

    def test_case_insensitive(self):
        assert _map_gold_stars("Practice Guideline") == 4

    def test_empty_string(self):
        assert _map_gold_stars("") == 0


class TestVersionMetadata:
    def test_load_nonexistent(self, tmp_path):
        result = _load_versions(tmp_path)
        assert result == {}

    def test_save_and_load(self, tmp_path):
        versions = {"clinvar": {"updated": "2024-01-01", "source": "test"}}
        _save_versions(tmp_path, versions)
        loaded = _load_versions(tmp_path)
        assert loaded == versions

    def test_file_created(self, tmp_path):
        _save_versions(tmp_path, {"test": True})
        assert (tmp_path / DATA_VERSIONS_FILE).exists()

    def test_load_valid_json(self, tmp_path):
        path = tmp_path / DATA_VERSIONS_FILE
        path.write_text('{"clinvar": {"updated": "2024-01-01"}}')
        result = _load_versions(tmp_path)
        assert "clinvar" in result


class TestValidatePharmgkb:
    def _write_tsv(self, path, headers, rows=None):
        """Helper to write a TSV file."""
        with open(path, "w", newline="") as f:
            writer = csv.writer(f, delimiter="\t")
            writer.writerow(headers)
            if rows:
                for row in rows:
                    writer.writerow(row)

    def test_missing_files(self, tmp_path):
        result = validate_pharmgkb(tmp_path)
        assert result is False

    def test_valid_files(self, tmp_path):
        self._write_tsv(
            tmp_path / "clinical_annotations.tsv",
            ["Clinical Annotation ID", "Variant/Haplotypes", "Gene", "Drug(s)",
             "Phenotype(s)", "Level of Evidence", "Phenotype Category"],
            [["ann1", "rs762551", "CYP1A2", "caffeine", "Metabolism", "1A", "PK"]],
        )
        self._write_tsv(
            tmp_path / "clinical_ann_alleles.tsv",
            ["Clinical Annotation ID", "Genotype/Allele", "Annotation Text"],
            [["ann1", "AA", "Fast metabolizer"]],
        )
        result = validate_pharmgkb(tmp_path)
        assert result is True
        # Check metadata was saved
        versions = _load_versions(tmp_path)
        assert "pharmgkb" in versions

    def test_missing_columns(self, tmp_path):
        # Write file with wrong columns
        self._write_tsv(
            tmp_path / "clinical_annotations.tsv",
            ["Wrong Column", "Another"],
        )
        self._write_tsv(
            tmp_path / "clinical_ann_alleles.tsv",
            ["Clinical Annotation ID", "Genotype/Allele", "Annotation Text"],
        )
        result = validate_pharmgkb(tmp_path)
        assert result is False


class TestShowStatus:
    def test_runs_without_error(self, tmp_path, capsys):
        show_status(tmp_path)
        captured = capsys.readouterr()
        assert "Data Status" in captured.out

    def test_shows_missing_files(self, tmp_path, capsys):
        show_status(tmp_path)
        captured = capsys.readouterr()
        assert "MISSING" in captured.out


class TestReviewStatusStarsMapping:
    def test_all_keys_have_valid_values(self):
        for key, stars in REVIEW_STATUS_STARS.items():
            assert 0 <= stars <= 4, f"Stars {stars} out of range for '{key}'"

    def test_known_entries_exist(self):
        assert "practice guideline" in REVIEW_STATUS_STARS
        assert "reviewed by expert panel" in REVIEW_STATUS_STARS
