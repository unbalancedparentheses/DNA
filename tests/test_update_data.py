"""Tests for data update module."""

import csv
import gzip
import json
from pathlib import Path
from unittest.mock import patch

from genetic_health.update_data import (
    _load_versions,
    _save_versions,
    _map_gold_stars,
    update_clinvar,
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


class TestUpdateClinvar:
    """Test ClinVar processing with a mock gzip file (mocks the download)."""

    VARIANT_SUMMARY_HEADER = [
        "#AlleleID", "Type", "Name", "GeneID", "GeneSymbol", "HGNC_ID",
        "ClinicalSignificance", "ClinSigSimple", "LastEvaluated",
        "RS# (dbSNP)", "nsv/esv (dbVar)", "RCVaccession", "PhenotypeIDS",
        "PhenotypeList", "Origin", "OriginSimple", "Assembly",
        "ChromosomeAccession", "Chromosome", "Start", "Stop",
        "ReferenceAllele", "AlternateAllele", "Cytogenetic",
        "ReviewStatus", "NumberSubmitters", "Guidelines",
        "TestedInGTR", "OtherIDs", "SubmitterCategories",
        "VariationID", "PositionVCF", "ReferenceAlleleVCF",
        "AlternateAlleleVCF", "RS#(dbSNP2)", "ProteinChange",
        "MolecularConsequence",
    ]

    def _make_variant_summary_gz(self, tmp_path, rows):
        """Create a mock variant_summary.txt.gz in a subfolder to avoid same-file conflict."""
        src_dir = tmp_path / "_src"
        src_dir.mkdir(exist_ok=True)
        gz_path = src_dir / "variant_summary.txt.gz"
        with gzip.open(str(gz_path), "wt", encoding="utf-8") as gz:
            gz.write("\t".join(self.VARIANT_SUMMARY_HEADER) + "\n")
            for row in rows:
                full = {h: "" for h in self.VARIANT_SUMMARY_HEADER}
                full.update(row)
                gz.write("\t".join(full[h] for h in self.VARIANT_SUMMARY_HEADER) + "\n")
        return gz_path

    def _mock_download(self, url, dest):
        """Mock that copies the pre-built gz file to the expected dest."""
        import shutil
        shutil.copy2(str(self._gz_path), str(dest))

    def test_processes_grch37_snps(self, tmp_path):
        """GRCh37 rows should appear in output; GRCh38 rows should be filtered."""
        rows = [
            {
                "GeneSymbol": "CFTR",
                "ClinicalSignificance": "Pathogenic",
                "Assembly": "GRCh37",
                "Chromosome": "7",
                "Start": "117199646",
                "ReferenceAlleleVCF": "G",
                "AlternateAlleleVCF": "A",
                "ReviewStatus": "criteria provided, multiple submitters, no conflicts",
                "PhenotypeList": "Cystic fibrosis",
                "OriginSimple": "germline",
                "MolecularConsequence": "missense",
            },
            {
                "GeneSymbol": "BRCA1",
                "ClinicalSignificance": "Pathogenic",
                "Assembly": "GRCh38",  # Should be filtered out
                "Chromosome": "17",
                "Start": "43044295",
                "ReferenceAlleleVCF": "T",
                "AlternateAlleleVCF": "C",
                "ReviewStatus": "reviewed by expert panel",
                "PhenotypeList": "Breast cancer",
            },
        ]
        self._gz_path = self._make_variant_summary_gz(tmp_path, rows)

        with patch("genetic_health.update_data.urllib.request.urlretrieve", self._mock_download):
            update_clinvar(tmp_path)

        output = tmp_path / "clinvar_alleles.tsv"
        assert output.exists()

        with open(output) as f:
            reader = list(csv.DictReader(f, delimiter="\t"))

        # Only GRCh37 row should be present
        assert len(reader) == 1
        assert reader[0]["symbol"] == "CFTR"
        assert reader[0]["chrom"] == "7"
        assert reader[0]["gold_stars"] == "3"

    def test_metadata_written(self, tmp_path):
        """data_versions.json should be created after processing."""
        self._gz_path = self._make_variant_summary_gz(tmp_path, [
            {
                "GeneSymbol": "TEST",
                "ClinicalSignificance": "Benign",
                "Assembly": "GRCh37",
                "Chromosome": "1",
                "Start": "100",
                "ReferenceAlleleVCF": "A",
                "AlternateAlleleVCF": "G",
                "ReviewStatus": "no assertion criteria provided",
            },
        ])

        with patch("genetic_health.update_data.urllib.request.urlretrieve", self._mock_download):
            update_clinvar(tmp_path)

        versions = _load_versions(tmp_path)
        assert "clinvar" in versions
        assert versions["clinvar"]["grch37_variants_written"] == 1

    def test_download_cleaned_up(self, tmp_path):
        """The downloaded gz file should be removed after processing."""
        self._gz_path = self._make_variant_summary_gz(tmp_path, [])

        with patch("genetic_health.update_data.urllib.request.urlretrieve", self._mock_download):
            update_clinvar(tmp_path)

        assert not (tmp_path / "variant_summary.txt.gz").exists()


class TestReviewStatusStarsMapping:
    def test_all_keys_have_valid_values(self):
        for key, stars in REVIEW_STATUS_STARS.items():
            assert 0 <= stars <= 4, f"Stars {stars} out of range for '{key}'"

    def test_known_entries_exist(self):
        assert "practice guideline" in REVIEW_STATUS_STARS
        assert "reviewed by expert panel" in REVIEW_STATUS_STARS
