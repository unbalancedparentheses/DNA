"""Tests for disease risk analysis (ClinVar logic in run_full_analysis.py)."""

import csv
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent / "scripts"))
from run_full_analysis import classify_zygosity, load_clinvar_and_analyze, DATA_DIR


def _make_finding(**overrides):
    """Create a minimal disease finding dict."""
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


class TestZygosityClassification:
    def test_homozygous(self):
        f = _make_finding(is_homozygous=True, is_heterozygous=False)
        status, desc = classify_zygosity(f)
        assert status == "AFFECTED"
        assert "Homozygous" in desc

    def test_heterozygous_recessive(self):
        f = _make_finding(
            is_homozygous=False,
            is_heterozygous=True,
            inheritance="Autosomal recessive",
        )
        status, desc = classify_zygosity(f)
        assert status == "CARRIER"
        assert "recessive" in desc.lower()

    def test_heterozygous_dominant(self):
        f = _make_finding(
            is_homozygous=False,
            is_heterozygous=True,
            inheritance="Autosomal dominant",
        )
        status, desc = classify_zygosity(f)
        assert status == "AFFECTED"
        assert "dominant" in desc.lower()

    def test_heterozygous_unknown_inheritance(self):
        f = _make_finding(
            is_homozygous=False,
            is_heterozygous=True,
            inheritance="",
        )
        status, desc = classify_zygosity(f)
        assert status == "HETEROZYGOUS"

    def test_no_variant(self):
        f = _make_finding(is_homozygous=False, is_heterozygous=False)
        status, _ = classify_zygosity(f)
        assert status == "UNKNOWN"


class TestVariantDetection:
    """Test ClinVar variant detection logic using a synthetic ClinVar TSV."""

    def _write_clinvar(self, rows: list[dict]) -> Path:
        fieldnames = [
            "chrom", "pos", "ref", "alt", "clinical_significance",
            "review_status", "gold_stars", "all_traits", "symbol",
            "inheritance_modes", "hgvs_p", "hgvs_c",
            "molecular_consequence", "xrefs",
        ]
        f = tempfile.NamedTemporaryFile(
            mode="w", suffix=".tsv", delete=False, dir=str(DATA_DIR),
            prefix="test_clinvar_",
        )
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            full = {k: "" for k in fieldnames}
            full.update(row)
            writer.writerow(full)
        f.close()
        return Path(f.name)

    def _run_analysis(self, clinvar_rows, genome_positions):
        """Run load_clinvar_and_analyze with synthetic data."""
        import run_full_analysis as mod

        clinvar_path = self._write_clinvar(clinvar_rows)
        orig = mod.DATA_DIR
        # Point DATA_DIR at the temp file's directory and rename
        # Instead, monkey-patch the path directly
        tmp_data = clinvar_path.parent
        clinvar_target = tmp_data / "clinvar_alleles.tsv"
        if clinvar_path != clinvar_target:
            clinvar_path.rename(clinvar_target)

        mod.DATA_DIR = tmp_data
        try:
            findings, stats = mod.load_clinvar_and_analyze(genome_positions)
        finally:
            mod.DATA_DIR = orig
            if clinvar_target.exists():
                clinvar_target.unlink()

        return findings, stats

    def test_pathogenic_snp_detected(self):
        findings, stats = self._run_analysis(
            clinvar_rows=[{
                "chrom": "7", "pos": "1000", "ref": "A", "alt": "G",
                "clinical_significance": "Pathogenic",
                "gold_stars": "3", "all_traits": "Cystic fibrosis",
                "symbol": "CFTR",
            }],
            genome_positions={
                "7:1000": {"rsid": "rs123", "genotype": "AG"},
            },
        )
        assert len(findings["pathogenic"]) == 1
        assert findings["pathogenic"][0]["gene"] == "CFTR"

    def test_ref_only_not_flagged(self):
        findings, _ = self._run_analysis(
            clinvar_rows=[{
                "chrom": "7", "pos": "1000", "ref": "A", "alt": "G",
                "clinical_significance": "Pathogenic",
                "gold_stars": "2", "all_traits": "Test",
                "symbol": "TEST",
            }],
            genome_positions={
                "7:1000": {"rsid": "rs123", "genotype": "AA"},
            },
        )
        assert len(findings["pathogenic"]) == 0

    def test_indels_skipped(self):
        findings, _ = self._run_analysis(
            clinvar_rows=[{
                "chrom": "7", "pos": "1000", "ref": "AG", "alt": "A",
                "clinical_significance": "Pathogenic",
                "gold_stars": "3", "all_traits": "Test",
                "symbol": "TEST",
            }],
            genome_positions={
                "7:1000": {"rsid": "rs123", "genotype": "AA"},
            },
        )
        assert len(findings["pathogenic"]) == 0

    def test_risk_factor_classified(self):
        findings, _ = self._run_analysis(
            clinvar_rows=[{
                "chrom": "1", "pos": "500", "ref": "C", "alt": "T",
                "clinical_significance": "risk factor",
                "gold_stars": "1", "all_traits": "Hypertension",
                "symbol": "AGT",
            }],
            genome_positions={
                "1:500": {"rsid": "rs456", "genotype": "CT"},
            },
        )
        assert len(findings["risk_factor"]) == 1

    def test_drug_response_classified(self):
        findings, _ = self._run_analysis(
            clinvar_rows=[{
                "chrom": "10", "pos": "200", "ref": "G", "alt": "A",
                "clinical_significance": "drug response",
                "gold_stars": "2", "all_traits": "Warfarin response",
                "symbol": "VKORC1",
            }],
            genome_positions={
                "10:200": {"rsid": "rs789", "genotype": "GA"},
            },
        )
        assert len(findings["drug_response"]) == 1

    def test_protective_classified(self):
        findings, _ = self._run_analysis(
            clinvar_rows=[{
                "chrom": "3", "pos": "300", "ref": "T", "alt": "C",
                "clinical_significance": "protective",
                "gold_stars": "1", "all_traits": "HIV resistance",
                "symbol": "CCR5",
            }],
            genome_positions={
                "3:300": {"rsid": "rs321", "genotype": "CC"},
            },
        )
        assert len(findings["protective"]) == 1

    def test_homozygous_detection(self):
        findings, _ = self._run_analysis(
            clinvar_rows=[{
                "chrom": "7", "pos": "1000", "ref": "A", "alt": "G",
                "clinical_significance": "Pathogenic",
                "gold_stars": "3", "all_traits": "Test",
                "symbol": "TEST",
            }],
            genome_positions={
                "7:1000": {"rsid": "rs123", "genotype": "GG"},
            },
        )
        assert findings["pathogenic"][0]["is_homozygous"] is True
