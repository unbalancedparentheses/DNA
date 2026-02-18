"""Tests for disease risk analysis (ClinVar logic) and disease report generation."""

from genetic_health.analysis import load_clinvar_and_analyze
from genetic_health.reports import classify_zygosity
from genetic_health.reports.markdown_reports import generate_unified_report


class TestZygosityClassification:
    def test_homozygous(self, make_finding):
        f = make_finding(is_homozygous=True, is_heterozygous=False)
        status, desc = classify_zygosity(f)
        assert status == "AFFECTED"
        assert "Homozygous" in desc

    def test_heterozygous_recessive(self, make_finding):
        f = make_finding(inheritance="Autosomal recessive")
        status, desc = classify_zygosity(f)
        assert status == "CARRIER"
        assert "recessive" in desc.lower()

    def test_heterozygous_dominant(self, make_finding):
        f = make_finding(inheritance="Autosomal dominant")
        status, desc = classify_zygosity(f)
        assert status == "AFFECTED"
        assert "dominant" in desc.lower()

    def test_heterozygous_unknown_inheritance(self, make_finding):
        f = make_finding(inheritance="")
        status, desc = classify_zygosity(f)
        assert status == "HETEROZYGOUS"

    def test_no_variant(self, make_finding):
        f = make_finding(is_homozygous=False, is_heterozygous=False)
        status, _ = classify_zygosity(f)
        assert status == "UNKNOWN"


class TestVariantDetection:
    """Test ClinVar variant detection logic using a synthetic ClinVar TSV."""

    def _run_analysis(self, clinvar_file, clinvar_rows, genome_positions):
        clinvar_path = clinvar_file(clinvar_rows)
        findings, stats = load_clinvar_and_analyze(
            genome_positions, data_dir=clinvar_path.parent
        )
        return findings, stats

    def test_pathogenic_snp_detected(self, clinvar_file):
        findings, stats = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "7", "pos": "1000", "ref": "A", "alt": "G",
                "clinical_significance": "Pathogenic",
                "gold_stars": "3", "all_traits": "Cystic fibrosis",
                "symbol": "CFTR",
            }],
            genome_positions={"7:1000": {"rsid": "rs123", "genotype": "AG"}},
        )
        assert len(findings["pathogenic"]) == 1
        assert findings["pathogenic"][0]["gene"] == "CFTR"

    def test_ref_only_not_flagged(self, clinvar_file):
        findings, _ = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "7", "pos": "1000", "ref": "A", "alt": "G",
                "clinical_significance": "Pathogenic",
                "gold_stars": "2", "all_traits": "Test", "symbol": "TEST",
            }],
            genome_positions={"7:1000": {"rsid": "rs123", "genotype": "AA"}},
        )
        assert len(findings["pathogenic"]) == 0

    def test_indels_skipped(self, clinvar_file):
        findings, _ = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "7", "pos": "1000", "ref": "AG", "alt": "A",
                "clinical_significance": "Pathogenic",
                "gold_stars": "3", "all_traits": "Test", "symbol": "TEST",
            }],
            genome_positions={"7:1000": {"rsid": "rs123", "genotype": "AA"}},
        )
        assert len(findings["pathogenic"]) == 0

    def test_risk_factor_classified(self, clinvar_file):
        findings, _ = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "1", "pos": "500", "ref": "C", "alt": "T",
                "clinical_significance": "risk factor",
                "gold_stars": "1", "all_traits": "Hypertension", "symbol": "AGT",
            }],
            genome_positions={"1:500": {"rsid": "rs456", "genotype": "CT"}},
        )
        assert len(findings["risk_factor"]) == 1

    def test_drug_response_classified(self, clinvar_file):
        findings, _ = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "10", "pos": "200", "ref": "G", "alt": "A",
                "clinical_significance": "drug response",
                "gold_stars": "2", "all_traits": "Warfarin response", "symbol": "VKORC1",
            }],
            genome_positions={"10:200": {"rsid": "rs789", "genotype": "GA"}},
        )
        assert len(findings["drug_response"]) == 1

    def test_protective_classified(self, clinvar_file):
        findings, _ = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "3", "pos": "300", "ref": "T", "alt": "C",
                "clinical_significance": "protective",
                "gold_stars": "1", "all_traits": "HIV resistance", "symbol": "CCR5",
            }],
            genome_positions={"3:300": {"rsid": "rs321", "genotype": "CC"}},
        )
        assert len(findings["protective"]) == 1

    def test_homozygous_detection(self, clinvar_file):
        findings, _ = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "7", "pos": "1000", "ref": "A", "alt": "G",
                "clinical_significance": "Pathogenic",
                "gold_stars": "3", "all_traits": "Test", "symbol": "TEST",
            }],
            genome_positions={"7:1000": {"rsid": "rs123", "genotype": "GG"}},
        )
        assert findings["pathogenic"][0]["is_homozygous"] is True

    def test_underscore_clinical_significance(self, clinvar_file):
        """clinical_significance with underscores should match (e.g. 'likely_pathogenic')."""
        findings, stats = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "1", "pos": "100", "ref": "A", "alt": "G",
                "clinical_significance": "Likely_pathogenic",
                "gold_stars": "2", "all_traits": "Test", "symbol": "TEST",
            }],
            genome_positions={"1:100": {"rsid": "rs100", "genotype": "AG"}},
        )
        assert len(findings["likely_pathogenic"]) == 1
        assert stats["likely_pathogenic_matched"] == 1

    def test_missing_gold_stars(self, clinvar_file):
        """Empty gold_stars should default to 0, not crash."""
        findings, _ = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "7", "pos": "1000", "ref": "A", "alt": "G",
                "clinical_significance": "Pathogenic",
                "gold_stars": "", "all_traits": "Test", "symbol": "TEST",
            }],
            genome_positions={"7:1000": {"rsid": "rs123", "genotype": "AG"}},
        )
        assert findings["pathogenic"][0]["gold_stars"] == 0

    def test_missing_chrom_skipped(self, clinvar_file):
        """Rows with empty chrom should be skipped."""
        findings, stats = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "", "pos": "1000", "ref": "A", "alt": "G",
                "clinical_significance": "Pathogenic",
                "gold_stars": "2", "all_traits": "Test", "symbol": "TEST",
            }],
            genome_positions={"7:1000": {"rsid": "rs123", "genotype": "AG"}},
        )
        assert len(findings["pathogenic"]) == 0
        assert stats["matched"] == 0

    def test_missing_optional_fields_default_empty(self, clinvar_file):
        """Optional fields like inheritance_modes should default to empty string."""
        findings, _ = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "7", "pos": "1000", "ref": "A", "alt": "G",
                "clinical_significance": "Pathogenic",
                "gold_stars": "3", "all_traits": "Cystic fibrosis", "symbol": "CFTR",
            }],
            genome_positions={"7:1000": {"rsid": "rs123", "genotype": "AG"}},
        )
        f = findings["pathogenic"][0]
        assert f["inheritance"] == ""
        assert f["hgvs_p"] == ""
        assert f["molecular_consequence"] == ""

    def test_no_other_significant_category(self, clinvar_file):
        """'association' variants should not appear — other_significant was removed."""
        findings, _ = self._run_analysis(
            clinvar_file,
            clinvar_rows=[{
                "chrom": "1", "pos": "100", "ref": "A", "alt": "G",
                "clinical_significance": "association",
                "gold_stars": "1", "all_traits": "Test", "symbol": "TEST",
            }],
            genome_positions={"1:100": {"rsid": "rs100", "genotype": "AG"}},
        )
        assert "other_significant" not in findings


class TestDiseaseReportFormatting:
    _MINIMAL_HEALTH = {
        "findings": [],
        "pharmgkb_findings": [],
        "summary": {"total_snps": 0, "analyzed_snps": 0,
                     "high_impact": 0, "moderate_impact": 0, "low_impact": 0},
    }

    def test_traits_whitespace_trimmed(self, tmp_path, make_finding):
        """Traits with leading/trailing whitespace after semicolon split should be trimmed."""
        disease_findings = {
            "pathogenic": [
                make_finding(
                    is_homozygous=True,
                    traits="  Cystic fibrosis ; CFTR-related disorders ",
                    gene="CFTR",
                    gold_stars=3,
                ),
            ],
            "likely_pathogenic": [],
            "risk_factor": [],
            "drug_response": [],
            "protective": [],
        }
        stats = {"total_clinvar": 100, "pathogenic_matched": 1, "likely_pathogenic_matched": 0}
        output = tmp_path / "report.md"
        generate_unified_report(self._MINIMAL_HEALTH, disease_findings, stats, output)
        content = output.read_text()
        # The header should have trimmed trait, not "  Cystic fibrosis "
        assert "CFTR" in content and "Cystic fibrosis" in content
        assert "  Cystic fibrosis " not in content

    def test_unified_report_with_none_disease_findings(self, tmp_path):
        """Report should generate without crashing when ClinVar is missing."""
        health_results = {
            "findings": [
                {"gene": "CYP1A2", "category": "Caffeine", "genotype": "AA",
                 "status": "fast", "description": "Fast metabolizer",
                 "magnitude": 2, "note": "", "rsid": "rs762551"},
            ],
            "pharmgkb_findings": [],
            "summary": {"total_snps": 1, "analyzed_snps": 1,
                         "high_impact": 0, "moderate_impact": 1, "low_impact": 0},
        }
        output = tmp_path / "report.md"
        generate_unified_report(health_results, None, None, output)
        content = output.read_text()
        assert "Genetic Health Report" in content
        assert "ClinVar data not available" in content

    def test_findings_dict_keeps_highest_magnitude(self, tmp_path):
        """When a gene has multiple rsIDs, the highest magnitude finding should be used."""
        health_results = {
            "findings": [
                {"gene": "MTHFR", "category": "Methylation", "genotype": "CT",
                 "status": "reduced", "description": "Reduced activity",
                 "magnitude": 1, "note": "", "rsid": "rs1801131"},
                {"gene": "MTHFR", "category": "Methylation", "genotype": "CT",
                 "status": "reduced", "description": "C677T heterozygous",
                 "magnitude": 3, "note": "", "rsid": "rs1801133"},
            ],
            "pharmgkb_findings": [],
            "summary": {"total_snps": 2, "analyzed_snps": 2,
                         "high_impact": 1, "moderate_impact": 0, "low_impact": 1},
        }
        output = tmp_path / "report.md"
        generate_unified_report(health_results, None, None, output)
        content = output.read_text()
        # With magnitude 3, MTHFR should trigger supplement recommendations
        assert "Methylfolate" in content
