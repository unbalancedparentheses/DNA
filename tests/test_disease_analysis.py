"""Tests for disease risk analysis (ClinVar logic)."""

from genetic_health.analysis import load_clinvar_and_analyze
from genetic_health.reports import classify_zygosity


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
