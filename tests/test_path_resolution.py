"""Tests to verify project structure is correct after restructure."""

from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent


class TestProjectStructure:
    def test_scripts_dir_exists(self):
        assert (PROJECT_ROOT / "scripts").is_dir()

    def test_main_entry_point_exists(self):
        assert (PROJECT_ROOT / "scripts" / "run_full_analysis.py").is_file()

    def test_snp_database_exists(self):
        assert (PROJECT_ROOT / "scripts" / "comprehensive_snp_database.py").is_file()

    def test_report_generator_exists(self):
        assert (PROJECT_ROOT / "scripts" / "generate_exhaustive_report.py").is_file()

    def test_wgs_pipeline_exists(self):
        assert (PROJECT_ROOT / "scripts" / "wgs_pipeline.py").is_file()

    def test_data_dir_exists(self):
        assert (PROJECT_ROOT / "data").is_dir()

    def test_reports_dir_exists(self):
        # reports/ may not exist yet if no analysis has run, but the parent should
        # allow creating it
        reports = PROJECT_ROOT / "reports"
        assert reports.exists() or reports.parent.exists()

    def test_claude_md_at_root(self):
        assert (PROJECT_ROOT / "CLAUDE.md").is_file()

    def test_gitignore_at_root(self):
        assert (PROJECT_ROOT / ".gitignore").is_file()

    def test_no_nested_scripts_dir(self):
        """The old scripts/scripts/ nesting should be gone."""
        assert not (PROJECT_ROOT / "scripts" / "scripts").exists()

    def test_legacy_duplicates_removed(self):
        scripts = PROJECT_ROOT / "scripts"
        assert not (scripts / "analyze_genome.py").exists()
        assert not (scripts / "full_health_analysis.py").exists()
        assert not (scripts / "disease_risk_analyzer.py").exists()

    def test_no_stale_pycache_in_nested_dir(self):
        """The old scripts/scripts/__pycache__ should not exist."""
        assert not (PROJECT_ROOT / "scripts" / "scripts" / "__pycache__").exists()

    def test_flake_nix_at_root(self):
        assert (PROJECT_ROOT / "flake.nix").is_file()

    def test_makefile_at_root(self):
        assert (PROJECT_ROOT / "Makefile").is_file()

    def test_path_resolution_consistent(self):
        """Verify SCRIPT_DIR.parent resolves to project root."""
        import sys
        sys.path.insert(0, str(PROJECT_ROOT / "scripts"))
        import run_full_analysis as mod

        script_dir = mod.SCRIPT_DIR
        base_dir = mod.BASE_DIR

        assert script_dir == PROJECT_ROOT / "scripts"
        assert base_dir == PROJECT_ROOT
        assert mod.DATA_DIR == PROJECT_ROOT / "data"
        assert mod.REPORTS_DIR == PROJECT_ROOT / "reports"
