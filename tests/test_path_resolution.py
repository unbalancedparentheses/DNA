"""Tests to verify project structure is correct after restructure."""

from pathlib import Path


class TestProjectStructure:
    def test_scripts_dir_exists(self, project_root):
        assert (project_root / "scripts").is_dir()

    def test_main_entry_point_exists(self, project_root):
        assert (project_root / "scripts" / "run_full_analysis.py").is_file()

    def test_snp_database_exists(self, project_root):
        assert (project_root / "scripts" / "comprehensive_snp_database.py").is_file()

    def test_report_generator_exists(self, project_root):
        assert (project_root / "scripts" / "generate_exhaustive_report.py").is_file()

    def test_report_generators_module_exists(self, project_root):
        assert (project_root / "scripts" / "report_generators.py").is_file()

    def test_wgs_pipeline_exists(self, project_root):
        assert (project_root / "scripts" / "wgs_pipeline.py").is_file()

    def test_data_dir_exists(self, project_root):
        assert (project_root / "data").is_dir()

    def test_reports_dir_exists(self, project_root):
        reports = project_root / "reports"
        assert reports.exists() or reports.parent.exists()

    def test_claude_md_at_root(self, project_root):
        assert (project_root / "CLAUDE.md").is_file()

    def test_gitignore_at_root(self, project_root):
        assert (project_root / ".gitignore").is_file()

    def test_no_nested_scripts_dir(self, project_root):
        assert not (project_root / "scripts" / "scripts").exists()

    def test_legacy_duplicates_removed(self, project_root):
        scripts = project_root / "scripts"
        assert not (scripts / "analyze_genome.py").exists()
        assert not (scripts / "full_health_analysis.py").exists()
        assert not (scripts / "disease_risk_analyzer.py").exists()

    def test_no_stale_pycache_in_nested_dir(self, project_root):
        assert not (project_root / "scripts" / "scripts" / "__pycache__").exists()

    def test_flake_nix_at_root(self, project_root):
        assert (project_root / "flake.nix").is_file()

    def test_makefile_at_root(self, project_root):
        assert (project_root / "Makefile").is_file()

    def test_pyproject_at_root(self, project_root):
        assert (project_root / "pyproject.toml").is_file()

    def test_path_resolution_consistent(self, project_root):
        import run_full_analysis as m
        assert m.SCRIPT_DIR == project_root / "scripts"
        assert m.BASE_DIR == project_root
        assert m.DATA_DIR == project_root / "data"
        assert m.REPORTS_DIR == project_root / "reports"
