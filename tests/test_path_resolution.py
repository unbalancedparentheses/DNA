"""Tests to verify project structure is correct after restructure."""

from pathlib import Path


class TestProjectStructure:
    def test_package_dir_exists(self, project_root):
        assert (project_root / "genetic_health").is_dir()

    def test_package_init_exists(self, project_root):
        assert (project_root / "genetic_health" / "__init__.py").is_file()

    def test_main_entry_point_exists(self, project_root):
        assert (project_root / "genetic_health" / "pipeline.py").is_file()

    def test_snp_database_exists(self, project_root):
        assert (project_root / "genetic_health" / "snp_database.py").is_file()

    def test_clinical_context_exists(self, project_root):
        assert (project_root / "genetic_health" / "clinical_context.py").is_file()

    def test_reports_subpackage_exists(self, project_root):
        assert (project_root / "genetic_health" / "reports" / "__init__.py").is_file()

    def test_wgs_pipeline_exists(self, project_root):
        assert (project_root / "genetic_health" / "wgs_pipeline.py").is_file()

    def test_setup_script_exists(self, project_root):
        assert (project_root / "scripts" / "setup_reference.sh").is_file()

    def test_data_dir_exists(self, project_root):
        assert (project_root / "data").is_dir()

    def test_reports_dir_exists(self, project_root):
        reports = project_root / "reports"
        assert reports.exists() or reports.parent.exists()

    def test_claude_md_at_root(self, project_root):
        assert (project_root / "CLAUDE.md").is_file()

    def test_gitignore_at_root(self, project_root):
        assert (project_root / ".gitignore").is_file()

    def test_no_stale_pycache_in_nested_dir(self, project_root):
        assert not (project_root / "scripts" / "scripts" / "__pycache__").exists()

    def test_flake_nix_at_root(self, project_root):
        assert (project_root / "flake.nix").is_file()

    def test_makefile_at_root(self, project_root):
        assert (project_root / "Makefile").is_file()

    def test_pyproject_at_root(self, project_root):
        assert (project_root / "pyproject.toml").is_file()

    def test_path_resolution_consistent(self, project_root):
        from genetic_health import config
        assert config.BASE_DIR == project_root
        assert config.DATA_DIR == project_root / "data"
        assert config.REPORTS_DIR == project_root / "reports"
