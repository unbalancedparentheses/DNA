"""Tests for Markdown-to-HTML report conversion."""

from report_generators import _md_to_html, _write_html, _inline


class TestInlineFormatting:
    def test_bold(self):
        assert "<strong>bold</strong>" in _inline("**bold** text")

    def test_code(self):
        assert "<code>AG</code>" in _inline("genotype `AG` found")

    def test_italic(self):
        assert "<em>note</em>" in _inline("*note* here")

    def test_combined(self):
        result = _inline("**CFTR** `AG` *het*")
        assert "<strong>CFTR</strong>" in result
        assert "<code>AG</code>" in result
        assert "<em>het</em>" in result


class TestMdToHtml:
    def test_headers(self):
        html = _md_to_html("# Title\n\n## Section\n\n### Sub")
        assert "<h1>Title</h1>" in html
        assert "<h2>Section</h2>" in html
        assert "<h3>Sub</h3>" in html

    def test_horizontal_rule(self):
        html = _md_to_html("text\n\n---\n\nmore")
        assert "<hr>" in html

    def test_unordered_list(self):
        html = _md_to_html("- first\n- second\n- third")
        assert "<ul>" in html
        assert "<li>first</li>" in html
        assert "<li>third</li>" in html

    def test_table(self):
        md = "| Gene | Status |\n|------|--------|\n| CFTR | Carrier |\n| BRCA | Clear |"
        html = _md_to_html(md)
        assert "<table>" in html
        assert "<th>Gene</th>" in html
        assert "<td>Carrier</td>" in html
        assert "<td>Clear</td>" in html

    def test_paragraph(self):
        html = _md_to_html("This is a paragraph.")
        assert "<p>This is a paragraph.</p>" in html

    def test_inline_in_table(self):
        md = "| Gene | Geno |\n|------|------|\n| **CFTR** | `AG` |"
        html = _md_to_html(md)
        assert "<strong>CFTR</strong>" in html
        assert "<code>AG</code>" in html

    def test_empty_input(self):
        assert _md_to_html("") == ""

    def test_full_report_snippet(self):
        md = """# Disease Risk Report

**Generated:** 2026-01-01

---

## Executive Summary

| Category | Count |
|----------|-------|
| **Pathogenic** | 3 |
| **Carrier** | 5 |

---

## Findings

- **CFTR**: Cystic fibrosis carrier
- **HBB**: Sickle cell trait

*Generated using ClinVar*
"""
        html = _md_to_html(md)
        assert "<h1>Disease Risk Report</h1>" in html
        assert "<h2>Executive Summary</h2>" in html
        assert "<table>" in html
        assert "<ul>" in html
        assert "<strong>CFTR</strong>" in html
        assert "<hr>" in html


class TestWriteHtml:
    def test_writes_html_file(self, tmp_path):
        md_path = tmp_path / "report.md"
        md_content = "# Test Report\n\nSome content.\n"
        md_path.write_text(md_content)

        _write_html(md_content, md_path)

        html_path = tmp_path / "report.html"
        assert html_path.exists()
        html = html_path.read_text()
        assert "<!DOCTYPE html>" in html
        assert "<title>Test Report</title>" in html
        assert "<h1>Test Report</h1>" in html
        assert "Some content." in html

    def test_html_has_css(self, tmp_path):
        md_path = tmp_path / "report.md"
        md_content = "# Report\n"
        md_path.write_text(md_content)

        _write_html(md_content, md_path)

        html = (tmp_path / "report.html").read_text()
        assert "<style>" in html
        assert "font-family" in html

    def test_html_dark_mode_support(self, tmp_path):
        md_path = tmp_path / "report.md"
        md_path.write_text("# R\n")

        _write_html("# R\n", md_path)

        html = (tmp_path / "report.html").read_text()
        assert "prefers-color-scheme: dark" in html
