"""PDF export for HTML reports.

Attempts multiple backends in order:
1. Chrome/Chromium headless (most commonly available)
2. weasyprint (Python library, requires system deps)
3. Falls back to opening browser for manual print-to-PDF
"""

import subprocess
import sys
from pathlib import Path


# Chrome/Chromium paths by platform
_CHROME_PATHS = [
    "/Applications/Google Chrome.app/Contents/MacOS/Google Chrome",
    "/Applications/Chromium.app/Contents/MacOS/Chromium",
    "google-chrome",
    "google-chrome-stable",
    "chromium",
    "chromium-browser",
]


def _try_chrome_headless(html_path: Path, pdf_path: Path) -> bool:
    """Try to convert HTML to PDF using Chrome headless mode."""
    for chrome in _CHROME_PATHS:
        try:
            result = subprocess.run(
                [chrome, "--headless", "--disable-gpu", "--no-sandbox",
                 f"--print-to-pdf={pdf_path}",
                 f"file://{html_path.resolve()}"],
                capture_output=True, timeout=30,
            )
            if result.returncode == 0 and pdf_path.exists():
                return True
        except (FileNotFoundError, subprocess.TimeoutExpired, OSError):
            continue
    return False


def _try_weasyprint(html_path: Path, pdf_path: Path) -> bool:
    """Try to convert HTML to PDF using weasyprint."""
    try:
        import weasyprint
        weasyprint.HTML(filename=str(html_path)).write_pdf(str(pdf_path))
        return True
    except ImportError:
        return False
    except Exception:
        return False


def export_pdf(html_path: Path, pdf_path: Path = None) -> Path | None:
    """Export HTML report to PDF.

    Tries Chrome headless, then weasyprint, then opens browser.

    Args:
        html_path: Path to the HTML report file.
        pdf_path: Output PDF path (default: same name with .pdf suffix).

    Returns:
        Path to generated PDF, or None if manual export needed.
    """
    html_path = Path(html_path)
    if pdf_path is None:
        pdf_path = html_path.with_suffix(".pdf")
    else:
        pdf_path = Path(pdf_path)

    if not html_path.exists():
        print(f"    ERROR: HTML file not found: {html_path}")
        return None

    # Try Chrome headless
    if _try_chrome_headless(html_path, pdf_path):
        return pdf_path

    # Try weasyprint
    if _try_weasyprint(html_path, pdf_path):
        return pdf_path

    # Fall back to opening browser
    import webbrowser
    webbrowser.open(f"file://{html_path.resolve()}")
    print("    PDF auto-generation not available.")
    print("    Opened HTML in browser — use File > Print > Save as PDF")
    print("    (For automatic PDF: install Chrome or `pip install weasyprint`)")
    return None
