"""Markdown -> HTML converter (no external dependencies)."""

import re
from pathlib import Path


_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>{title}</title>
<style>
  :root {{ --bg: #fff; --fg: #1a1a2e; --accent: #2563eb; --border: #e2e8f0;
           --code-bg: #f1f5f9; --table-stripe: #f8fafc; --warn: #dc2626; }}
  @media (prefers-color-scheme: dark) {{
    :root {{ --bg: #0f172a; --fg: #e2e8f0; --accent: #60a5fa; --border: #334155;
             --code-bg: #1e293b; --table-stripe: #1e293b; --warn: #f87171; }}
  }}
  *, *::before, *::after {{ box-sizing: border-box; }}
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
          line-height: 1.7; color: var(--fg); background: var(--bg);
          max-width: 54em; margin: 0 auto; padding: 2em 1.5em; }}
  h1 {{ border-bottom: 2px solid var(--accent); padding-bottom: .3em; }}
  h2 {{ border-bottom: 1px solid var(--border); padding-bottom: .2em; margin-top: 2em; }}
  h3 {{ margin-top: 1.5em; }}
  hr {{ border: none; border-top: 1px solid var(--border); margin: 2em 0; }}
  code {{ background: var(--code-bg); padding: .15em .35em; border-radius: 3px;
          font-size: .9em; }}
  table {{ border-collapse: collapse; width: 100%; margin: 1em 0; }}
  th, td {{ border: 1px solid var(--border); padding: .5em .75em; text-align: left; }}
  th {{ background: var(--code-bg); font-weight: 600; }}
  tr:nth-child(even) {{ background: var(--table-stripe); }}
  ul, ol {{ padding-left: 1.5em; }}
  li {{ margin: .3em 0; }}
  strong {{ color: var(--accent); }}
  .warn {{ color: var(--warn); font-weight: bold; }}
  @media print {{ body {{ max-width: 100%; font-size: 11pt; }} }}
</style>
</head>
<body>
{body}
</body>
</html>
"""


def _md_to_html(md: str) -> str:
    """Convert report Markdown to HTML. Handles the subset used by our reports."""
    lines = md.split("\n")
    html_parts = []
    i = 0

    while i < len(lines):
        line = lines[i]

        # Blank line
        if not line.strip():
            i += 1
            continue

        # Horizontal rule
        if re.match(r"^-{3,}\s*$", line) or re.match(r"^\*{3,}\s*$", line):
            html_parts.append("<hr>")
            i += 1
            continue

        # Headers
        m = re.match(r"^(#{1,6})\s+(.*)", line)
        if m:
            level = len(m.group(1))
            text = _inline(m.group(2))
            html_parts.append(f"<h{level}>{text}</h{level}>")
            i += 1
            continue

        # Table (starts with |)
        if line.strip().startswith("|"):
            table_lines = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                table_lines.append(lines[i])
                i += 1
            html_parts.append(_table_to_html(table_lines))
            continue

        # Unordered list
        if re.match(r"^\s*[-*]\s+", line):
            items = []
            while i < len(lines) and re.match(r"^\s*[-*]\s+", lines[i]):
                items.append(re.sub(r"^\s*[-*]\s+", "", lines[i]))
                i += 1
            html_parts.append("<ul>" + "".join(f"<li>{_inline(it)}</li>" for it in items) + "</ul>")
            continue

        # Paragraph (collect consecutive non-empty, non-special lines)
        para = []
        while i < len(lines) and lines[i].strip() and not re.match(r"^(#{1,6}\s|[-*]{3,}\s*$|\||\s*[-*]\s+)", lines[i]):
            para.append(lines[i])
            i += 1
        if para:
            html_parts.append(f"<p>{_inline(' '.join(para))}</p>")

    return "\n".join(html_parts)


def _inline(text: str) -> str:
    """Convert inline Markdown: bold, code, italic."""
    text = re.sub(r"`([^`]+)`", r"<code>\1</code>", text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", text)
    return text


def _table_to_html(lines: list[str]) -> str:
    """Convert Markdown table lines to HTML table."""
    rows = []
    for line in lines:
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        rows.append(cells)

    # Skip separator row (e.g. |---|---|)
    data_rows = [r for r in rows if not all(re.match(r"^[-:]+$", c) for c in r)]
    if not data_rows:
        return ""

    html = "<table>\n<thead><tr>"
    for cell in data_rows[0]:
        html += f"<th>{_inline(cell)}</th>"
    html += "</tr></thead>\n<tbody>\n"
    for row in data_rows[1:]:
        html += "<tr>" + "".join(f"<td>{_inline(c)}</td>" for c in row) + "</tr>\n"
    html += "</tbody></table>"
    return html


def _write_html(md_content: str, md_path: Path):
    """Write an HTML version alongside the Markdown file."""
    html_path = md_path.with_suffix(".html")
    title = "Genetic Health Report"
    # Extract title from first H1
    m = re.search(r"^#\s+(.+)", md_content, re.MULTILINE)
    if m:
        title = m.group(1)
    body = _md_to_html(md_content)
    html = _HTML_TEMPLATE.format(title=title, body=body)
    with open(html_path, "w") as f:
        f.write(html)
    print(f"    HTML:     {html_path}")
