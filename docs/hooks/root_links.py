"""Include repository-root Markdown files and rewrite their links for the site.

A page whose body is ``--8<-- "PATH"`` is replaced by the file at PATH (relative
to the repository root). Links inside it are rewritten relative to the page:
links into ``docs/`` point at the site page; links to root files that the site
also includes point at their page; links to other root files point at GitHub.
"""

import re
from pathlib import Path

INCLUDE = re.compile(r'^--8<-- "([^"]+)"\s*$')
ROOT_LINK = re.compile(r"\]\((?:\./)?docs/")
ROOT_FILE_LINK = re.compile(r"\]\((?:\./)?([A-Z][A-Z_]*\.(?:md|cff))(#[^)]*)?\)")
SITE_PAGES = {  # root file -> site page (relative to docs/)
    "INSTALL.md": "getting-started/installation.md",
    "QUICKSTART.md": "getting-started/quickstart.md",
    "QUICK_REFERENCE.md": "reference/quick-reference.md",
}
REPO_BLOB = "https://github.com/jwestrob/Sharur/blob/main/"


def on_page_markdown(markdown, page, config, **kwargs):
    match = INCLUDE.match(markdown.strip())
    if not match:
        return markdown
    root = Path(config["config_file_path"]).parent
    text = (root / match.group(1)).read_text()
    up = "../" * page.file.src_uri.count("/")
    text = ROOT_LINK.sub("](" + up, text)

    def root_file(m):
        name, anchor = m.group(1), m.group(2) or ""
        target = up + SITE_PAGES[name] if name in SITE_PAGES else REPO_BLOB + name
        return f"]({target}{anchor})"

    return ROOT_FILE_LINK.sub(root_file, text)
