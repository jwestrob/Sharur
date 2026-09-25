"""Include repository-root Markdown files and rewrite their links for the site.

A page whose body is ``--8<-- "PATH"`` is replaced by the file at PATH (relative
to the repository root); links into ``docs/`` are rewritten relative to the page.
"""

import re
from pathlib import Path

INCLUDE = re.compile(r'^--8<-- "([^"]+)"\s*$')
ROOT_LINK = re.compile(r"\]\((?:\./)?docs/")


def on_page_markdown(markdown, page, config, **kwargs):
    match = INCLUDE.match(markdown.strip())
    if not match:
        return markdown
    root = Path(config["config_file_path"]).parent
    text = (root / match.group(1)).read_text()
    return ROOT_LINK.sub("](" + "../" * page.file.src_uri.count("/"), text)
