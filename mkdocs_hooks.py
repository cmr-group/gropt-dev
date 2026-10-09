"""MkDocs hooks for the GrOpt docs."""

import posixpath
import re

_IPYNB_HREF = re.compile(r'href="([^":#?]+\.ipynb)(#[^"]*)?"')


def on_page_content(html, page, config, files):
    """Point relative .ipynb links inside a rendered notebook at the linked notebook's page.

    mkdocs-jupyter renders notebooks after MkDocs' own link handling, so a link like
    ``moments/flow_encoding.ipynb`` (correct on GitHub and in Jupyter) would 404 on the site.
    """
    if not page.file.src_uri.endswith('.ipynb'):
        return html

    def fix(m):
        target = posixpath.normpath(posixpath.join(posixpath.dirname(page.file.src_uri), m.group(1)))
        f = files.get_file_from_path(target)
        if f is None:
            return m.group(0)
        return f'href="{f.url_relative_to(page.file)}{m.group(2) or ""}"'

    return _IPYNB_HREF.sub(fix, html)
