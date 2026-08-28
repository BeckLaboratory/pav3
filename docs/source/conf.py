"""PAV documentation configuration."""

import pathlib
import re
import sys

_DOCS_SOURCE = pathlib.Path(__file__).parent
_REPO_ROOT = _DOCS_SOURCE.parents[1]

# Local extensions (schema table directive).
sys.path.insert(0, str(_DOCS_SOURCE / '_ext'))

project = 'PAV'
copyright = '2025, Peter Audano'
author = 'Peter Audano'


# Version
#
# Read from the source rather than importing so the documented version is the one in the tree
# being built, not whatever happens to be installed in the build environment.

def _read_version() -> str:
    """Extract ``__version__`` from the PAV package source."""
    init_file = _REPO_ROOT / 'src' / 'pav3' / '__init__.py'
    match = re.search(
        r"""^__version__\s*=\s*['"]([^'"]+)['"]""",
        init_file.read_text(),
        re.MULTILINE,
    )

    if match is None:
        raise RuntimeError(f'Could not find __version__ in {init_file}')

    return match.group(1)


release = _read_version()
version = '.'.join(release.split('.')[:2])

extensions = [
    'sphinx.ext.doctest',
    # 'sphinx.ext.autodoc',
    'sphinx.ext.intersphinx',
    'autoapi.extension',
    'sphinxarg.ext',
    'myst_parser',
    'schema_table',
    'sphinx_rtd_theme',
]

templates_path = ['_templates']
exclude_patterns = []

# AutoAPI
autoapi_dirs = [
    '../../src/pav3',
]

# "special-members" and "imported-members" are deliberately omitted: they pull dunder methods and
# re-exported names onto every page, which buries the actual API.
autoapi_options = [
    'members',
    'undoc-members',
    # 'private-members',
    'show-inheritance',
    'show-module-summary',
]

autoapi_keep_files = False

autoapi_member_order = 'groupwise'

# Intersphinx
#
# Each entry maps a package to its published objects.inv, which lets type annotations such as
# pl.DataFrame resolve to the upstream documentation instead of rendering as plain text.
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'polars': ('https://docs.pola.rs/api/python/stable', None),
    'numpy': ('https://numpy.org/doc/stable', None),
}

# Resolving type annotations to upstream documentation
#
# Two gaps stop intersphinx from linking annotations on its own:
#
# 1. Polars does not publish its core classes in objects.inv. Its documentation registers the
#    methods (polars.DataFrame.filter) but not the classes themselves, so "polars.DataFrame" and
#    friends cannot be resolved from the inventory at all. They are mapped to their documentation
#    pages explicitly below.
#
# 2. Some targets are published under a different object type than the one a signature asks for.
#    numpy.typing.NDArray is a "py:data" entry, so a "py:class" cross-reference misses it. Retrying
#    as a generic object matches any type.
#
# Verify the URL map against a nitpicky build after a Polars upgrade:
#     sphinx-build -n -b html docs/source <outdir>

_TYPE_ALIASES = {
    'pl.': 'polars.',
    'np.': 'numpy.',
}

_POLARS_DOCS = 'https://docs.pola.rs/api/python/stable/'

_EXTERNAL_TYPE_URLS = {
    'polars.DataFrame': _POLARS_DOCS + 'reference/dataframe/index.html',
    'polars.LazyFrame': _POLARS_DOCS + 'reference/lazyframe/index.html',
    'polars.Series': _POLARS_DOCS + 'reference/series/index.html',
    'polars.Expr': _POLARS_DOCS + 'reference/expressions/index.html',
    'polars.DataType': _POLARS_DOCS + 'reference/datatypes.html',
}


def _resolve_external_type(app, env, node, contnode):
    """Resolve annotations that intersphinx leaves unresolved.

    :param app: Sphinx application.
    :param env: Build environment.
    :param node: Pending cross-reference that did not resolve.
    :param contnode: Node holding the reference text.

    :return: Resolved reference node, or None to leave the reference unresolved.
    """
    from docutils import nodes
    from sphinx.ext.intersphinx import missing_reference

    if node.get('refdomain') != 'py':
        return None

    target = node.get('reftarget', '')

    for alias, module in _TYPE_ALIASES.items():
        if target.startswith(alias):
            target = module + target[len(alias):]
            node['reftarget'] = target
            break

    url = _EXTERNAL_TYPE_URLS.get(target)

    if url is not None:
        reference = nodes.reference('', '', internal=False, refuri=url, reftitle=target)
        reference += contnode

        return reference

    # Retry as a generic object so entries published as py:data (numpy.typing.NDArray) match a
    # py:class cross-reference.
    if node.get('reftype') != 'obj':
        node['reftype'] = 'obj'

        return missing_reference(app, env, node, contnode)

    return None


def setup(app):
    """Register documentation build hooks.

    :param app: Sphinx application.
    """
    app.connect('missing-reference', _resolve_external_type)


# Theme and build
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
