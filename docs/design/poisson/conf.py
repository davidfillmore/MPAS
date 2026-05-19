# Configuration for the Poisson design-narrative Sphinx tree.

project = "Poisson solver for MPAS-A"
author = "David Fillmore"
copyright = "2026, David Fillmore"

master_doc = "index"

extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
]

myst_enable_extensions = [
    "dollarmath",
    "amsmath",
    "colon_fence",
    "deflist",
]

source_suffix = {
    ".md": "markdown",
    ".rst": "restructuredtext",
}

bibtex_bibfiles = ["refs.bib"]
bibtex_default_style = "plain"

html_theme = "sphinx_rtd_theme"
html_title = project
