# -*- coding: utf-8 -*-
#
# Sphinx configuration for the SHDOM-in-MPAS-A design narrative.
# This is an independent Sphinx project under the MPAS fork's
# `docs/design/` subtree, deployed as its own Read the Docs project
# separate from the upstream-port docs at `docs/`.

# -- Project --------------------------------------------------------
project   = "Embedding SHDOM in MPAS-A"
author    = "David Fillmore"
copyright = "2026, David Fillmore"

# -- General --------------------------------------------------------
extensions = [
    "myst_parser",
    "sphinx.ext.mathjax",
    "sphinxcontrib.bibtex",
]

myst_enable_extensions = [
    "dollarmath",   # $...$ inline and $$...$$ display math
    "amsmath",      # \begin{equation}, \begin{align}, ...
    "colon_fence",  # ::: directives
    "deflist",      # definition lists
]

source_suffix = {
    ".md":  "markdown",
    ".rst": "restructuredtext",
}

master_doc = "index"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Bibliography ---------------------------------------------------
bibtex_bibfiles       = ["refs.bib"]
bibtex_default_style  = "plain"
bibtex_reference_style = "author_year"

# -- HTML -----------------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_title = "SHDOM in MPAS-A — design narrative"
html_static_path  = []
html_show_sphinx  = True
html_show_copyright = True
