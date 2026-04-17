# -*- coding: utf-8 -*-
#
# Sphinx configuration for the TRiSK-based Poisson solver design
# narrative. Independent Sphinx project under the MPAS fork's
# `docs/design/` subtree, deployed as its own Read the Docs project
# separate from the upstream-port docs at `docs/`.

# -- Project --------------------------------------------------------
project   = "A TRiSK-based Poisson solver for MPAS-A"
author    = "David Fillmore"
copyright = "2026, David Fillmore"

# -- General --------------------------------------------------------
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
    ".md":  "markdown",
    ".rst": "restructuredtext",
}

master_doc = "index"

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Bibliography ---------------------------------------------------
bibtex_bibfiles        = ["refs.bib"]
bibtex_default_style   = "plain"
bibtex_reference_style = "author_year"

# -- HTML -----------------------------------------------------------
html_theme = "sphinx_rtd_theme"
html_title = "Poisson solver for MPAS-A — design narrative"

# -- MathJax custom macros ------------------------------------------
# Mirrors the \newcommand definitions in the LaTeX paper preamble so
# that the ported derivations can use them verbatim.
mathjax3_config = {
    "tex": {
        "macros": {
            "bphi":   r"\boldsymbol{\varphi}",
            "brho":   r"\boldsymbol{\rho}",
            "bE":     r"\boldsymbol{E}",
            "bx":     r"\boldsymbol{x}",
            "bn":     r"\boldsymbol{n}",
            "vu":     r"\boldsymbol{u}",
            "vv":     r"\boldsymbol{v}",
            "bb":     r"\boldsymbol{b}",
            "bc":     r"\boldsymbol{c}",
            "Rthree": r"\mathbb{R}^{3}",
            "epsz":   r"\varepsilon_{0}",
            "Aedge":  r"A_{e}",
            "dcen":   r"d_{e}",
            "elen":   r"\ell_{e}",
            "vol":    r"V_{i,k}",
            "Vol":    [r"V_{#1}", 1],
            "dz":     [r"\Delta z_{#1}", 1],
            "hlf":    r"\tfrac{1}{2}",
        }
    }
}
html_static_path  = []
html_show_sphinx  = True
html_show_copyright = True
