# Configuration file for the Sphinx documentation builder.
# Build: sphinx-build -b html docs/ docs/_build/html

import os
import sys

# -- Project information -----------------------------------------------------
project   = 'hftest'
copyright = '2026, Wang Shiyu'
author    = 'Wang Shiyu'
release   = '0.1.0'

# -- General configuration ---------------------------------------------------
extensions = [
    'myst_parser',          # Markdown support (used for supplementary pages)
    'sphinx.ext.mathjax',   # LaTeX-rendered equations in browser
    'sphinx.ext.autosectionlabel',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',  # "View source" link stubs
]

source_suffix = {
    '.rst': 'restructuredtext',
    '.md':  'markdown',
}

templates_path   = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- MathJax -----------------------------------------------------------------
mathjax3_config = {
    'tex': {
        'packages': {'[+]': ['ams']},
        'macros': {
            'bra': [r'\langle #1 |', 1],
            'ket': [r'| #1 \rangle', 1],
            'braket': [r'\langle #1 | #2 \rangle', 2],
        },
    }
}

# -- Cross-referencing -------------------------------------------------------
autosectionlabel_prefix_document = True
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
}

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'navigation_depth':  4,
    'titles_only':       False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': True,
}
html_static_path    = ['_static']
html_show_sourcelink = False
html_show_sphinx     = True

# -- Options for LaTeX (PDF) output ------------------------------------------
latex_elements = {
    'preamble': r'\usepackage{amsmath}\usepackage{amssymb}',
}
