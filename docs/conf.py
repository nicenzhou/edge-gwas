# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
project = 'edge-gwas'
copyright = '2025, edge-gwas Contributors'
author = 'Jiayan Zhou'

# The full version, including alpha/beta/rc tags
release = '0.1.2'
version = '0.1.2'

# -- General configuration ---------------------------------------------------
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
    'sphinx.ext.mathjax',
    'sphinx.ext.todo',
    'sphinx.ext.coverage',
    'sphinx.ext.ifconfig',
    'sphinx_rtd_theme',
    'myst_parser',
]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# -- Options for HTML output -------------------------------------------------
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

html_theme_options = {
    'logo_only': False,
    'display_version': True,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'style_nav_header_background': '#2980B9',
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}

# Custom CSS
html_css_files = [
    'custom.css',
]

# HTML context
html_context = {
    "display_github": True,
    "github_user": "nicenzhou",
    "github_repo": "edge-gwas",
    "github_version": "main",
    "conf_py_path": "/docs/",
}

# Add custom sidebar
html_sidebars = {
    '**': [
        'globaltoc.html',
        'relations.html',
        'sourcelink.html',
        'searchbox.html',
    ]
}

# Logo and favicon (if you have them)
# html_logo = '_static/logo.png'
# html_favicon = '_static/favicon.ico'

# -- Napoleon settings -------------------------------------------------------
napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = True
napoleon_use_admonition_for_notes = True
napoleon_use_admonition_for_references = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_use_ivar = True
napoleon_preprocess_types = True
napoleon_type_aliases = None
napoleon_attr_annotations = True

# -- Intersphinx mapping -----------------------------------------------------
intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable/', None),
    'pandas': ('https://pandas.pydata.org/docs/', None),
    'scipy': ('https://docs.scipy.org/doc/scipy/', None),
    'sklearn': ('https://scikit-learn.org/stable/', None),
    'matplotlib': ('https://matplotlib.org/stable/', None),
    'statsmodels': ('https://www.statsmodels.org/stable/', None),
}

# -- Source file parsers -----------------------------------------------------
source_suffix = {
    '.rst': 'restructuredtext',
    '.md': 'markdown',
}

# -- MyST settings -----------------------------------------------------------
myst_enable_extensions = [
    "colon_fence",
    "deflist",
    "dollarmath",
    "fieldlist",
    "html_admonition",
    "html_image",
    "replacements",
    "smartquotes",
    "strikethrough",
    "substitution",
    "tasklist",
]

# -- Autodoc settings --------------------------------------------------------
autodoc_default_options = {
    'members': True,
    'member-order': 'bysource',
    'special-members': '__init__',
    'undoc-members': True,
    'exclude-members': '__weakref__',
    'show-inheritance': True,
}

autodoc_typehints = 'description'
autodoc_typehints_description_target = 'documented'

# -- Todo extension settings -------------------------------------------------
todo_include_todos = True

# -- Math settings -----------------------------------------------------------
mathjax_path = 'https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-mml-chtml.js'

mathjax3_config = {
    'tex': {
        'inlineMath': [['$', '$'], ['\\(', '\\)']],
        'displayMath': [['$$', '$$'], ['\\[', '\\]']],
    },
}

# -- Options for LaTeX output ------------------------------------------------
latex_elements = {
    'papersize': 'letterpaper',
    'pointsize': '10pt',
    'preamble': r'''
        \usepackage{amsmath}
        \usepackage{amssymb}
    ''',
}

latex_documents = [
    ('index', 'edge-gwas.tex', 'edge-gwas Documentation',
     'Jiayan Zhou', 'manual'),
]

# -- Options for manual page output ------------------------------------------
man_pages = [
    ('index', 'edge-gwas', 'edge-gwas Documentation',
     [author], 1)
]

# -- Options for Texinfo output ----------------------------------------------
texinfo_documents = [
    ('index', 'edge-gwas', 'edge-gwas Documentation',
     author, 'edge-gwas', 
     'Efficient Detection of Genetic Effects in Genome-Wide Association Studies',
     'Miscellaneous'),
]

# -- Options for Epub output -------------------------------------------------
epub_title = project
epub_author = author
epub_publisher = author
epub_copyright = copyright

# Epub options
epub_exclude_files = ['search.html']

# -- Extension configuration -------------------------------------------------

# Add any paths that contain custom static files (such as style sheets)
# html_static_path = ['_static']

# Add version info
rst_prolog = f"""
.. |version| replace:: {version}
.. |release| replace:: {release}
"""

# Add common substitutions
rst_epilog = """
.. _GitHub: https://github.com/nicenzhou/edge-gwas
.. _Issues: https://github.com/nicenzhou/edge-gwas/issues
.. _PyPI: https://pypi.org/project/edge-gwas/
.. _PLINK2: https://www.cog-genomics.org/plink/2.0/
.. _GCTA: https://yanglab.westlake.edu.cn/software/gcta/
.. _GENESIS: https://bioconductor.org/packages/GENESIS/
"""

# -- Options for linkcheck ---------------------------------------------------
linkcheck_ignore = [
    r'http://localhost:\d+/',
]

# -- Nitpicky mode -----------------------------------------------------------
# Be strict about warnings
nitpicky = True

# Suppress certain warnings
nitpick_ignore = [
    ('py:class', 'optional'),
    ('py:class', 'callable'),
]

# -- Version warning ---------------------------------------------------------
# Add version warning for older versions
if version != '0.1.2':
    html_theme_options['announcement'] = (
        f'<b>Warning:</b> You are viewing documentation for version {version}. '
        f'The latest version is 0.1.2. '
        f'<a href="https://edge-gwas.readthedocs.io/en/latest/">View latest docs</a>.'
    )

# -- Custom setup ------------------------------------------------------------
def setup(app):
    """Custom Sphinx setup."""
    # Add custom CSS
    app.add_css_file('custom.css')
    
    # Add custom JavaScript (if needed)
    # app.add_js_file('custom.js')
    
    # Add custom domains (if needed)
    # app.add_domain(MyCustomDomain)
