# File: docs/conf.py

import os
import sys

from multiproject.utils import get_project

sys.path.insert(0, os.path.abspath('.'))
sys.path.insert(0, os.path.abspath('../doc/breathe'))
sys.path.insert(0, os.path.abspath('../doc/breathe'))
sys.path.insert(0, os.path.abspath('../interfaces/Python/'))
sys.path.insert(0, os.path.abspath('../interfaces/Python/'))


# -- Project information -----------------------------------------------------

project = 'ViennaRNA'
copyright = '1994 - 2023, Ronny Lorenz, Ivo L. Hofacker, et al.'
author = 'Ronny Lorenz, Ivo L. Hofacker, et al.'

# The full version, including alpha/beta/rc tags
release = '2.6.4'

# -- General configuration ---------------------------------------------------

# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.
extensions = [
      'multiproject',
      'sphinx.ext.napoleon',
      'sphinx.ext.autodoc',
      'sphinx.ext.mathjax',
      'sphinx_rtd_theme',
      'myst_parser',
]

# Define the projects that will share this configuration file.
multiproject_projects = {
    "master": {
        "use_config_file": False,
        "path": "source",
        "config": {
            "html_title" : "ViennaRNA Package",
            "html_favicon" : "source/gfx/vrna_32.png",
            "html_css_files" : ['css/custom.css',],
        },
    },
    "python": {
        "use_config_file": False,
        "path": "python",
        "config": {
            "html_title" : "ViennaRNA Package - Python",
        },
    },
}

current_project = get_project(multiproject_projects)

if current_project == "master":
    extensions += [
      'sphinx.ext.todo',
      "sphinx.ext.autosectionlabel",
      "sphinx.ext.imgconverter",
      'sphinxcontrib.bibtex',
      'breathe',
      'sphinx_copybutton',
    ]
elif current_project == "python":
    extensions += [
      'sphinx.ext.viewcode',
      'sphinx.ext.autosummary',
      'sphinx.ext.graphviz',
    ]

# Common options.
mathjax_path = "js/mathjax/tex-chtml.js"
mathjax3_config = {
    'tex': {
        'inlineMath': [['$', '$'], ['\\(', '\\)']],
        'displayMath': [["\\[", "\\]"]]
    }
}

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_include_private_with_doc = False
napoleon_include_special_with_doc = True
napoleon_use_admonition_for_examples = False
napoleon_use_admonition_for_notes = False
napoleon_use_admonition_for_references = False
napoleon_use_ivar = False
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True
napoleon_type_aliases = None
napoleon_attr_annotations = True
napoleon_type_aliases = {
    "PRIVATE int": "int",
    "PRIVATE FLT_OR_DBL" : "double",
    "unsigned int *": "list-like(unsigned int)",
    "unsigned int **": "list-like(list-like(unsigned int))",
    "short *": "list-like(int)",
    "char *": "string",
    "const char *": "string",
    "float *": "list-like(double)",
    "double *": "list-like(double)",
    "double **": "list-like(list-like(double))",
    "vrna_fold_compound_t *" : "fold_compound",
    "vrna_param_t *" : "param",
    "vrna_exp_param_t *" : "exp_param",
    "vrna_md_t *" : "md",
    "std::string": "string",
    "FLT_OR_DBL" : "double",
    "FLT_OR_DBL *" : "list-like(double)",
    "std::vector<FLT_OR_DBL>" : "list-like(double)"
}


# Add any paths that contain templates here, relative to this directory.
templates_path = [f"{multiproject_projects[current_project]['path']}/_templates"]

# The encoding of source files.
source_encoding = "utf-8"

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = []

# Make sure the target is unique
autosectionlabel_prefix_document = True

pygments_style = 'default'
option_emphasise_placeholders = True

# -- Options for HTML output -------------------------------------------------

#html_logo = "gfx/vrna_logo.png"

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'logo_only': False,
    'display_version': False,
}

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = [f"{multiproject_projects[current_project]['path']}/_static"]
html_show_sourcelink = False
html_context = {
    "sidebar_external_links_caption": "Links",
    "sidebar_external_links": [
        (
            '<i class="fa fa-github fa-fw"></i> Source code',
            "https://github.com/ViennaRNA/ViennaRNA",
        ),
        (
            '<i class="fa fa-bug fa-fw"></i> Issue tracker',
            "https://github.com/ViennaRNA/ViennaRNA/issues",
        ),
    ],
}

latex_elements = {
    # The paper size ('letterpaper' or 'a4paper').
    'papersize': 'a4paper',

    # The font size ('10pt', '11pt' or '12pt').
    #'pointsize': '10pt',

    # Additional stuff for the LaTeX preamble.
    #'preamble':'',

    #Figure placement within LaTeX paper NOT WORKING
    'figure_align': 'H',
}

# Bibtex settings
bibtex_bibfiles = ['../viennarna.bib']
bibtex_default_style = 'unsrt'
bibtex_reference_style = 'author_year'

# Breathe settings
breathe_projects = {"ViennaRNA": "doxygen/xml/"}
breathe_default_project = "ViennaRNA"
breathe_domain_by_extension = {
      "h" : "c",
      "c" : "c",
      "inc" : "c"
  }
breathe_show_define_initializer = False
breathe_show_include = True
breathe_show_decl_file_include = True
breathe_default_members = ('members', 'undoc-members')
