# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
import os

from pathlib import Path

# Generate a small file containing the autodoc/autoclass directive
# for one particular class. This file can then be included from
# the .rst files of this documentation. The detour is needed
# because the class names only become available while conf.py itself
# runs (via the JIT compiler we invoke below). Therefore the
# class names cannot be hard-coded directly in the .rst files.
def generateSnippet(path, file_name_without_suffix, raw_class_name):
    # The raw class name (as given by type(...)) is of the form <class '...'>.
    # To be able to give this to autodoc as a string, we have to remove everything
    # except what is within the quotation marks.
    class_name = raw_class_name.split("'")[1]

    (path / (file_name_without_suffix + ".snippet")).write_text(
        f".. autoclass:: {class_name}\n"
        "   :members:\n"
    )

# Large parts of the Python interface is just-in-time compiled
# from C++ code. Therefore, before we can let Sphinx document
# the API, we have to trigger to JIT-compiler to actually
# build the API (yes, we actually call Dune here!)

# A grid object
from dune.grid import structuredGrid
import dune.functions as functions

grid = structuredGrid([0,0],[1,1],[1,1])
basis = functions.defaultGlobalBasis(grid, functions.Lagrange(order=2))
lv = basis.localView()

print(lv.cppIncludes)

# Generate a directory for auto-generated reST files
generated_dir = Path(__file__).parent / "_generated"
generated_dir.mkdir(exist_ok=True)

# Generate short reST files that contain the autoclass directive
# for the various API classes. These are later included from the
# actual documentation reST files.
generateSnippet(generated_dir, "globalbasis", str(type(basis)))
generateSnippet(generated_dir, "localview", str(type(lv)))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'dune-functions Python API'
copyright = '2026, The Dune core developers'
author = 'The Dune core developers'
release = '2.12'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    'sphinx.ext.autodoc',
    "sphinx.ext.todo",
]

todo_include_todos = True

add_module_names = False
autodoc_typehints_format = "short"

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
}

templates_path = ['_templates']
exclude_patterns = []



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
