# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys
sys.path.insert(0, os.path.abspath('..'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SLOTH'
copyright = '2023, Niklas WAGNER'
author = 'Niklas WAGNER'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx.ext.autodoc',
        'sphinx.ext.viewcode',
        'sphinx.ext.napoleon',
        'sphinx.ext.mathjax',
        'sphinx_copybutton',
        'myst_nb',
        ]

# nb settings
# Tun off execution of notebooks. 
# You can rendere the notebook by running them, but this is more complex as you
# have to manage different kernals. If you do not run the noteboke while
# rendering, sphinx is taking the last results from the last execution with a
# jupyter-notebook environment, stored witht the .ipnb file. 
nb_execution_mode = "off"
# If you still want to run your notebook while rendering, you will run into the 
# problem of different kernals. Different contributors use different Jupyter 
# notebook environments and thus different kernals. 
# [MyST-NB](https://myst-nb.readthedocs.io/en/latest/computation/execute.html#execute-intro) 
# can map different kernals to a single one using regex. So uncomment the line 
# below to force MyST-NB to use the kernel named `python3` (most likely the 
# default) instead of any other kernel name (`.*`).
#nb_kernel_rgx_aliases = {'.*': 'python3'}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
source_suffix = ['.rst', '.md']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
