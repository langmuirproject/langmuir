# Needs path for autodoc
import os
import sys
sys.path.insert(0, os.path.abspath('../../langmuir'))
sys.path.insert(0, os.path.abspath('../..'))
sys.path.insert(0, os.path.abspath('../langmuir'))

project = 'Langmuir'
copyright = '2019, Sigvald Marholm and Diako Darian'
author = 'Sigvald Marholm and Diako Darian'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.imgmath',
]

templates_path = ['_templates']
source_suffix = '.rst'
exclude_patterns = []
html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
