#!/usr/bin/env python
"""
Copyright 2018
    Sigvald Marholm <marholm@marebakken.com>
    Diako Darian <diakod@math.uio.no>

This file is part of langmuir.

langmuir is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

langmuir is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with langmuir  If not, see <http://www.gnu.org/licenses/>.
"""

from setuptools import setup
from io import open # Necessary for Python 2.7

with open('README.rst', encoding='utf-8') as f:
    long_description = f.read()

with open('.version') as f:
    version = f.read().strip()

setup(name='langmuir',
      version=version,
      description='Functions for work on Langmuir probe theory',
      long_description=long_description,
      author='Sigvald Marholm',
      author_email='marholm@marebakken.com',
      url='https://github.com/sigvaldm/langmuir.git',
      packages=['langmuir'],
      install_requires=['numpy', 'scipy'],
      license='LGPL',
      classifiers=[
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7'
        ]
     )

