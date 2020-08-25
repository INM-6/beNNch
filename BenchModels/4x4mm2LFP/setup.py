#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
setup.py for mesocircuit_analysis.

Copyright (C) 2014-2018, 4x4mm2LFP model authors.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

'''

from setuptools import setup, Extension
import numpy
from Cython.Distutils import build_ext
cmdclass = { 'build_ext' : build_ext}
ext_modules = [
    Extension('mesocircuit_analysis.helperfun', 
    ['mesocircuit_analysis/helperfun.pyx'],
    include_dirs=[numpy.get_include()]),
    ]


with open('README.md') as file:
    long_description = file.read()


setup(
    name = 'mesocircuit_analysis',
    version = '0.1',
    maintainer = ['Espen Hagen', 'Johanna Senk'],
    maintainer_email = ['espehage@fys.uio.no', 'j.senk@fz-juelich.de'],
    url = 'github.com/INM6',
    packages = ['mesocircuit_analysis'],
    provides = ['mesocircuit_analysis'],
    cmdclass = cmdclass, 
    ext_modules = ext_modules,
    description = 'methods to analyze cortical mesocircuit network model',
    long_description = long_description,
    license='LICENSE',
)
