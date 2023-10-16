#!/usr/bin/env python

# Copyright (C) 2015 University of California, Los Angeles (UCLA)
# Yu-Ting Tseng, Yi Xing
#
# Authors: Yu-Ting Tseng, Yi Xing
#
# This program is free software: you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.

# This program is distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

# You should have received a copy of the GNU General Public License along with
# this program. If not, see http://www.gnu.org/licenses/.

from os import path
from setuptools import setup, find_packages

# Get the long description from the README file
here = path.abspath(path.dirname(__file__))
with open(path.join(here, 'README.md')) as f:
    long_description = f.read()

setup(
    name='rmats2sashimiplot',
    version='3.0.0',
    # 'find_packages' finds just the package name without the full path.
    # 'package_dir' is needed in combination with 'packages' to get the
    # full path of each package.
    packages=find_packages('src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'rmats2sashimiplot=rmats2sashimiplot.rmats2sashimiplot:main',
            'index_gff=MISO.misopy.index_gff:main',
            'sashimi_plot=MISO.misopy.sashimi_plot.sashimi_plot:main'
        ]
    },
    description='rmats2sashimiplot',
    long_description=long_description,
    long_description_content_type='text/markdown',
    licence='GNU GPLv2',

    author='Yu-Ting Tseng, Emad Bahrami-Samani, Zhijie Xie, Yukai Jiang',
    author_email='shiehshiehzhijie@gmail.com',
    url='https://github.com/Xinglab/rmats2sashimiplot',
    download_url='https://github.com/Xinglab/rmats2sashimiplot',
    classifiers=[
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: GNU General Public License v2 (GPLv2)'
    ],
)
