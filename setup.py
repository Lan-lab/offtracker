#!/usr/bin/env python
# -*- coding: utf-8 -*-

import io
import os
import sys
from shutil import rmtree
from setuptools import find_packages, setup, Command

#
NAME = 'offtracker'
DESCRIPTION = 'Tracking-seq data analysis'
AUTHOR = 'Runda Xu'
EMAIL = 'xrd18@tsinghua.org.cn'
URL = 'https://github.com/Lan-lab/offtracker'
REQUIRES_PYTHON = '>=3.6.0'

here = os.path.abspath(os.path.dirname(__file__))

package_folder = NAME.lower().replace("-", "_").replace(" ", "_")
with open(os.path.join(here, package_folder, '_version.py'),'r',encoding='utf-8') as f:
    for line in f:
        if line.startswith("__version__"):
            VERSION = line.strip().split("=")[1].strip().replace('"', '')
            break

# requirements
REQUIRED = [
    'biopython', 'pybedtools', 'pyyaml', 'pandas', 'numpy',
]
## pybedtools may be not supported in Windows

try:
    with io.open(os.path.join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION


setup(
    name=NAME,
    version=VERSION,
    description=DESCRIPTION,
    long_description=long_description,
    long_description_content_type='text/markdown',
    author=AUTHOR,
    author_email=EMAIL,
    url=URL,
    python_requires=REQUIRES_PYTHON,
    packages=['offtracker'],
    package_data={'offtracker': ['snakefile/*','utility/*']},
    scripts = ['scripts/offtracker_qc.py',
               'scripts/offtracker_config.py',
               'scripts/offtracker_candidates.py',
               'scripts/offtracker_analysis.py',
               'scripts/offtracker_plot.py'],
    install_requires=REQUIRED,
    include_package_data=True
)
