#! /usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from setuptools import setup, find_namespace_packages

if sys.argv[-1] == "publish":
    os.system("python setup.py sdist upload")
    sys.exit()

setup(
    name='ntsz',
    version='1.2',
    author='Vyoma Muralidhara',
    author_email='vyoma1993@gmail.com',
    description='A python 3 implementation for calculation of various Sunyaev\
        -Zeldovich effects.',
    long_description=open('README.md').read(),
    packages=find_namespace_packages(),
    package_data={
        "ntsz": ["LICENSE"],
        "ntsz.data": [
            "ntsz_grid_single.fits",
            "ntsz_grid_broken.fits",
            "planck_bandpass.fits"
        ]
    },
    include_package_data=True,
    long_description_content_type='text/markdown',
    url='https://github.com/Vyoma-M/ntsz',
    classifiers=[
        'Programming Language :: Python :: 3',
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=[
        'healpy',
        'numpy',
        'scipy',
        'astropy',
    ],
    test_suite='tests',
    zip_safe=False,
)
