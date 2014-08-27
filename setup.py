#!/usr/bin/env python

# Setup script for bestmsm package

import os
from setuptools import setup,find_packages

def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
		name='BestMSM',
		version='0.2dev',
		description='the Best way to make an msm',
		url='http://github.com/daviddesancho/bestmsm',
		author='David De Sancho',
		author_email='daviddesancho.at.gmail.com',
		license='GPL',
		packages=find_packages(),
		keywords= "markov state model",
		long_description=read('README.md'),
		classifiers = ["""\
				Development Status :: 1 - Planning
				Operating System :: POSIX :: Linux
				Operating System :: MacOS
				Programming Language :: Python :: 2.7
				Topic :: Scientific/Engineering :: Bio-Informatics
				Topic :: Scientific/Engineering :: Chemistry
				"""]
		)
