#!/usr/bin/env python

#import os
import setuptools
from distutils.core import setup

setup(
  name='iRep',
  packages = ['iRep'],
  scripts = ['bin/iRep', 'bin/bPTR', 'bin/gc_skew', 'iRep/iRep_filter.py'],
  version='1.1.3',
  description='calculate iRep replication rates from metagenome sequencing',
  author='Chris Brown',
  author_email='ctb@berkeley.edu',
  license='MIT',
  url='https://github.com/christophertbrown/iRep',
  install_requires=['lmfit', 'numpy', 'scipy', 'pandas', \
                    'seaborn', 'matplotlib'],
  classifiers=[
          'Programming Language :: Python',
          'Programming Language :: Python :: 3'
  ]
)
