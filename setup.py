#!/usr/bin/env python

from distutils.core import setup

setup(name='smeshing',
      version='0.0.0',
      description='Simple mesh generator.',
      author='Radovan Bast',
      author_email='radovan.bast@uit.no',
      url='https://github.com/bast/smeshing',
      scripts=['smesh'],
      packages=['smeshing'],
      license='GPL-3.0',
      install_requires=['cffi', 'click', 'pyyaml'])
