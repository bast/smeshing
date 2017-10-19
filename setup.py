from setuptools import setup
import os

_here = os.path.abspath(os.path.dirname(__file__))

with open(os.path.join(_here, 'README.rst'), encoding='utf-8') as f:
    long_description = f.read()

version = {}
with open(os.path.join(_here, 'smeshing', 'version.py')) as f:
    exec(f.read(), version)

setup(
    name='smeshing',
    version=version['__version__'],
    description=('Simple mesh generator.'),
    long_description=long_description,
    author='Radovan Bast',
    author_email='radovan.bast@uit.no',
    url='https://github.com/bast/smeshing',
    license='GPL-3.0',
    packages=['smeshing'],
    install_requires=['cffi', 'click', 'pyyaml'],
#   install_requires=[
#       'click==6.7',
#       'numpy==1.13.1',
#       'pyyaml==3.12',
#   ],
#   scripts=['bin/cooley'],
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.6'],
    )
