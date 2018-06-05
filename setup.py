from setuptools import setup
import os
import sys

_here = os.path.abspath(os.path.dirname(__file__))

if sys.version_info[0] < 3:
    with open(os.path.join(_here, 'README.rst')) as f:
        long_description = f.read()
else:
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
    install_requires=[
        'cffi==1.11.2',
        'click==6.7',
        'pyyaml==3.12',
        'polygons==0.0.0',
        'delaunay==0.0.0',
        'scipy',
    ],
    dependency_links=[
        "git+https://github.com/bast/polygons.git@d7b2c7e706742e194014d3006267dfcf97f80d3e#egg=polygons-0.0.0",
        "git+https://github.com/bast/delaunay.git@radovan/python-interface#egg=delaunay-0.0.0",
    ],
    scripts=['smesh'],
    include_package_data=True,
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3.5',        
        'Programming Language :: Python :: 3.6',
    ],
)
