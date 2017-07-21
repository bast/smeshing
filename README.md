[![Build Status](https://travis-ci.org/bast/smeshing.svg?branch=master)](https://travis-ci.org/bast/smeshing/builds)
[![Coverage Status](https://coveralls.io/repos/github/bast/smeshing/badge.svg?branch=master)](https://coveralls.io/github/bast/smeshing?branch=master)


# Simple mesh generator in Python (fork of PyDistMesh)

<img src="https://github.com/bast/smeshing/raw/master/img/example.jpg" width="600">

This version of the code is based on
[PyDistMesh](https://github.com/bfroehle/pydistmesh) developed by [Bradley M.
Froehle](https://github.com/bfroehle).
This code generates unstructured triangular and tetrahedral meshes using signed
distance functions.
Like [DistMesh](http://persson.berkeley.edu/distmesh/) and
[PyDistMesh](https://github.com/bfroehle/pydistmesh) upon which this work is
based, this code is distributed under the [GNU GPL](../master/LICENSE).


## Status and roadmap for this fork

Under heavy development. The plan is to be able to do large grids. We wish to
use the code to prepare grids for simulations on the Norwegian coast with
thousands of islands (polygons).

Prototyping work is currently done in Python. Later we will probably move to
Fortran or C(++) and introduce parallelization but provide a Python interface.

Currently the code is a bit hardcoded for islands and coastline and specific
resolution functions but later the code will be generalized and the hardcoded
aspects will be abstracted out.


## Installing dependencies for development

```
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt
```


## Running tests

```
py.test -vv smeshing/*.py
```


## Plotting

Example:
```
python smeshing/plot.py data/fiction/result.txt example.png
```


## Example installation and run

```
virtualenv venv
source venv/bin/activate

pip install git+https://github.com/bast/polygons.git
pip install git+https://github.com/bast/flanders.git
pip install git+https://github.com/bast/delaunay.git@radovan/python-interface
pip install git+https://github.com/bast/smeshing.git

smesh --help
```


## Installation on [Stallo](https://www.sigma2.no/content/stallo)

```
module load foss/2016b
module load CMake
module load libffi

cd ${SLURM_SUBMIT_DIR}

virtualenv venv
source venv/bin/activate

export CC=gcc
export CXX=g++
export FC=gfortran

   pip install git+https://github.com/bast/polygons.git \
&& pip install git+https://github.com/bast/flanders.git \
&& pip install git+https://github.com/bast/delaunay.git@radovan/python-interface \
&& pip install git+https://github.com/bast/smeshing.git
```


## Pros

- Individual components live in separate libraries.
- A lot of effort was invested in avoiding quadratic scaling.
- Optimization is fully relaxed - no interpolation is done.
- Delaunay is performed at every step.
- Good memory profile.


## Known issues

- Code uses shared-memory parallelization but the load leveling is not optimal and the scaling on Stallo has not been studied.
- Sometimes a step generates crazy meshes. Possibly there is a numerical instability for small forces or division by very small number.
- One needs to set a `seeding_speed`. This controls in how many steps the initial point distribution is set up.
  The higher the number, the fewer steps, and the more uniform distribution. The more steps, the more time it takes, but also the more
  it will reflect the resolution function. Currently it needs some experimentation. Later we need something more black box.
- Coastline resolution is currently hard-coded to be inversely proportional to the distance to nearest neighbor in view divided by 6.
  In future versions this will be read from input.
- Currently no lower and upper bounds on resolution can be set.
- There is no stop criterion, it will run as many iterations as you ask it to.


## Restart

It is possible to restart a calculation if you provide `--restart=/path/to/restart/file`.


## References

The DistMesh algorithm is described in the following two references.
If you use the algorithm in a program or publication, please
acknowledge its authors by adding a reference to the first paper
below.

- [P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB, SIAM Review, Volume 46 (2), pp. 329-345, June 2004](http://persson.berkeley.edu/distmesh/persson04mesh.pdf)
- [P.-O. Persson, Mesh Generation for Implicit Geometries, Ph.D. thesis, Department of Mathematics, MIT, Dec 2004](http://persson.berkeley.edu/thesis/persson-thesis-color.pdf)
