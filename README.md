[![Build Status](https://travis-ci.org/bast/smeshing.svg?branch=master)](https://travis-ci.org/bast/smeshing/builds)
[![Coverage Status](https://coveralls.io/repos/github/bast/smeshing/badge.svg?branch=master)](https://coveralls.io/github/bast/smeshing?branch=master)


# Simple mesh generator in Python (fork of PyDistMesh)

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


## Installing dependencies

```
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt
```


## Running tests

```
py.test -vv smeshing/test.py
```


## Plotting

Example:
```
python smeshing/plot.py test/circle.txt circle.png
```


## References

The DistMesh algorithm is described in the following two references.
If you use the algorithm in a program or publication, please
acknowledge its authors by adding a reference to the first paper
below.

- [P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB, SIAM Review, Volume 46 (2), pp. 329-345, June 2004](http://persson.berkeley.edu/distmesh/persson04mesh.pdf)
- [P.-O. Persson, Mesh Generation for Implicit Geometries, Ph.D. thesis, Department of Mathematics, MIT, Dec 2004](http://persson.berkeley.edu/thesis/persson-thesis-color.pdf)
