# Fork of PyDistMesh: A Simple Mesh Generator in Python

This version of the code is based on
[PyDistMesh](https://github.com/bfroehle/pydistmesh) developed by [Bradley M.
Froehle](https://github.com/bfroehle).
This code generates unstructured triangular and tetrahedral meshes using signed
distance functions.
Like [DistMesh](http://persson.berkeley.edu/distmesh/) and
[PyDistMesh](https://github.com/bfroehle/pydistmesh) upon which this work is
based, this code is distributed under the [GNU GPL](../master/LICENSE).


## Installing the code and its dependencies

```
virtualenv venv
source venv/bin/activate
pip install -r requirements.txt
python setup.py install
```


## Running tests

```
python distmesh/demo2d.py
```

If the tests fail, it will stop with an AssertionError. Soon we will switch to
pytest (and rename files) but for this we need to get rid of Cython and switch to CFFI.


## References

The DistMesh algorithm is described in the following two references.
If you use the algorithm in a program or publication, please
acknowledge its authors by adding a reference to the first paper
below.

- P.-O. Persson, G. Strang, **A Simple Mesh Generator in MATLAB**.
  *SIAM Review*, Volume 46 (2), pp. 329-345, June 2004 [PDF](http://persson.berkeley.edu/distmesh/persson04mesh.pdf)

- P.-O. Persson, **Mesh Generation for Implicit Geometries**.
  Ph.D. thesis, *Department of Mathematics, MIT*, Dec 2004 [PDF](http://persson.berkeley.edu/thesis/persson-thesis-color.pdf)