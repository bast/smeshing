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


## 2-D Examples

- Uniform Mesh on Unit Circle:
```python
>>> import distmesh as dm
>>> import numpy as np
>>> fd = lambda p: np.sqrt((p**2).sum(1))-1.0
>>> p, t = dm.distmesh2d(fd, dm.huniform, 0.2, (-1,-1,1,1))
```

- Rectangle with circular hole, refined at circle boundary:
```python
>>> import distmesh as dm
>>> fd = lambda p: dm.ddiff(dm.drectangle(p,-1,1,-1,1),
...                         dm.dcircle(p,0,0,0.5))
>>> fh = lambda p: 0.05+0.3*dm.dcircle(p,0,0,0.5)
>>> p, t = dm.distmesh2d(fd, fh, 0.05, (-1,-1,1,1),
...                      [(-1,-1),(-1,1),(1,-1),(1,1)])
```


## Demos

For a quick demonstration, run:
```
$ python -m distmesh.demo2d
```

or:
```
$ python -m distmesh.demond
```


## References

The DistMesh algorithm is described in the following two references.
If you use the algorithm in a program or publication, please
acknowledge its authors by adding a reference to the first paper
below.

- P.-O. Persson, G. Strang, **A Simple Mesh Generator in MATLAB**.
  *SIAM Review*, Volume 46 (2), pp. 329-345, June 2004 [PDF](http://persson.berkeley.edu/distmesh/persson04mesh.pdf)

- P.-O. Persson, **Mesh Generation for Implicit Geometries**.
  Ph.D. thesis, *Department of Mathematics, MIT*, Dec 2004 [PDF](http://persson.berkeley.edu/thesis/persson-thesis-color.pdf)
