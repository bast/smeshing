.. image:: img/example.jpg


Simple mesh generator in Python (fork of PyDistMesh)
====================================================

.. image:: https://travis-ci.org/bast/smeshing.svg?branch=master
   :target: https://travis-ci.org/bast/smeshing/builds

.. image:: https://coveralls.io/repos/github/bast/smeshing/badge.svg?branch=master
   :target: https://coveralls.io/github/bast/smeshing?branch=master

.. image:: https://img.shields.io/badge/license-%20GPL--v3.0-blue.svg
   :target: https://github.com/bast/smeshing/blob/master/LICENSE


This version of the code is based on
`PyDistMesh <https://github.com/bfroehle/pydistmesh>`__ developed by
`Bradley M. Froehle <https://github.com/bfroehle>`__. This code
generates unstructured triangular and tetrahedral meshes using signed
distance functions. Like
`DistMesh <http://persson.berkeley.edu/distmesh/>`__ and
`PyDistMesh <https://github.com/bfroehle/pydistmesh>`__ upon which this
work is based, this code is distributed under the `GNU
GPL <../master/LICENSE>`__.


Status and roadmap for this fork
--------------------------------

Under heavy development. The plan is to be able to do large grids. We
wish to use the code to prepare grids for simulations on the Norwegian
coast with thousands of islands (polygons).

Prototyping work is currently done in Python. Later we will probably
move to Fortran or C(++) and introduce parallelization but provide a
Python interface.

Currently the code is a bit hardcoded for islands and coastline and
specific resolution functions but later the code will be generalized and
the hardcoded aspects will be abstracted out.


Installing dependencies for development
---------------------------------------

::

    virtualenv venv
    source venv/bin/activate
    pip install -r requirements.txt


Running tests
-------------

::

    py.test -vv smeshing/*.py


Plotting
--------

Example::

    python plot.py data/happy-bear/result.txt example.png


Example installation and run
----------------------------

::

    virtualenv venv
    source venv/bin/activate

    pip install --process-dependency-links git+https://github.com/bast/smeshing.git

    smesh --help


Example config.yml
------------------

.. code-block:: yaml

  # number of grid points
  num_grid_points: 5000

  # number of all boundary and coastline interpolation points
  # these will not be part of the grid points
  # instead of num_interpolation_points you can also provide
  # interpolation_step_length using the same units as the coordinates of your data
  num_interpolation_points: 1000

  # number of iterations
  num_iterations: 100

  # view angle for nearest coastline point computations
  view_angle: 90.0


Installation on `Stallo <https://www.sigma2.no/content/stallo>`__
-----------------------------------------------------------------

::

    module load foss/2016b
    module load CMake/3.7.1-foss-2016b
    module load libffi/3.2.1-foss-2016b

    cd ${SLURM_SUBMIT_DIR}

    virtualenv venv
    source venv/bin/activate

    export CC=gcc
    export CXX=g++
    export FC=gfortran

    pip install --process-dependency-links git+https://github.com/bast/smeshing.git


Pros
----

-  Individual components live in separate libraries.
-  A lot of effort was invested in avoiding quadratic scaling.
-  Optimization is fully relaxed.
-  Delaunay is performed at every step.
-  Good memory profile (hopefully, please report if not).


Known issues
------------

-  Code uses shared-memory parallelization but the load leveling is not
   optimal and the scaling has not been studied in detail.


Restart
-------

It is possible to restart a calculation if you provide
``--restart=/path/to/restart/file``.


Why do we need to provide islands and the boundary separately?
--------------------------------------------------------------

For two reasons:

- We compute view vectors for nearest neighbor polygon points in view. For the boundary
  they point to the "inside". For islands they point to the "outside".
- During the computation we need to figure out whether points are inside or outside of polygons.
  We want grid points to be inside the boundary but outside islands.


Why not using GeoJSON?
----------------------

GeoJSON is a nice and standard format but the choice was to prefer a custom format
for the following reasons:

- Meshing should not be restricted to geospatial data
- Meshing should not be restricted to longitude and
  latitude units of decimal degrees but operate on arbitrary units


References
----------

The DistMesh algorithm is described in the following two references. If
you use the algorithm in a program or publication, please acknowledge
its authors by adding a reference to the first paper below.

-  `P.-O. Persson, G. Strang, A Simple Mesh Generator in MATLAB, SIAM
   Review, Volume 46 (2), pp. 329-345, June
   2004 <http://persson.berkeley.edu/distmesh/persson04mesh.pdf>`__
-  `P.-O. Persson, Mesh Generation for Implicit Geometries, Ph.D.
   thesis, Department of Mathematics, MIT, Dec
   2004 <http://persson.berkeley.edu/thesis/persson-thesis-color.pdf>`__
