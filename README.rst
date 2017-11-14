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

    python smeshing/plot.py data/fiction/result.txt example.png


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

  # the smaller this number, the better starting grid
  seeding_speed: 20.0
  # number of grid points
  num_grid_points: 5000
  # number of boundary points
  num_boundary_points: 1000
  # number of iterations
  num_iterations: 100
  # magic factor, needs documentation
  scale_factor: 0.995792
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
-  Optimization is fully relaxed - no interpolation is done.
-  Delaunay is performed at every step.
-  Good memory profile.


Known issues
------------

-  Code uses shared-memory parallelization but the load leveling is not
   optimal and the scaling on Stallo has not been studied.
-  Sometimes a step generates crazy meshes. Possibly there is a
   numerical instability for small forces or division by very small
   number.
-  One needs to set a ``seeding_speed``. This controls in how many steps
   the initial point distribution is set up. The higher the number, the
   fewer steps, and the more uniform distribution. The more steps, the
   more time it takes, but also the more it will reflect the resolution
   function. Currently it needs some experimentation. Later we need
   something more black box.
-  Coastline resolution is currently hard-coded to be inversely
   proportional to the distance to nearest neighbor in view divided by
   6. In future versions this will be read from input.
-  Currently no lower and upper bounds on resolution can be set.
-  There is no stop criterion, it will run as many iterations as you ask
   it to.


Restart
-------

It is possible to restart a calculation if you provide
``--restart=/path/to/restart/file``.


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
