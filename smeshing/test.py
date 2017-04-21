# encoding: utf-8
"""Distmesh 2D examples."""

#-----------------------------------------------------------------------------
#  Copyright (C) 2004-2012 Per-Olof Persson
#  Copyright (C) 2012 Bradley Froehle

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import numpy as np

# Local imports.
from main import distmesh2d

from file_io import read_data, write_data

#-----------------------------------------------------------------------------
# Demo functions
#-----------------------------------------------------------------------------


def is_same_point(p1, p2):
    tol = 1.0e-12
    if abs(p1[0] - p2[0]) > tol:
        return False
    if abs(p1[1] - p2[1]) > tol:
        return False
    return True


def get_triangle_midpoints(xs, ys, ts):
    ms = []
    for (t1, t2, t3) in ts:
        mx = xs[t1] + xs[t2] + xs[t3]
        my = ys[t1] + ys[t2] + ys[t3]
        ms.append((mx, my))
    return ms


def matches_with_reference(ps, ts, file_name):
    """
    We compare triangle midpoints.
    """
    xs = []
    ys = []
    for p in ps:
        xs.append(p[0])
        ys.append(p[1])

    ms = get_triangle_midpoints(xs, ys, ts)

    xs_ref, ys_ref, ts_ref = read_data(file_name)
    for i, t in enumerate(ts):
        assert t[0] == ts_ref[i][0]
        assert t[1] == ts_ref[i][1]
        assert t[2] == ts_ref[i][2]

    for i in range(len(xs)):
        assert abs((xs[i] - xs_ref[i])/xs[i]) < 1.0e-5
        assert abs((ys[i] - ys_ref[i])/ys[i]) < 1.0e-5


def huniform(p):
    """Implements the trivial uniform mesh size function h=1."""
    return np.ones(p.shape[0])


def polygon(file_name, benchmark=False):

    pv = []
    with open(file_name, 'r') as f:
        for line in f:
            x = float(line.split()[0])
            y = float(line.split()[1])
            pv.append((x, y))
    pv = np.array(pv)

    if benchmark:
         f = huniform
         h0 = 0.03
         _p, _t = distmesh2d(pv, f, h0, (-1, -1, 2, 1), pv, max_num_iterations=100)
    else:
         f = huniform
         h0 = 0.1
         _p, _t = distmesh2d(pv, f, h0, (-1, -1, 2, 1), pv)
    return _p, _t



generate_tests = False
np.random.seed(1) # Always the same results


def test_polygon():
    p, t = polygon('test/polygon.txt')
    if generate_tests:
        write_data(p, t, 'test/result.txt')
    matches_with_reference(p, t, 'test/result.txt')


def test_bench():
    p, t = polygon('test/polygon.txt', benchmark=True)
    if generate_tests:
        write_data(p, t, 'test/result-bench.txt')
    matches_with_reference(p, t, 'test/result-bench.txt')
