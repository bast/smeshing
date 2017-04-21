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


def matches_with_reference(ps, ts, file_name):
    """
    We compare triangle midpoints.
    """
    xs = []
    ys = []
    for p in ps:
        xs.append(p[0])
        ys.append(p[1])

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
    return np.ones(len(p))


def dcircle(p, xc, yc, r):
    """Signed distance to circle centered at xc, yc with radius r."""
    return np.sqrt(((p-np.array([xc,yc]))**2).sum(-1))-r


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
       # f = lambda p: 0.05 + 0.3*dcircle(p,0,0,0.01)
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
