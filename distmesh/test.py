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
from _distmesh2d import distmesh2d
from distance_functions import huniform, dcircle, dpoly
from utils import simpqual

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
    ms_ref = get_triangle_midpoints(xs_ref, ys_ref, ts_ref)

    for m in ms:
        # we search m in m_ref, as soon as we find it, we pop it
        for i, m_ref in enumerate(ms_ref):
            if is_same_point(m, m_ref):
                ms_ref.pop(i)
                break

    # the points match if ms_ref is emptied
    return len(ms_ref) == 0


def polygon():
    """Polygon"""
    pv = np.array([(-0.4,-0.5),(0.4,-0.2),(0.4,-0.7),(1.5,-0.4),(0.9,0.1),
                   (1.6,0.8),(0.5,0.5),(0.2,1.0),(0.1,0.4),(-0.7,0.7),
                   (-0.4,-0.5)])
    fd = lambda p: dpoly(p, pv)
    fh = lambda p: 0.05 + 0.3*dcircle(p,0,0,0.01)
    return distmesh2d(fd, fh, 0.1, (-1,-1, 2,1), pv)
  # return distmesh2d(fd, huniform, 0.1, (-1,-1, 2,1), pv)



generate_tests = False
np.random.seed(1) # Always the same results


def fstats(p, t):
    print('%d nodes, %d elements, min quality %.2f'
          % (len(p), len(t), simpqual(p,t).min()))


def test_polygon():
    p, t = polygon()
    fstats(p, t)
    if generate_tests:
        write_data(p, t, 'test/result.txt')
    assert matches_with_reference(p, t, 'test/result.txt')
