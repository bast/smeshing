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
from distance_functions import huniform, drectangle, dcircle, ddiff, dpoly
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


def uniform_mesh_on_unit_circle():
    """Uniform Mesh on Unit Circle"""
    fd = lambda p: np.sqrt((p**2).sum(1))-1.0
    return distmesh2d(fd, huniform, 0.2, (-1,-1,1,1))

def rectangle_with_circular_hole():
    """Rectangle with circular hole, refined at circle boundary"""
    fd = lambda p: ddiff(drectangle(p,-1,1,-1,1), dcircle(p,0,0,0.5))
    fh = lambda p: 0.05+0.3*dcircle(p,0,0,0.5)
    return distmesh2d(fd, fh, 0.05, (-1,-1,1,1),
                         [(-1,-1),(-1,1),(1,-1),(1,1)])

def polygon():
    """Polygon"""
    pv = np.array([(-0.4,-0.5),(0.4,-0.2),(0.4,-0.7),(1.5,-0.4),(0.9,0.1),
                   (1.6,0.8),(0.5,0.5),(0.2,1.0),(0.1,0.4),(-0.7,0.7),
                   (-0.4,-0.5)])
    fd = lambda p: dpoly(p, pv)
    return distmesh2d(fd, huniform, 0.1, (-1,-1, 2,1), pv)

def ellipse():
    """Ellipse"""
    fd = lambda p: p[:,0]**2/2**2 + p[:,1]**2/1**2 - 1
    return distmesh2d(fd, huniform, 0.2, (-2,-1, 2,1))

def square():
    """Square, with size function point and line sources"""
    fd = lambda p: drectangle(p,0,1,0,1)
    fh = lambda p: np.minimum(np.minimum(
        0.01+0.3*abs(dcircle(p,0,0,0)),
        0.025+0.3*abs(dpoly(p,[(0.3,0.7),(0.7,0.5)]))), 0.15)
    return distmesh2d(fd, fh, 0.01, (0,0,1,1), [(0,0), (1,0), (0,1), (1,1)])

def naca0012_airfoil():
    """NACA0012 airfoil"""
    hlead=0.01; htrail=0.04; hmax=2; circx=2; circr=4
    a=.12/.2*np.array([0.2969,-0.1260,-0.3516,0.2843,-0.1036])
    a0=a[0]; a1=np.hstack((a[5:0:-1], 0.0))

    fd = lambda p: ddiff(
        dcircle(p,circx,0,circr),
        (abs(p[:,1])-np.polyval(a1, p[:,0]))**2-a0**2*p[:,0])
    fh = lambda p: np.minimum(np.minimum(
        hlead+0.3*dcircle(p,0,0,0),
        htrail+0.3*dcircle(p,1,0,0)),hmax)

    fixx = 1.0-htrail*np.cumsum(1.3**np.arange(5))
    fixy = a0*np.sqrt(fixx)+np.polyval(a1, fixx)
    fix = np.vstack((
        np.array([(circx-circr,0),(circx+circr,0),
                  (circx,-circr),(circx,circr),
                  (0,0),(1,0)]),
        np.vstack((fixx, fixy)).T,
        np.vstack((fixx, -fixy)).T))
    box = (circx-circr,-circr, circx+circr,circr)
    h0 = min(hlead, htrail, hmax)
    return distmesh2d(fd, fh, h0, box, fix)


generate_tests = False
np.random.seed(1) # Always the same results


def fstats(p, t):
    print('%d nodes, %d elements, min quality %.2f'
          % (len(p), len(t), simpqual(p,t).min()))


def test_circle():
    p, t = uniform_mesh_on_unit_circle()
    fstats(p,t)
    if generate_tests:
        write_data(p, t, 'test/circle.txt')
    assert matches_with_reference(p, t, 'test/circle.txt')


def test_rectangle():
    p, t = rectangle_with_circular_hole()
    fstats(p, t)
    if generate_tests:
        write_data(p, t, 'test/rectangle.txt')
    assert matches_with_reference(p, t, 'test/rectangle.txt')


def test_polygon():
    p, t = polygon()
    fstats(p, t)
    if generate_tests:
        write_data(p, t, 'test/polygon.txt')
    assert matches_with_reference(p, t, 'test/polygon.txt')


def test_ellipse():
    p, t = ellipse()
    fstats(p, t)
    if generate_tests:
        write_data(p, t, 'test/ellipse.txt')
    assert matches_with_reference(p, t, 'test/ellipse.txt')


def test_square():
    p, t = square()
    fstats(p, t)
    if generate_tests:
        write_data(p, t, 'test/square.txt')
    assert matches_with_reference(p, t, 'test/square.txt')


def test_airfoil():
    p, t = naca0012_airfoil()
    fstats(p, t)
    if generate_tests:
        write_data(p, t, 'test/airfoil.txt')
    assert matches_with_reference(p, t, 'test/airfoil.txt')
