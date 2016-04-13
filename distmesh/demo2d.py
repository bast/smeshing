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
import matplotlib.pyplot as plt

# Local imports.
import distmesh as dm

#-----------------------------------------------------------------------------
# Demo functions
#-----------------------------------------------------------------------------


def save_p_t(p, t, file_name):
    with open(file_name, 'w') as f:
        f.write('{}\n'.format(len(p)))
        for px, py in p:
            f.write('{} {}\n'.format(px, py))
        f.write('{}\n'.format(len(t)))
        for t1, t2, t3 in t:
            f.write('{} {} {}\n'.format(t1, t2, t3))


def matches_with_reference(p, t, file_name):

    tol = 1.0e-12

    with open(file_name, 'r') as f:
        lines = f.readlines()

    offset = 0

    p_len_ref = int(lines[offset])
    assert len(p) == p_len_ref
    offset += 1

    for i, (px, py) in enumerate(p):
        px_ref = float(lines[i + offset].split()[0])
        py_ref = float(lines[i + offset].split()[1])
        assert abs(px - px_ref) < tol
        assert abs(py - py_ref) < tol
    offset += len(p)

    t_len_ref = int(lines[offset])
    assert len(t) == t_len_ref
    offset += 1

    for i, (t1, t2, t3) in enumerate(t):
        t1_ref = int(lines[i + offset].split()[0])
        t2_ref = int(lines[i + offset].split()[1])
        t3_ref = int(lines[i + offset].split()[2])
        assert t1 == t1_ref
        assert t2 == t2_ref
        assert t3 == t3_ref
    offset += len(t)

    return True


def uniform_mesh_on_unit_circle():
    """Uniform Mesh on Unit Circle"""
    fd = lambda p: np.sqrt((p**2).sum(1))-1.0
    return dm.distmesh2d(fd, dm.huniform, 0.2, (-1,-1,1,1))

def rectangle_with_circular_hole():
    """Rectangle with circular hole, refined at circle boundary"""
    fd = lambda p: dm.ddiff(dm.drectangle(p,-1,1,-1,1), dm.dcircle(p,0,0,0.5))
    fh = lambda p: 0.05+0.3*dm.dcircle(p,0,0,0.5)
    return dm.distmesh2d(fd, fh, 0.05, (-1,-1,1,1),
                         [(-1,-1),(-1,1),(1,-1),(1,1)])

def polygon():
    """Polygon"""
    pv = np.array([(-0.4,-0.5),(0.4,-0.2),(0.4,-0.7),(1.5,-0.4),(0.9,0.1),
                   (1.6,0.8),(0.5,0.5),(0.2,1.0),(0.1,0.4),(-0.7,0.7),
                   (-0.4,-0.5)])
    fd = lambda p: dm.dpoly(p, pv)
    return dm.distmesh2d(fd, dm.huniform, 0.1, (-1,-1, 2,1), pv)

def ellipse():
    """Ellipse"""
    fd = lambda p: p[:,0]**2/2**2 + p[:,1]**2/1**2 - 1
    return dm.distmesh2d(fd, dm.huniform, 0.2, (-2,-1, 2,1))

def square():
    """Square, with size function point and line sources"""
    fd = lambda p: dm.drectangle(p,0,1,0,1)
    fh = lambda p: np.minimum(np.minimum(
        0.01+0.3*abs(dm.dcircle(p,0,0,0)),
        0.025+0.3*abs(dm.dpoly(p,[(0.3,0.7),(0.7,0.5)]))), 0.15)
    return dm.distmesh2d(fd, fh, 0.01, (0,0,1,1), [(0,0), (1,0), (0,1), (1,1)])

def naca0012_airfoil():
    """NACA0012 airfoil"""
    hlead=0.01; htrail=0.04; hmax=2; circx=2; circr=4
    a=.12/.2*np.array([0.2969,-0.1260,-0.3516,0.2843,-0.1036])
    a0=a[0]; a1=np.hstack((a[5:0:-1], 0.0))

    fd = lambda p: dm.ddiff(
        dm.dcircle(p,circx,0,circr),
        (abs(p[:,1])-np.polyval(a1, p[:,0]))**2-a0**2*p[:,0])
    fh = lambda p: np.minimum(np.minimum(
        hlead+0.3*dm.dcircle(p,0,0,0),
        htrail+0.3*dm.dcircle(p,1,0,0)),hmax)

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
    return dm.distmesh2d(fd, fh, h0, box, fix)

def meshdemo2d():
    """Run all Distmesh 2D examples."""

    generate_tests = False

    plt.ion()
    np.random.seed(1) # Always the same results

    def fstats(p, t):
        print('%d nodes, %d elements, min quality %.2f'
              % (len(p), len(t), dm.simpqual(p,t).min()))

    print('Uniform Mesh on Unit Circle')
    p, t = uniform_mesh_on_unit_circle()
    fstats(p,t)
    if generate_tests:
        save_p_t(p, t, 'test/circle.txt')
    assert matches_with_reference(p, t, 'test/circle.txt')
    print('')

    print('Rectangle with circular hole, refined at circle boundary')
    p, t = rectangle_with_circular_hole()
    fstats(p, t)
    if generate_tests:
        save_p_t(p, t, 'test/rectangle.txt')
    assert matches_with_reference(p, t, 'test/rectangle.txt')
    print('')

    print('Polygon')
    p, t = polygon()
    fstats(p, t)
    if generate_tests:
        save_p_t(p, t, 'test/polygon.txt')
    assert matches_with_reference(p, t, 'test/polygon.txt')
    print('')

    print('Ellipse')
    p, t = ellipse()
    fstats(p, t)
    if generate_tests:
        save_p_t(p, t, 'test/ellipse.txt')
    assert matches_with_reference(p, t, 'test/ellipse.txt')
    print('')

    print('Square, with size function point and line sources')
    p, t = square()
    fstats(p, t)
    if generate_tests:
        save_p_t(p, t, 'test/square.txt')
    assert matches_with_reference(p, t, 'test/square.txt')
    print('')

    print('NACA0012 airfoil')
    p, t = naca0012_airfoil()
    fstats(p, t)
    if generate_tests:
        save_p_t(p, t, 'test/airfoil.txt')
    assert matches_with_reference(p, t, 'test/airfoil.txt')

if __name__ == "__main__":
    meshdemo2d()
