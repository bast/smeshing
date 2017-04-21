# encoding: utf-8
"""DistMesh 2D"""

#-----------------------------------------------------------------------------
#  Copyright (C) 2004-2012 Per-Olof Persson
#  Copyright (C) 2012 Bradley Froehle

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

import numpy as np
import scipy.spatial as spspatial
import polygons
import math
import sys

# Local imports
import mlcompat as ml


def density_control(p, count, densityctrlfreq, L, L0, bars, nfix):
    """
    Density control - remove points that are too close.
    """
    apply_density_control = False
    points_to_remove = []
    if count % densityctrlfreq == 0:
        for i in range(len(L0)):
            if L0[i] > 2.0*L[i]:
                apply_density_control = True
                for k in [0, 1]:
                    ip = bars[i][k]
                    if ip > nfix:
                        points_to_remove.append(ip)
    if apply_density_control:
        points_to_remove = list(set(points_to_remove))
        p = np.array([p[ip] for ip in range(len(p)) if ip not in points_to_remove])
    return apply_density_control, p


def compute_forces(L0, L, bars, barvec, p):
    # Bar forces (scalars)
    F = L0 - L

    # we ignore attractive forces and only keep repulsion
    F[F < 0.0] = 0.0

    # compute forces along bars
    Fvec = np.zeros((len(F), 2))
    for i in range(len(F)):
        r = F[i]/L[i]
        Fvec[i][0] = r*barvec[i][0]
        Fvec[i][1] = r*barvec[i][1]

    # compute resulting force on each point from all adjacent bars
    Ftot = np.zeros((len(p), 2))
    for k, bar in enumerate(bars):
        i, j = bar
        Ftot[i][0] += Fvec[k][0]
        Ftot[i][1] += Fvec[k][1]
        Ftot[j][0] -= Fvec[k][0]
        Ftot[j][1] -= Fvec[k][1]

    return Ftot


def bring_outside_points_back_to_boundary(p, contains, deps, polygons_context):
    for i in range(len(p)):
        if not contains[i]:
            px = p[i] + [deps, 0]
            py = p[i] + [0, deps]
            d0 = polygons.get_distances(polygons_context, [p[i]])[0]
            dx = polygons.get_distances(polygons_context, [px])[0]
            dy = polygons.get_distances(polygons_context, [py])[0]
            dgradx = (dx - d0)/deps
            dgrady = (dy - d0)/deps
            dgrad2 = dgradx**2.0 + dgrady**2.0
            p[i][0] -= d0*dgradx/dgrad2
            p[i][1] -= d0*dgrady/dgrad2
    return p


def movement_below_threshold(p, delta_t, Ftot, dptol, h0, contains):
    s = []
    for i in range(len(p)):
        if contains[i]:
            tmp = delta_t*Ftot[i]**2.0
            s.append(tmp[0] + tmp[1])
    return max(s) < (dptol*h0)**2.0


def create_initial_distribution(bbox, h0):
    """
    Create initial distribution in bounding box (equilateral triangles).
    """
    xmin, ymin, xmax, ymax = bbox

    xs = [xmin]
    while xs[-1] <= xmax:
        xs.append(xs[-1] + h0)
    ys = [ymin]
    while ys[-1] <= ymax:
        ys.append(ys[-1] + h0*math.sqrt(3.0)/2.0)

    points = []
    for x in xs:
        for row, y in enumerate(ys):
            if row%2 != 0:
                # shift every second row to the right
                points.append([x + h0/2.0, y])
            else:
                points.append([x, y])
    return points


def form_bars(triangles):
    bars = []
    for triangle in triangles:
        bars.append((triangle[0], triangle[1]))
        bars.append((triangle[1], triangle[2]))
        bars.append((triangle[2], triangle[0]))

    _bars = []
    for bar in bars:
        if bar[1] < bar[0]:
            _bars.append((bar[1], bar[0]))
        else:
            _bars.append(bar)
    _bars = set(_bars)
    _bars = list(_bars)
    _bars = sorted(_bars)
    return _bars


def get_bar_lengths(p, bars, fh, Fscale):
    barvec = p[bars[:,0]] - p[bars[:,1]]
    _barvec = []
    for bar in bars:
        vx = p[bar[0]][0] - p[bar[1]][0]
        vy = p[bar[0]][1] - p[bar[1]][1]
        _barvec.append([vx, vy])
    _L = []
    for bar in _barvec:
        _L.append(math.sqrt(bar[0]**2.0 + bar[1]**2.0))
    L = np.array(_L)

    bar_midpoints = []
    for bar in bars:
        _px = p[bar[0]][0] + p[bar[1]][0]
        _py = p[bar[0]][1] + p[bar[1]][1]
        bar_midpoints.append([_px/2.0, _py/2.0])
    bar_midpoints = np.array(bar_midpoints)
    hbars = fh(bar_midpoints)

    _L0 = []
    l2sum = 0.0
    hbars2sum = 0.0
    for i in range(len(_L)):
        l2sum += _L[i]**2.0
        hbars2sum += hbars[i]**2.0
    for i in range(len(_L)):
        _L0.append(hbars[i]*Fscale*math.sqrt(l2sum/hbars2sum))
    L0 = np.array(_L0)

    return L, L0, barvec


def large_movement(p, pold, ttol, h0):
    _temp = 0.0
    for i in range(len(p)):
        _d = (p[i][0] - pold[i][0])**2.0 + (p[i][1] - pold[i][1])**2.0
        _temp = max(_temp, _d)
    return _temp > (ttol*h0)**2.0


def delaunay(p, polygons_context):
    """
    Retriangulation by the Delaunay algorithm.
    """
    pold = p.copy()                          # Save current positions
    _triangles = spspatial.Delaunay(p).vertices       # List of triangles
    pmid = p[_triangles].sum(1)/3                     # Compute centroids

    # Keep interior triangles
    # in contrast to original implementation we do not use a tolerance at boundary
    # to avoid a distance computation
    contains = polygons.contains_points(polygons_context, pmid)
    t = []
    for i, triangle in enumerate(_triangles):
        if contains[i]:
            t.append(triangle)
    t = np.array(t)

    bars = form_bars(t)
    bars = np.array(bars)

    return pold, bars, t


def distmesh2d(pv, fh, h0, bbox, pfix=None, max_num_iterations=None):
    """
    distmesh2d: 2-D Mesh Generator using Distance Functions.

    Parameters
    ----------
    pv:        list of polygon coordinate tuples
    fh:        Scaled edge length function h(x,y)
    h0:        Initial edge length
    bbox:      Bounding box, (xmin, ymin, xmax, ymax)
    pfix:      Fixed node positions, shape (nfix, 2)

    Returns
    -------
    p:         Node positions (Nx2)
    t:         Triangle indices (NTx3)
    """

    polygons_context = polygons.new_context()

    polygons.add_polygon(polygons_context, pv)

    dptol = 0.001
    ttol = .1
    Fscale = 1.2
    delta_t = 0.2
    geps = .001*h0
    smallest_representable_double_float = np.finfo(np.double).eps
    deps = math.sqrt(smallest_representable_double_float)*h0
    densityctrlfreq = 30

    if pfix is not None:
        pfix = np.array(pfix, dtype='d')

    _points = np.array(create_initial_distribution(bbox, h0))

    # 2. Remove points outside the region, apply the rejection method

    contains = polygons.contains_points(polygons_context, _points)
    p = []
    for i, point in enumerate(_points):
        if contains[i]:
            p.append(point)
    p = np.array(p)

    # Probability to keep point
    r0 = 1/fh(p)**2.0
    r0_max = max(r0)
    randoms = np.random.random(p.shape[0])
    _p = []
    for i, point in enumerate(p):
        if randoms[i] < r0[i]/r0_max:
            _p.append(point)
    p = np.array(_p)

    if pfix is not None:
        p = ml.setdiff_rows(p, pfix)                 # Remove duplicated nodes

        _pfix = pfix.tolist()
        _pfix = [tuple(x) for x in _pfix]
        _pfix = set(_pfix)
        _pfix = list(_pfix)
        pfix = np.array(_pfix)

        nfix = pfix.shape[0]
        p = np.vstack((pfix, p))                     # Prepend fix points
    else:
        nfix = 0

    shift = 100.0*ttol*h0**2.0  # this shift is so that the first movement is large enough to trigger delaunay
    pold = np.array([[point[0] + shift, point[1] + shift] for point in p])

    count = 0
    while True:
        count += 1

        if max_num_iterations is not None:
            if count > max_num_iterations:
                break

        if large_movement(p, pold, ttol, h0):
            pold, bars, t = delaunay(p, polygons_context)

        L, L0, barvec = get_bar_lengths(p, bars, fh, Fscale)

        apply_density_control, p = density_control(p, count, densityctrlfreq, L, L0, bars, nfix)
        if apply_density_control:
            pold = np.array([[point[0] + shift, point[1] + shift] for point in p])
            continue

        F = compute_forces(L0, L, bars, barvec, p)

        # set force to zero at fixed points
        F[:nfix] = 0.0

        # update node positions
        p += delta_t*F

        contains = polygons.contains_points(polygons_context, p)

        p = bring_outside_points_back_to_boundary(p, contains, deps, polygons_context)

        if movement_below_threshold(p, delta_t, F, dptol, h0, contains):
            break

    polygons.free_context(polygons_context)

    return p, t
