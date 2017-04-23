# encoding: utf-8
"""DistMesh 2D"""

# Copyright (C) 2004-2012 Per-Olof Persson
# Copyright (C) 2012 Bradley Froehle

# Distributed under the terms of the GNU General Public License. You should
# have received a copy of the license along with this program. If not,
# see <http://www.gnu.org/licenses/>.

import scipy.spatial as spspatial
import polygons
import math
import sys
import random


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
        p = [p[ip] for ip in range(len(p)) if ip not in points_to_remove]
    return apply_density_control, p


def compute_forces(L0, L, bars, barvec, p):
    # Bar forces (scalars)
    F = [L0[i] - L[i] for i in range(len(L))]

    # we ignore attractive forces and only keep repulsion
    for i in range(len(F)):
        if F[i] < 0.0:
            F[i] = 0.0

    # compute forces along bars
    Fvec = []
    for i in range(len(F)):
        r = F[i]/L[i]
        Fvec.append((r*barvec[i][0], r*barvec[i][1]))

    # compute resulting force on each point from all adjacent bars
    Ftot = [[0.0, 0.0] for _ in p]
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
            px = [p[i][0] + deps, p[i][1]]
            py = [p[i][0], p[i][1] + deps]
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
            s.append(delta_t*Ftot[i][0]**2.0 + delta_t*Ftot[i][1]**2.0)
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
            if row % 2 != 0:
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
    barvec = []
    for bar in bars:
        vx = p[bar[0]][0] - p[bar[1]][0]
        vy = p[bar[0]][1] - p[bar[1]][1]
        barvec.append([vx, vy])
    L = []
    for bar in barvec:
        L.append(math.sqrt(bar[0]**2.0 + bar[1]**2.0))

    bar_midpoints = []
    for bar in bars:
        _px = p[bar[0]][0] + p[bar[1]][0]
        _py = p[bar[0]][1] + p[bar[1]][1]
        bar_midpoints.append([_px/2.0, _py/2.0])
    hbars = [fh(x, y) for (x, y) in bar_midpoints]

    L0 = []
    l2sum = 0.0
    hbars2sum = 0.0
    for i in range(len(L)):
        l2sum += L[i]**2.0
        hbars2sum += hbars[i]**2.0
    for i in range(len(L)):
        L0.append(hbars[i]*Fscale*math.sqrt(l2sum/hbars2sum))

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
    # save current positions
    pold = []
    for _p in p:
        pold.append([_p[0], _p[1]])
    _triangles = spspatial.Delaunay(p).vertices       # List of triangles

    triangle_centroids = []
    for triangle in _triangles:
        x = 0.0
        y = 0.0
        for ip in triangle:
            x += p[ip][0]
            y += p[ip][1]
        triangle_centroids.append((x/3.0, y/3.0))

    # Keep interior triangles
    # in contrast to original implementation we do not use a tolerance at boundary
    # to avoid a distance computation
    contains = polygons.contains_points(polygons_context, triangle_centroids)
    t = []
    for i, triangle in enumerate(_triangles):
        if contains[i]:
            t.append(triangle)

    bars = form_bars(t)

    return pold, bars, t


def remove_points_outside_region(polygons_context, points):
    contains = polygons.contains_points(polygons_context, points)
    return [point for i, point in enumerate(points) if contains[i]]


def apply_rejection_method(fh, p):
    r0 = [1.0/fh(x, y)**2.0 for (x, y) in p]
    r0_max = max(r0)

    random.seed(1)
    _p = []
    for i, point in enumerate(p):
        if random.uniform(0.0, 1.0) < r0[i]/r0_max:
            _p.append(point)
    return _p


def prepend_fix_points(pfix, p):
    if pfix is not None:
        nfix = len(pfix)
        p = [point for point in pfix] + [point for point in p]
    else:
        nfix = 0
    return nfix, p


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
    geps = 0.001*h0
    epsilon = sys.float_info.epsilon
    deps = math.sqrt(epsilon)*h0
    densityctrlfreq = 30

    _points = create_initial_distribution(bbox, h0)

    p = remove_points_outside_region(polygons_context, _points)

    p = apply_rejection_method(fh, p)

    nfix, p = prepend_fix_points(pfix, p)

    shift = 100.0*ttol*h0**2.0  # this shift is so that the first movement is large enough to trigger delaunay
    pold = [[point[0] + shift, point[1] + shift] for point in p]

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
            pold = [[point[0] + shift, point[1] + shift] for point in p]
            continue

        F = compute_forces(L0, L, bars, barvec, p)

        # set force to zero at fixed points
        for i in range(nfix):
            F[i] = [0.0, 0.0]

        # update node positions
        for i in range(len(p)):
            p[i][0] += delta_t*F[i][0]
            p[i][1] += delta_t*F[i][1]

        contains = polygons.contains_points(polygons_context, p)

        p = bring_outside_points_back_to_boundary(p, contains, deps, polygons_context)

        if movement_below_threshold(p, delta_t, F, dptol, h0, contains):
            break

    polygons.free_context(polygons_context)

    print('num points: {0}, num triangles: {1}'.format(len(p), len(t)))
    return p, t
