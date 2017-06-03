# encoding: utf-8
"""DistMesh 2D"""

# Copyright (C) 2004-2012 Per-Olof Persson
# Copyright (C) 2012 Bradley Froehle

# Distributed under the terms of the GNU General Public License. You should
# have received a copy of the license along with this program. If not,
# see <http://www.gnu.org/licenses/>.

import scipy.spatial as spspatial
import math
import sys
import random
import time
from .bbox import get_bbox


def density_control(p, L, L0, bars, nfix):
    """
    Density control - remove points that are too close.
    """
    apply_density_control = False
    points_to_remove = []

    for i in range(len(L0)):
        if L0[i] > 2.0 * L[i]:
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
        r = F[i] / L[i]
        Fvec.append((r * barvec[i][0], r * barvec[i][1]))

    # compute resulting force on each point from all adjacent bars
    Ftot = [[0.0, 0.0] for _ in p]
    for k, bar in enumerate(bars):
        i, j = bar
        Ftot[i][0] += Fvec[k][0]
        Ftot[i][1] += Fvec[k][1]
        Ftot[j][0] -= Fvec[k][0]
        Ftot[j][1] -= Fvec[k][1]

    return Ftot


def bring_outside_points_back_to_boundary(p, within_bounds, deps, distance_function):
    for i in range(len(p)):
        if not within_bounds[i]:
            px = [p[i][0] + deps, p[i][1]]
            py = [p[i][0], p[i][1] + deps]
            d0, dx, dy = tuple(distance_function([p[i], px, py]))
            dgradx = (dx - d0) / deps
            dgrady = (dy - d0) / deps
            dgrad2 = dgradx**2.0 + dgrady**2.0
            p[i][0] -= d0 * dgradx / dgrad2
            p[i][1] -= d0 * dgrady / dgrad2
    return p


def movement_below_threshold(p, delta_t, Ftot, dptol, h0, within_bounds):
    s = []
    for i in range(len(p)):
        if within_bounds[i]:
            s.append((delta_t * Ftot[i][0])**2.0 + (delta_t * Ftot[i][1])**2.0)
    return max(s) < (dptol * h0)**2.0


def create_initial_distribution(r0_max, points_polygon, num_points, within_bounds_function, fh):
    """
    Create initial distribution in bounding box (equilateral triangles).
    """
    xmin, xmax, ymin, ymax = get_bbox(points_polygon)

    random.seed(1)

    count = 0

    _p = []
    while True:
        print('number of initial points {0} out of {1}'.format(count, num_points))
        _points = []
        for _ in range(10000):
            x = random.uniform(xmin, xmax)
            y = random.uniform(ymin, ymax)
            _points.append([x, y])
        within_bounds = within_bounds_function(_points)
        _points = [point for i, point in enumerate(_points) if within_bounds[i]]
    
        fh_applied = fh(_points)
        r0 = [1.0 / fh_applied[i]**2.0 for i in range(len(_points))]
    
        for i, point in enumerate(_points):
            if random.uniform(0.0, 1.0) < r0[i] / r0_max:
                _p.append(point)
                count += 1
                if count == num_points:
                    return _p


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
        bar_midpoints.append([_px / 2.0, _py / 2.0])
    hbars = fh(bar_midpoints)

    L0 = []
    l2sum = 0.0
    hbars2sum = 0.0
    for i in range(len(L)):
        l2sum += L[i]**2.0
        hbars2sum += hbars[i]**2.0
    for i in range(len(L)):
        L0.append(hbars[i] * Fscale * math.sqrt(l2sum / hbars2sum))

    return L, L0, barvec


def delaunay(p, within_bounds_function):
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
        triangle_centroids.append((x / 3.0, y / 3.0))

    # Keep interior triangles
    # in contrast to original implementation we do not use a tolerance at boundary
    # to avoid a distance computation
    within_bounds = within_bounds_function(triangle_centroids)
    t = []
    for i, triangle in enumerate(_triangles):
        if within_bounds[i]:
            t.append(triangle)

    bars = form_bars(t)

    return pold, bars, t


def prepend_fix_points(pfix, p):
    if pfix is not None:
        nfix = len(pfix)
        p = [point for point in pfix] + [point for point in p]
    else:
        nfix = 0
    return nfix, p


def distmesh2d(config, pv, fh, distance_function, within_bounds_function, h0, pfix=None):
    """
    distmesh2d: 2-D Mesh Generator using Distance Functions.

    Parameters
    ----------
    pv:        list of polygon coordinate tuples
    fh:        Scaled edge length function h(x,y)
    h0:        Initial edge length
    pfix:      Fixed node positions, shape (nfix, 2)

    Returns
    -------
    p:         Node positions (Nx2)
    t:         Triangle indices (NTx3)
    """

    dptol = 0.001
    ttol = 0.1
    Fscale = 1.2
    delta_t = 0.2
    epsilon = sys.float_info.epsilon
    deps = math.sqrt(epsilon) * h0
    density_control_frequency = 30

    t0 = time.time()
    p = create_initial_distribution(config['r0_max'], pv, config['num_points'], within_bounds_function, fh)
    print('time spent in create_initial_distribution: {0:.2f}'.format(time.time() - t0))

    nfix, p = prepend_fix_points(pfix, p)

    # this shift was chosen so that the first movement is large enough to trigger delaunay
    shift = (100.0 * ttol * h0)**2.0
    pold = [[point[0] + shift, point[1] + shift] for point in p]

    count = 0
    while True:
        print('iteration', count)
        t0_iter = time.time()
        count += 1

        if count > config['max_num_iterations']:
            break

        t0 = time.time()
        pold, bars, t = delaunay(p, within_bounds_function)
        print('time spent in delaunay: {0:.2f}'.format(time.time() - t0))

        t0 = time.time()
        L, L0, barvec = get_bar_lengths(p, bars, fh, Fscale)
        print('time spent in get_bar_lengths: {0:.2f}'.format(time.time() - t0))

        t0 = time.time()
        if count % density_control_frequency == 0:
            apply_density_control, p = density_control(p, L, L0, bars, nfix)
            if apply_density_control:
                pold = [[point[0] + shift, point[1] + shift] for point in p]
                continue
        print('time spent in density_control: {0:.2f}'.format(time.time() - t0))

        t0 = time.time()
        F = compute_forces(L0, L, bars, barvec, p)
        print('time spent in compute_forces: {0:.2f}'.format(time.time() - t0))

        t0 = time.time()
        # set force to zero at fixed points
        for i in range(nfix):
            F[i] = [0.0, 0.0]

        # update node positions
        for i in range(len(p)):
            p[i][0] += delta_t * F[i][0]
            p[i][1] += delta_t * F[i][1]
        print('time spent in forces and update: {0:.2f}'.format(time.time() - t0))

        t0 = time.time()
        within_bounds = within_bounds_function(p)
        print('time spent in within_bounds: {0:.2f}'.format(time.time() - t0))

        t0 = time.time()
        p = bring_outside_points_back_to_boundary(p, within_bounds, deps, distance_function)
        print('time spent in bring_outside_points_back_to_boundary: {0:.2f}'.format(time.time() - t0))

      # for the moment not considered
      # if movement_below_threshold(p, delta_t, F, dptol, h0, within_bounds):
      #     break
        print('time spent in iter: {0:.2f}'.format(time.time() - t0_iter))

    print('num points: {0}, num triangles: {1}'.format(len(p), len(t)))
    return p, t
