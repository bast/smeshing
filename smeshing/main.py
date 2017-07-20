# encoding: utf-8
"""DistMesh 2D"""

# Copyright (C) 2004-2012 Per-Olof Persson
# Copyright (C) 2012 Bradley Froehle

# Distributed under the terms of the GNU General Public License. You should
# have received a copy of the license along with this program. If not,
# see <http://www.gnu.org/licenses/>.

import math
import sys
import random
import yaml
import time
import polygons
import flanders
import delaunay
from .bbox import get_bbox
from .clockwise import edges_sum
from .file_io import read_data


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
        for _ in range(min(num_points, 1000000)):
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


def solve_delaunay(p, within_bounds_function):
    """
    Retriangulation by the Delaunay algorithm.
    """
    # save current positions
    pold = []
    for _p in p:
        pold.append([_p[0], _p[1]])

    _triangles = delaunay.solve(p)

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


def distmesh2d(config,
               pv,
               fh,
               distance_function,
               within_bounds_function,
               h0,
               pfix=None,
               restart_file_name=None):
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

    if restart_file_name is None:
        t0 = time.time()
        p = create_initial_distribution(config['r0_max'], pv, config['num_points'], within_bounds_function, fh)
        print('time spent in create_initial_distribution: {0:.2f}'.format(time.time() - t0))
        nfix, p = prepend_fix_points(pfix, p)
    else:
        nfix = len(pfix)
        p, _ = read_data(restart_file_name)

    # remove points which are on top of other points
    p_tuples = [(point[0], point[1]) for point in p]
    p = [[point[0], point[1]] for point in set(p_tuples)]

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
        pold, bars, t = solve_delaunay(p, within_bounds_function)
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

        print('time spent in iter: {0:.2f}'.format(time.time() - t0_iter))

    print('num points: {0}, num triangles: {1}'.format(len(p), len(t)))
    return p, t


def run(boundary_file_name,
        island_file_names,
        config_file_name,
        restart_file_name=None):

    with open(config_file_name, 'r') as f:
        try:
            config = yaml.load(f, yaml.SafeLoader)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit(-1)

    plot_nearest_in_view = False
    if plot_nearest_in_view:
        import matplotlib.pyplot as plt

    boundary_points = read_points(boundary_file_name)
    view_vectors = compute_view_vectors(boundary_points, scale=-1.0)
    all_points = boundary_points
    for island_file in island_file_names:
        islands_points = read_points(island_file)
        view_vectors += compute_view_vectors(islands_points, scale=1.0)
        all_points += islands_points
    num_points = len(all_points)
    flanders_context = flanders.new_context(num_points, all_points)
    angles_deg = [90.0 for _ in range(num_points)]
    flanders_indices = flanders.search_neighbor(context=flanders_context,
                                                ref_indices=list(range(num_points)),
                                                view_vectors=view_vectors,
                                                angles_deg=angles_deg)

    nearest_distance_at_coastline_point = []
    for i in range(len(all_points)):
        nearest_distance_at_coastline_point.append(get_distance(all_points[i], all_points[flanders_indices[i]]) / 6.0)  # FIXME 6.0 is hardcoded

    flanders.free_context(flanders_context)

    all_polygons_context = polygons.new_context()
    boundary_context = polygons.new_context()
    islands_context = polygons.new_context()

    boundary_points = read_points(boundary_file_name)
    polygons.add_polygon(all_polygons_context, boundary_points, nearest_distance_at_coastline_point[0:len(boundary_points)])
    polygons.add_polygon(boundary_context, boundary_points, nearest_distance_at_coastline_point[0:len(boundary_points)])
    counter = len(boundary_points)
    all_points = boundary_points
    if plot_nearest_in_view:
        for i in range(len(boundary_points) - 1):
            plt.plot([boundary_points[i][0], boundary_points[i + 1][0]],
                     [boundary_points[i][1], boundary_points[i + 1][1]],
                     'r-')

    for island_file in island_file_names:
        islands_points = read_points(island_file)
        polygons.add_polygon(all_polygons_context, islands_points, nearest_distance_at_coastline_point[counter:counter + len(islands_points)])
        polygons.add_polygon(islands_context, islands_points, nearest_distance_at_coastline_point[counter:counter + len(islands_points)])
        counter += len(islands_points)
        all_points += islands_points
        if plot_nearest_in_view:
            for i in range(len(islands_points) - 1):
                plt.plot([islands_points[i][0], islands_points[i + 1][0]],
                         [islands_points[i][1], islands_points[i + 1][1]],
                         'b-')

    def distance_function(points):
        return polygons.get_distances_edge(all_polygons_context, points)

    def within_bounds_function(points):
        within_boundary = polygons.contains_points(boundary_context, points)
        within_islands = polygons.contains_points(islands_context, points)  # FIXME we don't need to check all points
        within_bounds = []
        for i, point in enumerate(points):
            if not within_boundary[i]:
                within_bounds.append(False)
            else:
                if within_islands[i]:
                    within_bounds.append(False)
                else:
                    within_bounds.append(True)
        return within_bounds

    xmin, xmax, ymin, ymax = get_bbox(boundary_points)
    h0 = (xmax - xmin) / 500.0

    def _r(points):
        # FIXME slope is hardcoded
        return polygons.get_distances_vertex_weighted(all_polygons_context, points, [0.995792] * len(points))
    h_function = _r

    if plot_nearest_in_view:
        for i in range(len(all_points)):
            if flanders_indices[i] > -1:
                plt.plot([all_points[i][0], all_points[flanders_indices[i]][0]],
                         [all_points[i][1], all_points[flanders_indices[i]][1]],
                         'k-')
            else:
                plt.plot([all_points[i][0], 0.0],
                         [all_points[i][1], 0.0],
                         'g-')
                print('-1 distance found for x={0} y={1}'.format(all_points[i][0], all_points[i][1]))
        plt.savefig('foo.png')

    points, triangles = distmesh2d(config,
                                   all_points,
                                   h_function,
                                   distance_function,
                                   within_bounds_function,
                                   h0,
                                   all_points,
                                   restart_file_name)

    polygons.free_context(all_polygons_context)
    polygons.free_context(boundary_context)
    polygons.free_context(islands_context)

    return points, triangles


def read_points(file_name):
    points = []
    with open(file_name, 'r') as f:
        for line in f:
            x = float(line.split()[0])
            y = float(line.split()[1])
            points.append([x, y])
    return points


def compute_view_vectors(points, scale):
    """
    If scale is negative, then view vectors are towards inside.
    """

    # we figure out whether polygon is clockwise or anticlockwise
    if edges_sum(points) < 0.0:
        s = +1.0 * scale
    else:
        s = -1.0 * scale

    # we remove the last point since it repeats the first
    # we assume clock-wise closed loops
    _points = [points[i] for i in range(len(points) - 1)]

    if len(_points) > 2:
        vectors = get_normal_vectors(_points, s)
    else:
        # this is a "linear" island consisting of two points
        vectors = []
        vector = (_points[1][0] - _points[0][0], _points[1][1] - _points[0][1])
        vectors.append(normalize(vector, -s))
        vectors.append(normalize(vector, s))

    # we add the removed point again
    return vectors + [vectors[0]]


def get_normal_vectors(points, s):
    num_points = len(points)
    vectors = []
    for i in range(num_points):
        i_before = i - 1
        i_after = (i + 1) % num_points
        vector = (points[i_after][1] - points[i_before][1], -(points[i_after][0] - points[i_before][0]))
        vector = normalize(vector, s)
        vectors.append(vector)
    return vectors


def normalize(vector, s):
    norm = math.sqrt(vector[0]**2.0 + vector[1]**2.0)
    return (s * vector[0] / norm, s * vector[1] / norm)


def get_distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2.0 + (p2[1] - p1[1])**2.0)
