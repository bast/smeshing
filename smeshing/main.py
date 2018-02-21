# encoding: utf-8
"""DistMesh 2D"""

# Copyright (C) 2004-2012 Per-Olof Persson
# Copyright (C) 2012 Bradley Froehle

# Distributed under the terms of the GNU General Public License. You should
# have received a copy of the license along with this program. If not,
# see <http://www.gnu.org/licenses/>.

import math
import sys
import os
import random
import yaml
import time
import ntpath
import numpy as np
from scipy.interpolate import RectBivariateSpline, SmoothBivariateSpline

import polygons

from .bbox import get_bbox
from .file_io import read_data
from .version import __version__
from .triangulate import solve_delaunay


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


def minimize_over_resolution_fields(distance_fields,
                                    distance_field_bounds,
                                    distances,
                                    points):
    '''
    Function receives distances and distance fields.
    It will return an array of distances which minimize
    distances and interpolated resolution fields.
    '''
    _distances = []
    for i, (x, y) in enumerate(points):
        d = distances[i]
        for j, distance_field in enumerate(distance_fields):
            d_interpolated = distance_field(x, y)[0][0]
            min_bound, max_bound = distance_field_bounds[j]
            d_interpolated = max(d_interpolated, min_bound)
            d_interpolated = min(d_interpolated, max_bound)
        _distances.append(min(d, d_interpolated))
    return _distances


def get_random_points(num_points, xmin, xmax, ymin, ymax):
    points = []
    for _ in range(num_points):
        x = random.uniform(xmin, xmax)
        y = random.uniform(ymin, ymax)
        points.append([x, y])
    return points


def filter_points(points, resolutions, max_resolution, damping_factor, max_num_filtered_points):
    filtered_points = []
    for point, resolution in zip(points, resolutions):
        if random.uniform(0.0, 1.0) < damping_factor * resolution / max_resolution:
            filtered_points.append(point)
    num_filtered_points = len(filtered_points)
    if num_filtered_points > max_num_filtered_points:
        return filtered_points[:max_num_filtered_points]
    else:
        return filtered_points


def create_initial_distribution(polygon_points, num_grid_points, within_bounds_function, resolution_function):
    """
    Create initial distribution in bounding box (equilateral triangles).
    """
    xmin, xmax, ymin, ymax = get_bbox(polygon_points)

    random.seed(1)

    count = 0

    grid_points = []

    # consider taking the diagonal as starting point
    damping_factor = max(xmax - xmin, ymax - ymin)
    damping_factor_calibrated = False

    print('\ngenerating starting grid:')
    _num_steps = 0
    while True:
        # get a larger batch of points
        points = get_random_points(min(num_grid_points, 1000000), xmin, xmax, ymin, ymax)

        # only keep those which are within bounds
        within_bounds = within_bounds_function(points)
        points = [point for i, point in enumerate(points) if within_bounds[i]]

        resolution_function_applied = resolution_function(points)
        resolutions = [1.0 / resolution_function_applied[i]**2.0 for i in range(len(points))]
        max_resolution = max(resolutions)

        max_num_filtered_points = num_grid_points - count

        if not damping_factor_calibrated:
            # here we try to find a good enough damping factor
            target_percentage = 1.0
            while True:
                _filtered_points = filter_points(points, resolutions, max_resolution, damping_factor, max_num_filtered_points)
                _num_filtered_points = len(_filtered_points)
                _percentage = 100.0 * _num_filtered_points / num_grid_points
                if _percentage < target_percentage:
                    damping_factor_calibrated = True
                    break
                else:
                    damping_factor /= 2.0

        filtered_points = filter_points(points, resolutions, max_resolution, damping_factor, max_num_filtered_points)
        count += len(filtered_points)
        grid_points += filtered_points

        print('    generated {0} points out of {1}'.format(count, num_grid_points))
        _num_steps += 1
        if count == num_grid_points:
            break
    print('\ngenerated initial distribution in {0} steps'.format(_num_steps))
    return grid_points


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


def distmesh2d(config,
               pv,
               fh,
               distance_function,
               within_bounds_function,
               h0,
               boundary_points,
               restart_file_name=None):
    """
    distmesh2d: 2-D Mesh Generator using Distance Functions.

    Parameters
    ----------
    pv:        list of polygon coordinate tuples
    fh:        Scaled edge length function h(x,y)
    h0:        Initial edge length

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

    print_timing = False
    if 'print_timing' in config:
        if config['print_timing']:
            print_timing = True

    if restart_file_name is None:
        t0 = time.time()
        p = create_initial_distribution(pv, config['num_grid_points'], within_bounds_function, fh)
        print('time spent in creating initial distribution: {0:.2f}'.format(time.time() - t0))
    else:
        p, _ = read_data(restart_file_name)

    print('\nnumber of grid points: {0}'.format(len(p)))

    count = 0
    print('\nrunning iterations:'.format(count))
    while True:
        t0_iter = time.time()

        t0 = time.time()
        bars, t = solve_delaunay(p, within_bounds_function)
        if print_timing:
            print('time spent in delaunay: {0:.2f}'.format(time.time() - t0))

        t0 = time.time()
        L, L0, barvec = get_bar_lengths(p, bars, fh, Fscale)
        if print_timing:
            print('time spent in get_bar_lengths: {0:.2f}'.format(time.time() - t0))

        t0 = time.time()
        F = compute_forces(L0, L, bars, barvec, p)
        if print_timing:
            print('time spent in compute_forces: {0:.2f}'.format(time.time() - t0))

        # update node positions
        t0 = time.time()
        for i in range(len(p)):
            p[i][0] += delta_t * F[i][0]
            p[i][1] += delta_t * F[i][1]
        if print_timing:
            print('time spent in forces and update: {0:.2f}'.format(time.time() - t0))

        t0 = time.time()
        within_bounds = within_bounds_function(p)
        if print_timing:
            print('time spent in within_bounds: {0:.2f}'.format(time.time() - t0))

        t0 = time.time()
        p = bring_outside_points_back_to_boundary(p, within_bounds, deps, distance_function)
        if print_timing:
            print('time spent in bring_outside_points_back_to_boundary: {0:.2f}'.format(time.time() - t0))

        if print_timing:
            print('time spent in iter: {0:.2f}'.format(time.time() - t0_iter))

        count += 1
        num_iterations_left = config['num_iterations'] - count

        # estimate the time left
        time_spent_iteration = time.time() - t0_iter
        time_left = time_spent_iteration * num_iterations_left
        if time_left > 60.0:
            time_left_string = '{0:.1f} minutes'.format(time_left / 60.0)
        else:
            time_left_string = '{0:.1f} seconds'.format(time_left)

        print('    completed iteration {0} ({1} to go, approx. {2} left)'.format(count, num_iterations_left, time_left_string))

        if num_iterations_left < 1:
            break

    print('\nnumber of points: {0}'.format(len(p)))
    print('\nnumber of triangles: {0}'.format(len(t)))

    return p, t


def get_polygon_length(polygon):
    l = 0.0
    for (p1, p2) in zip(polygon, polygon[1:]):
        l += get_distance(p1, p2)
    return l


def get_boundary_length(boundary_file_name, island_file_names):
    l = 0.0
    boundary_points = read_points(boundary_file_name)[0]  # there is only one boundary (right?)
    l += get_polygon_length(boundary_points)
    for island_file in island_file_names:
        for islands_points in read_points(island_file):
            l += get_polygon_length(islands_points)
    return l


def run(boundary_file_name,
        island_file_names,
        config_file_name,
        resolution_file_names=None,
        restart_file_name=None):

    with open(config_file_name, 'r') as f:
        try:
            config = yaml.load(f, yaml.SafeLoader)
        except yaml.YAMLError as exc:
            print(exc)
            sys.exit(-1)

    print('this output generated by smeshing v{0}'.format(__version__))

    if 'interpolation_step_length' in config:
        interpolation_step_length = config['interpolation_step_length']
        if 'num_interpolation_points' in config:
            sys.stderr.write('ERROR: options interpolation_step_length and num_interpolation_points conflict\n')
            sys.exit(1)
    else:
        boundary_length = get_boundary_length(boundary_file_name, island_file_names)
        interpolation_step_length = boundary_length / config['num_interpolation_points']

    _boundary_points = read_points(boundary_file_name)[0]  # there is only one boundary (right?)
    boundary_points = interpolate_polygon(_boundary_points, interpolation_step_length)
    all_points = boundary_points
    for island_file in island_file_names:
        for _islands_points in read_points(island_file):
            islands_points = interpolate_polygon(_islands_points, interpolation_step_length)
            all_points += islands_points
    num_points = len(all_points)

    print('\nnumber of boundary interpolation points: {0}'.format(len(all_points)))

    # pass only x and y coordinates and skip all the rest
    all_points_xy = [[t[0], t[1]] for t in all_points]

    all_polygons_context = polygons.new_context(1)
    boundary_context = polygons.new_context(0)
    islands_context = polygons.new_context(0)

    _boundary_points = read_points(boundary_file_name)[0]  # there is only one boundary (right?)
    boundary_points = interpolate_polygon(_boundary_points, interpolation_step_length)
    index_off = 0
    indices = list(range(index_off, index_off + len(boundary_points)))
    index_off += len(boundary_points)
    # pass only x and y coordinates and skip all the rest
    boundary_points_xy = [[t[0], t[1]] for t in boundary_points]

    distance_fields = []
    distance_field_bounds = []
    if resolution_file_names is not None:
        for file_name in resolution_file_names:
            xs, ys, zs = [], [], []
            for points in read_points(file_name):
                for (x, y, z) in points:
                    xs.append(x)
                    ys.append(y)
                    zs.append(z)
            distance_field_bounds.append((min(zs), max(zs)))
            distance_fields.append(SmoothBivariateSpline(xs, ys, zs))

    # FIXME generalize to arbitrary number of coefficients
    # now only one per line
    coefficients = [t[2] for t in boundary_points]

    polygons.add_polygon(all_polygons_context, boundary_points_xy, indices, coefficients)
    polygons.add_polygon(boundary_context, boundary_points_xy, indices, [])
    counter = len(boundary_points)
    all_points = boundary_points

    index_off_islands = 0
    for island_file in island_file_names:
        for _islands_points in read_points(island_file):
            islands_points = interpolate_polygon(_islands_points, interpolation_step_length)
            indices = list(range(index_off, index_off + len(islands_points)))
            index_off += len(islands_points)
            # pass only x and y coordinates and skip all the rest
            islands_points_xy = [[t[0], t[1]] for t in islands_points]

            # FIXME generalize to arbitrary number of coefficients
            # now only one per line
            coefficients = [t[2] for t in islands_points]

            polygons.add_polygon(all_polygons_context, islands_points_xy, indices, coefficients)
            indices = list(range(index_off_islands, index_off_islands + len(islands_points)))
            index_off_islands += len(islands_points)
            polygons.add_polygon(islands_context, islands_points_xy, indices, [])
            counter += len(islands_points)
            all_points += islands_points

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

    use_interpolation = False
    if use_interpolation:
        interpolation_step = h0
        xs = [xmin + i * interpolation_step for i in range(int((xmax - xmin) / interpolation_step))]
        ys = [ymin + i * interpolation_step for i in range(int((ymax - ymin) / interpolation_step))]
        ps = []
        zs = []
        for x in xs:
            ps = [(x, y) for y in ys]
            zs.append(polygons.get_distances_vertex_custom(all_polygons_context, ps))
        interp_spline = RectBivariateSpline(xs, ys, np.array(zs))

        def h_function(points):
            _distances = [interp_spline(x, y)[0][0] for (x, y) in points]
            if len(distance_fields) > 0:
                _distances = minimize_over_resolution_fields(distance_fields,
                                                             distance_field_bounds,
                                                             _distances,
                                                             points)
            return _distances
    else:
        def h_function(points):
            _distances = polygons.get_distances_vertex_custom(all_polygons_context, points)
            if len(distance_fields) > 0:
                _distances = minimize_over_resolution_fields(distance_fields,
                                                             distance_field_bounds,
                                                             _distances,
                                                             points)
            return _distances

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


def interpolate_points(p1, p2, distance_from_p1):
    d = get_distance(p1, p2)
    d1 = distance_from_p1
    d2 = d - d1
    new_point = []
    for i in range(len(p1)):
        x = d1 / d * p2[i] + d2 / d * p1[i]
        new_point.append(x)
    return new_point


def interpolate_polygon(points, step_length):
    '''
    Walk along the polygon and interpolate it with points
    which are step_length apart.
    '''
    l_total = get_polygon_length(points)
    if step_length > l_total:
        return points

    l = 0.0
    lengths = [l]
    for (p1, p2) in zip(points, points[1:]):
        l += get_distance(p1, p2)
        lengths.append(l)

    interpolated_points = [points[0]]
    current_step = 0.0
    while True:
        current_step += step_length
        if current_step > l_total:
            break
        i = 0
        for (p1, p2) in zip(points, points[1:]):
            i += 1
            if current_step < lengths[i]:
                l = current_step - lengths[i - 1]
                new_point = interpolate_points(p1, p2, l)
                interpolated_points.append(new_point)
                break
    return interpolated_points + [points[0]]


def read_points(file_name):
    polygons = []
    with open(file_name, 'r') as f:
        for line in f:
            points = []
            num_points = int(line)
            for _ in range(num_points):
                line = next(f)
                point = [float(element) for element in line.split()]
                points.append(point)
            polygons.append(points)
    return polygons


def get_distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2.0 + (p2[1] - p1[1])**2.0)
