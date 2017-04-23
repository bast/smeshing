# Copyright (C) 2004-2012 Per-Olof Persson
# Copyright (C) 2012 Bradley Froehle

# Distributed under the terms of the GNU General Public License. You should
# have received a copy of the license along with this program. If not,
# see <http://www.gnu.org/licenses/>.

import sys
import os
import math
import polygons
import flanders

from .main import distmesh2d
from .file_io import read_data, write_data
from .bbox import get_bbox
from .clockwise import edges_sum
from .resolution import linear_function, tanh_function


def normalize(vector, s):
    norm = math.sqrt(vector[0]**2.0 + vector[1]**2.0)
    return (s*vector[0]/norm, s*vector[1]/norm)


def get_normal_vectors(points, s):
    num_points = len(points)
    vectors = []
    for i in range(num_points):
        i_before = i - 1
        i_after = (i + 1)%num_points
        vector = (points[i_after][1] - points[i_before][1], -(points[i_after][0] - points[i_before][0]))
        vector = normalize(vector, s)
        vectors.append(vector)
    return vectors


def compute_view_vectors(points, scale):
    """
    If scale is negative, then view vectors are towards inside.
    """

    # we figure out whether polygon is clockwise or anticlockwise
    if edges_sum(points) < 0.0:
        s = +1.0*scale
    else:
        s = -1.0*scale

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
        assert abs((xs[i] - xs_ref[i]) / xs[i]) < 1.0e-5
        assert abs((ys[i] - ys_ref[i]) / ys[i]) < 1.0e-5


def get_distance(p1, p2):
    return math.sqrt((p2[0] - p1[0])**2.0 + (p2[1] - p1[1])**2.0)


def huniform(x, y):
    """
    Implements the trivial uniform mesh size function h=1.
    """
    return 1.0


def read_points(file_name):
    points = []
    with open(file_name, 'r') as f:
        for line in f:
            x = float(line.split()[0])
            y = float(line.split()[1])
            points.append([x, y])
    return points


def sub(file_name, benchmark=False):

    plot_nearest_in_view = False
    if plot_nearest_in_view:
        import matplotlib.pyplot as plt

    all_polygons_context = polygons.new_context()
    boundary_context = polygons.new_context()
    islands_context = polygons.new_context()

    boundary_points = read_points('test/boundary.txt')
    polygons.add_polygon(all_polygons_context, boundary_points)
    polygons.add_polygon(boundary_context, boundary_points)
    view_vectors = compute_view_vectors(boundary_points, scale=-1.0)
    all_points = boundary_points
    if plot_nearest_in_view:
        for i in range(len(boundary_points) - 1):
            plt.plot([boundary_points[i][0], boundary_points[i + 1][0]],
                     [boundary_points[i][1], boundary_points[i + 1][1]],
                     'r-')

    for island_file in ['test/island1.txt', 'test/island2.txt', 'test/island3.txt']:
        islands_points = read_points(island_file)
        polygons.add_polygon(all_polygons_context, islands_points)
        polygons.add_polygon(islands_context, islands_points)
        all_points += islands_points
        view_vectors += compute_view_vectors(islands_points, scale=1.0)
        if plot_nearest_in_view:
            for i in range(len(islands_points) - 1):
                plt.plot([islands_points[i][0], islands_points[i + 1][0]],
                         [islands_points[i][1], islands_points[i + 1][1]],
                         'b-')

    def distance_function(points):
        return polygons.get_distances(all_polygons_context, points)

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

    if benchmark:
        h_function = huniform
        h0 = (xmax - xmin) / 80.0
    else:
        h_function = huniform
        h0 = (xmax - xmin) / 25.0

    # currently not used
    num_points = len(all_points)
    flanders_context = flanders.new_context(num_points, all_points)
    angles_deg = [90.0 for _ in range(num_points)]
    flanders_indices = flanders.search_neighbor(context=flanders_context,
                                                ref_indices=list(range(num_points)),
                                                view_vectors=view_vectors,
                                                angles_deg=angles_deg)
#   print(flanders_indices)
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

    points, triangles = distmesh2d(all_points,
                                   h_function,
                                   distance_function,
                                   within_bounds_function,
                                   h0,
                                   all_points,
                                   max_num_iterations=100)

    polygons.free_context(all_polygons_context)
    polygons.free_context(boundary_context)
    polygons.free_context(islands_context)

    flanders.free_context(flanders_context)

    if os.getenv('GENERATE_TESTS', False):
        write_data(points, triangles, file_name)
    matches_with_reference(points, triangles, file_name)


def test_polygon():
    sub('test/result.txt')


def test_bench():
    sub('test/result-bench.txt', benchmark=True)


def test_resolution():
    plot_nearest_in_view = False

    if plot_nearest_in_view:
        import matplotlib.pyplot as plt

    points = []
    view_vectors = []

    for island_file in ['data/1.txt', 'data/2.txt', 'data/3.txt', 'data/4.txt']:
        islands_points = read_points(island_file)
        points += islands_points
        view_vectors += compute_view_vectors(islands_points, scale=1.0)
        if plot_nearest_in_view:
            for i in range(len(islands_points) - 1):
                plt.plot([islands_points[i][0], islands_points[i + 1][0]],
                         [islands_points[i][1], islands_points[i + 1][1]],
                         'b-')

    num_points = len(points)
    flanders_context = flanders.new_context(num_points, points)

    angles_deg = [90.0 for _ in range(num_points)]
    flanders_indices = flanders.search_neighbor(context=flanders_context,
                                                ref_indices=list(range(num_points)),
                                                view_vectors=view_vectors,
                                                angles_deg=angles_deg)

    for (x, y, r, f) in [(69.731182, 70.688529, 10.55650049085928, linear_function),
                         (69.731182, 70.688529, 4.3657, tanh_function),
                         (29.7312, 41.3754, 2.5481, tanh_function)]:
        _r = sys.float_info.max
        for i in range(len(points)):
            nearest_distance_at_coastline_point = get_distance(points[i], points[flanders_indices[i]])
            distace_to_coastline_point = get_distance(points[i], (x, y))
            _r = min(_r, f(nearest_distance_at_coastline_point, distace_to_coastline_point))
        assert abs(_r - r) < 1.0e-4

    if plot_nearest_in_view:
        for i in range(len(points)):
            if flanders_indices[i] > -1:
                plt.plot([points[i][0], points[flanders_indices[i]][0]],
                         [points[i][1], points[flanders_indices[i]][1]],
                         'k-')
                print(points[i], flanders_indices[i], get_distance(points[i], points[flanders_indices[i]]))
        plt.savefig('foo.png')
    flanders.free_context(flanders_context)
