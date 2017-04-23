# Copyright (C) 2004-2012 Per-Olof Persson
# Copyright (C) 2012 Bradley Froehle

# Distributed under the terms of the GNU General Public License. You should
# have received a copy of the license along with this program. If not,
# see <http://www.gnu.org/licenses/>.

from .main import distmesh2d
from .file_io import read_data, write_data
from .bbox import get_bbox
import os
import polygons


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

    all_polygons_context = polygons.new_context()
    boundary_context = polygons.new_context()
    islands_context = polygons.new_context()

    boundary_points = read_points('test/boundary.txt')
    polygons.add_polygon(all_polygons_context, boundary_points)
    polygons.add_polygon(boundary_context, boundary_points)

    for island_file in ['test/island1.txt', 'test/island2.txt', 'test/island3.txt']:
        islands_points = read_points(island_file)
        polygons.add_polygon(all_polygons_context, islands_points)
        polygons.add_polygon(islands_context, islands_points)
        all_points = boundary_points + islands_points

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

    if os.getenv('GENERATE_TESTS', False):
        write_data(points, triangles, file_name)
    matches_with_reference(points, triangles, file_name)


def test_polygon():
    sub('test/result.txt')


def test_bench():
    sub('test/result-bench.txt', benchmark=True)
