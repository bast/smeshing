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
    """Implements the trivial uniform mesh size function h=1."""
    return 1.0


def solve(scale=1.0, benchmark=False):

    polygons_context = polygons.new_context()

    def distance_function(points):
        return polygons.get_distances(polygons_context, points)

    def contains_function(points):
        return polygons.contains_points(polygons_context, points)

    points = []
    with open('test/boundary.txt', 'r') as f:
        for line in f:
            x = float(line.split()[0])
            y = float(line.split()[1])
            points.append([scale * x, scale * y])
    polygons.add_polygon(polygons_context, points)

    xmin, xmax, ymin, ymax = get_bbox(points)

    if benchmark:
        f = huniform
        h0 = (xmax - xmin) / 80.0
        _p, _t = distmesh2d(points, f, distance_function, contains_function, scale * h0, points, max_num_iterations=100)
    else:
        f = huniform
        # f = lambda p: 0.05 + 0.3*dcircle(p, 0, 0, 0.01)
        h0 = (xmax - xmin) / 25.0
        _p, _t = distmesh2d(points, f, distance_function, contains_function, scale * h0, points, max_num_iterations=100)

    polygons.free_context(polygons_context)

    return _p, _t


def test_polygon():
    p, t = solve()
    if os.getenv('GENERATE_TESTS', False):
        write_data(p, t, 'test/result.txt')
    matches_with_reference(p, t, 'test/result.txt')


def test_bench():
    p, t = solve(benchmark=True)
    if os.getenv('GENERATE_TESTS', False):
        write_data(p, t, 'test/result-bench.txt')
    matches_with_reference(p, t, 'test/result-bench.txt')
