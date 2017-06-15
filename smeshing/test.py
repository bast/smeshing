# Copyright (C) 2004-2012 Per-Olof Persson
# Copyright (C) 2012 Bradley Froehle

# Distributed under the terms of the GNU General Public License. You should
# have received a copy of the license along with this program. If not,
# see <http://www.gnu.org/licenses/>.

import sys
import os
import math
import flanders
import glob
from smeshing import get_resolution

from .main import distmesh2d, run, read_points, compute_view_vectors, get_distance
from .file_io import read_data, write_data
from .bbox import get_bbox


def matches_with_reference(ps, ts, file_name):
    """
    We compare triangle midpoints.
    """
    xs = []
    ys = []
    for p in ps:
        xs.append(p[0])
        ys.append(p[1])

    points_ref, ts_ref = read_data(file_name)
    for i, t in enumerate(ts):
        assert t[0] == ts_ref[i][0]
        assert t[1] == ts_ref[i][1]
        assert t[2] == ts_ref[i][2]

    for i in range(len(xs)):
        assert abs((xs[i] - points_ref[i][0]) / xs[i]) < 1.0e-5
        assert abs((ys[i] - points_ref[i][1]) / ys[i]) < 1.0e-5


if os.getenv('ONLY_LOFOTEN', False):
    reference_file_name = 'data/lofoten/result.txt'
    def test_lofoten():
        points, triangles = run(boundary_file_name='data/lofoten/boundary.txt',
                                island_file_names=glob.glob('data/lofoten/islands/*.txt'),
                                config_file_name='data/lofoten/config.yml')
        write_data(points, triangles, reference_file_name)
    def dont_test_lofoten_tiny():
        points, triangles = run(boundary_file_name='data/lofoten/simple-boundary.txt',
                                island_file_names=['data/lofoten/islands/{0}.txt'.format(i) for i in [275, 209, 38, 154, 19,
                                                                                                      247, 173, 210, 39, 95]],
                                config_file_name='data/lofoten/config.yml')
else:
    def test_polygon():
        reference_file_name = 'data/fiction/result.txt'
        points, triangles = run(boundary_file_name='data/fiction/boundary.txt',
                                island_file_names=['data/fiction/island1.txt', 'data/fiction/island2.txt', 'data/fiction/island3.txt'],
                                config_file_name='data/fiction/config.yml')
        if os.getenv('GENERATE_TESTS', False):
            write_data(points, triangles, reference_file_name)
        matches_with_reference(points, triangles, reference_file_name)


def test_resolution():
    plot_nearest_in_view = False

    if plot_nearest_in_view:
        import matplotlib.pyplot as plt

    points = []
    view_vectors = []

    for island_file in ['data/resolution/1.txt', 'data/resolution/2.txt', 'data/resolution/3.txt', 'data/resolution/4.txt']:
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

    nearest_distance_at_coastline_point = []
    for i in range(len(points)):
        nearest_distance_at_coastline_point.append(get_distance(points[i], points[flanders_indices[i]]))

    for (x, y, r) in [(69.731182, 70.688529, 10.55650049085928)]:
        _r = get_resolution([(x, y)], False, points, nearest_distance_at_coastline_point, flanders_indices)
        assert abs(_r[0] - r) < 1.0e-4

    for (x, y, r) in [(69.731182, 70.688529, 4.3657),
                      (29.7312, 41.3754, 2.5481)]:
        _r = get_resolution([(x, y)], True, points, nearest_distance_at_coastline_point, flanders_indices)
        assert abs(_r[0] - r) < 1.0e-4

    if plot_nearest_in_view:
        for i in range(len(points)):
            if flanders_indices[i] > -1:
                plt.plot([points[i][0], points[flanders_indices[i]][0]],
                         [points[i][1], points[flanders_indices[i]][1]],
                         'k-')
                print(points[i], flanders_indices[i], get_distance(points[i], points[flanders_indices[i]]))
        plt.savefig('foo.png')
    flanders.free_context(flanders_context)
