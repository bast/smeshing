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

from .main import run, read_points, compute_view_vectors, get_distance
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


def test_polygon():
    reference_file_name = 'data/happy-bear/result.txt'
    points, triangles = run(boundary_file_name='data/happy-bear/boundary.txt',
                            island_file_names=['data/happy-bear/island1.txt',
                                               'data/happy-bear/island2.txt',
                                               'data/happy-bear/island3.txt',
                                               'data/happy-bear/island4.txt'],
                            config_file_name='data/happy-bear/config.yml',
                            resolution_function_file_name='data/happy-bear/resolution_function.py')
    if os.getenv('GENERATE_TESTS', False):
        write_data(points, triangles, reference_file_name)
    matches_with_reference(points, triangles, reference_file_name)
