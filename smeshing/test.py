# encoding: utf-8
"""Distmesh 2D examples."""

# Copyright (C) 2004-2012 Per-Olof Persson
# Copyright (C) 2012 Bradley Froehle

# Distributed under the terms of the GNU General Public License. You should
# have received a copy of the license along with this program. If not,
# see <http://www.gnu.org/licenses/>.

from .main import distmesh2d
from .file_io import read_data, write_data
import sys
import os


def get_bbox(points):

    xmin = sys.float_info.max
    xmax = -xmin
    ymin = xmin
    ymax = -ymin

    for point in points:
        xmin = min(xmin, point[0])
        xmax = max(xmax, point[0])
        ymin = min(ymin, point[1])
        ymax = max(ymax, point[1])

    return xmin, xmax, ymin, ymax


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
        assert abs((xs[i] - xs_ref[i])/xs[i]) < 1.0e-5
        assert abs((ys[i] - ys_ref[i])/ys[i]) < 1.0e-5


def huniform(x, y):
    """Implements the trivial uniform mesh size function h=1."""
    return 1.0


def read_polygon(file_name, scale=1.0, benchmark=False):

    pv = []
    with open(file_name, 'r') as f:
        for line in f:
            x = float(line.split()[0])
            y = float(line.split()[1])
            pv.append([scale*x, scale*y])

    xmin, xmax, ymin, ymax = get_bbox(pv)

    if benchmark:
        f = huniform
        h0 = (xmax - xmin)/80.0
        _p, _t = distmesh2d(pv, f, scale*h0, pv, max_num_iterations=100)
    else:
        f = huniform
        # f = lambda p: 0.05 + 0.3*dcircle(p, 0, 0, 0.01)
        h0 = (xmax - xmin)/25.0
        _p, _t = distmesh2d(pv, f, scale*h0, pv)
    return _p, _t


def test_polygon():
    p, t = read_polygon('test/input.txt', scale=1.0)
    if os.getenv('GENERATE_TESTS', False):
        write_data(p, t, 'test/result.txt')
    matches_with_reference(p, t, 'test/result.txt')


def test_bench():
    p, t = read_polygon('test/input.txt', benchmark=True)
    if os.getenv('GENERATE_TESTS', False):
        write_data(p, t, 'test/result-bench.txt')
    matches_with_reference(p, t, 'test/result-bench.txt')
