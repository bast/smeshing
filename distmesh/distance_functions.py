# encoding: utf-8
"""Distance functions."""

#-----------------------------------------------------------------------------
#  Copyright (C) 2004-2012 Per-Olof Persson
#  Copyright (C) 2012 Bradley Froehle

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import numpy as np
import inpoly
from c_interface import lib
from cffi import FFI

#-----------------------------------------------------------------------------
# Signed distance functions
#-----------------------------------------------------------------------------

def dcircle(p,xc,yc,r):
    """Signed distance to circle centered at xc, yc with radius r."""
    return np.sqrt(((p-np.array([xc,yc]))**2).sum(-1))-r


def dpoly(p, pv, context):
    """Signed distance function for polygon with vertices pv.

    Usually pv should also be provided as fixed points in distmesh2d.

    pv should be provided as a list of coordinates [(x0,y0), (x1,y1), ...]
    or an array of shape (nv, 2).
    """
    contains = inpoly.contains_points(context, p)

    result = dsegment(p, pv)
    for i in range(len(result)):
        if contains[i]:
            result[i] *= -1.0

    return result


def dsegment(ps, vs):
    """
    d = dsegment(p, v)

    Parameters
    ----------
    p : array, shape (np, 2)
        points
    v : array, shape (nv, 2)
        vertices of a closed array, whose edges are v[0]..v[1],
        ... v[nv-2]..v[nv-1]

    Output
    ------
    ds : array, shape (np, nv-1)
        distance from each point to each edge
    """

    ffi = FFI()
    distances_np = np.zeros(len(ps), dtype=np.float64)
    distances_p = ffi.cast("double *", distances_np.ctypes.data)
    ps_x, ps_y = zip(*ps)
    vs_x, vs_y = zip(*vs)
    # FIXME we should send in numpy arrays and not normal arrays
    # otherwise there can be memory issues for very large arrays
    lib.vdsegment(len(ps), ps_x, ps_y, len(vs), vs_x, vs_y, distances_p)
    return distances_np


def huniform(p):
    """Implements the trivial uniform mesh size function h=1."""
    return np.ones(p.shape[0])
