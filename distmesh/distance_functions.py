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


def dpoly(p, pv, inpoly_context, polygons_context):
    """Signed distance function for polygon with vertices pv.

    Usually pv should also be provided as fixed points in distmesh2d.

    pv should be provided as a list of coordinates [(x0,y0), (x1,y1), ...]
    or an array of shape (nv, 2).
    """
    contains = inpoly.contains_points(inpoly_context, p)

    ffi = FFI()
    distances_np = np.zeros(len(p), dtype=np.float64)
    distances_p = ffi.cast("double *", distances_np.ctypes.data)
    ps_x, ps_y = zip(*p)
    vs_x, vs_y = zip(*pv)
    # FIXME we should send in numpy arrays and not normal arrays
    # otherwise there can be memory issues for very large arrays
    lib.vdsegment(len(p), ps_x, ps_y, len(pv), vs_x, vs_y, distances_p)

    for i in range(len(p)):
        if contains[i]:
            distances_np[i] *= -1.0

    return distances_np


def huniform(p):
    """Implements the trivial uniform mesh size function h=1."""
    return np.ones(p.shape[0])
