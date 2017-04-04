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

__all__ = [
    # Distance functions:
    'dblock',
    'dcircle',
    'ddiff',
    'dintersect',
    'dmatrix3d',
    'dmatrix',
    'dpoly',
    'drectangle',
    'drectangle0',
    'dsegment',
    'dsphere',
    'dunion',

    # Mesh size functions:
    'hmatrix3d',
    'hmatrix',
    'huniform',

    # Generic node manipulation:
    'protate',
    'pshift',
    ]

# These are used very often:
min = np.minimum
max = np.maximum

#-----------------------------------------------------------------------------
# Signed distance functions
#-----------------------------------------------------------------------------

def dcircle(p,xc,yc,r):
    """Signed distance to circle centered at xc, yc with radius r."""
    return np.sqrt(((p-np.array([xc,yc]))**2).sum(-1))-r


def dpoly(p,pv):
    """Signed distance function for polygon with vertices pv.

    Usually pv should also be provided as fixed points in distmesh2d.

    pv should be provided as a list of coordinates [(x0,y0), (x1,y1), ...]
    or an array of shape (nv, 2).
    """
    context = inpoly.new_context()  # FIXME move this up
    inpoly.add_polygon(context, pv)
    contains = inpoly.contains_points(context, p)

    d = dsegment(p, pv)
    result = d.min(1)
    for i in range(len(result)):
        if contains[i]:
            result[i] *= -1.0

    inpoly.free_context(context)

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
    from c_interface import lib

    ds = []
    for iv in range(len(vs)- 1):
        p1x = vs[iv][0]
        p1y = vs[iv][1]
        p2x = vs[iv+1][0]
        p2y = vs[iv+1][1]

        d_row = []
        for i in range(len(ps)):
            d = lib.dsegment(ps[i][0], ps[i][1], p1x, p1y, p2x, p2y)
            d_row.append(d)
        ds.append(d_row)
    return np.transpose(np.asarray(ds))

def huniform(p):
    """Implements the trivial uniform mesh size function h=1."""
    return np.ones(p.shape[0])
