# encoding: utf-8
"""Utilities for manipulating meshes."""

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

# Local imports.
import mlcompat as ml

__all__ = [
    'simpvol',
    ]

def simpvol(p, t):
    """Signed volumes of the simplex elements in the mesh."""
    dim = p.shape[1]
    if dim == 1:
        d01 = p[t[:,1]]-p[t[:,0]]
        return d01
    elif dim == 2:
        d01 = p[t[:,1]]-p[t[:,0]]
        d02 = p[t[:,2]]-p[t[:,0]]
        return (d01[:,0]*d02[:,1]-d01[:,1]*d02[:,0])/2
    else:
        raise NotImplementedError

def fixmesh(p, t, ptol=2e-13):
    """Remove duplicated/unused nodes and fix element orientation.

    Parameters
    ----------
    p : array, shape (np, dim)
    t : array, shape (nt, nf)

    Usage
    -----
    p, t = fixmesh(p, t, ptol)
    """
    snap = (p.max(0)-p.min(0)).max()*ptol
    _, ix, jx = ml.unique_rows(np.round(p/snap)*snap, True, True)

    p = p[ix]
    t = jx[t]

    flip = simpvol(p,t)<0
    t[flip, :2] = t[flip, 1::-1]

    return p, t
