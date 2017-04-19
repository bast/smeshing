# encoding: utf-8
"""DistMesh 2D"""

#-----------------------------------------------------------------------------
#  Copyright (C) 2004-2012 Per-Olof Persson
#  Copyright (C) 2012 Bradley Froehle

#  Distributed under the terms of the GNU General Public License. You should
#  have received a copy of the license along with this program. If not,
#  see <http://www.gnu.org/licenses/>.
#-----------------------------------------------------------------------------

import numpy as np
import scipy.spatial as spspatial
import inpoly
import polygons

# Local imports
import mlcompat as ml
import utils as dmutils

from distance_functions import dpoly

__all__ = ['distmesh2d']

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def distmesh2d(pv, fh, h0, bbox, pfix=None):
    """
    distmesh2d: 2-D Mesh Generator using Distance Functions.

    Parameters
    ----------
    pv:        list of polygon coordinate tuples
    fh:        Scaled edge length function h(x,y)
    h0:        Initial edge length
    bbox:      Bounding box, (xmin, ymin, xmax, ymax)
    pfix:      Fixed node positions, shape (nfix, 2)

    Returns
    -------
    p:         Node positions (Nx2)
    t:         Triangle indices (NTx3)
    """

    inpoly_context = inpoly.new_context()
    polygons_context = polygons.new_context()

    inpoly.add_polygon(inpoly_context, pv)
    polygons.add_polygon(polygons_context, pv)

    dptol = .001
    ttol = .1
    Fscale = 1.2
    delta_t = 0.2
    geps = .001*h0
    deps = np.sqrt(np.finfo(np.double).eps)*h0
    densityctrlfreq = 30

    # Extract bounding box
    xmin, ymin, xmax, ymax = bbox
    if pfix is not None:
        pfix = np.array(pfix, dtype='d')

    # 1. Create initial distribution in bounding box (equilateral triangles)
    x, y = np.mgrid[xmin:(xmax+h0):h0,
                    ymin:(ymax+h0*np.sqrt(3)/2):h0*np.sqrt(3)/2]
    x[:, 1::2] += h0/2                               # Shift even rows
    p = np.vstack((x.flat, y.flat)).T                # List of node coordinates

    # 2. Remove points outside the region, apply the rejection method
    p = p[dpoly(p, pv, inpoly_context, polygons_context) < geps]                                # Keep only d<0 points
    r0 = 1/fh(p)**2                                  # Probability to keep point
    p = p[np.random.random(p.shape[0])<r0/r0.max()]  # Rejection method
    if pfix is not None:
        p = ml.setdiff_rows(p, pfix)                 # Remove duplicated nodes
        pfix = ml.unique_rows(pfix); nfix = pfix.shape[0]
        p = np.vstack((pfix, p))                     # Prepend fix points
    else:
        nfix = 0
    N = p.shape[0]                                   # Number of points N

    count = 0
    pold = float('inf')                              # For first iteration

    while True:
        count += 1

        # 3. Retriangulation by the Delaunay algorithm
        dist = lambda p1, p2: np.sqrt(((p1-p2)**2).sum(1))
        if (dist(p, pold)/h0).max() > ttol:          # Any large movement?
            pold = p.copy()                          # Save current positions
            t = spspatial.Delaunay(p).vertices       # List of triangles
            pmid = p[t].sum(1)/3                     # Compute centroids
            t = t[dpoly(pmid, pv, inpoly_context, polygons_context) < -geps]                  # Keep interior triangles
            # 4. Describe each bar by a unique pair of nodes
            bars = np.vstack((t[:, [0,1]],
                              t[:, [1,2]],
                              t[:, [2,0]]))          # Interior bars duplicated
            bars.sort(axis=1)
            bars = ml.unique_rows(bars)              # Bars as node pairs

        # 5. Move mesh points based on bar lengths L and forces F
        barvec = p[bars[:,0]] - p[bars[:,1]]         # List of bar vectors
        L = np.sqrt((barvec**2).sum(1))              # L = Bar lengths
        hbars = fh(p[bars].sum(1)/2)
        L0 = (hbars*Fscale
              *np.sqrt((L**2).sum()/(hbars**2).sum()))  # L0 = Desired lengths

        # Density control - remove points that are too close
        if (count % densityctrlfreq) == 0 and (L0 > 2*L).any():
            ixdel = np.setdiff1d(bars[L0 > 2*L].reshape(-1), np.arange(nfix))
            p = p[np.setdiff1d(np.arange(N), ixdel)]
            N = p.shape[0]; pold = float('inf')
            continue

        # Bar forces (scalars)
        F = L0 - L

        # we ignore attractive forces and only keep repulsion
        F[F < 0.0] = 0.0

        # compute forces along bars
        Fvec = np.zeros((len(F), 2))
        for i in range(len(F)):
            r = F[i]/L[i]
            Fvec[i][0] = r*barvec[i][0]
            Fvec[i][1] = r*barvec[i][1]

        # compute resulting force on each point from all adjacent bars
        Ftot = np.zeros((len(p), 2))
        for k, bar in enumerate(bars):
            i, j = bar
            Ftot[i][0] += Fvec[k][0]
            Ftot[i][1] += Fvec[k][1]
            Ftot[j][0] -= Fvec[k][0]
            Ftot[j][1] -= Fvec[k][1]

        # set force to zero at fixed points
        Ftot[:nfix] = 0.0

        # update node positions
        p += delta_t*Ftot

        # 6. Bring outside points back to the boundary
        d = dpoly(p, pv, inpoly_context, polygons_context); ix = d>0                          # Find points outside (d>0)
        if ix.any():
            dgradx = (dpoly(p[ix]+[deps,0], pv, inpoly_context, polygons_context)-d[ix])/deps # Numerical
            dgrady = (dpoly(p[ix]+[0,deps], pv, inpoly_context, polygons_context)-d[ix])/deps # gradient
            dgrad2 = dgradx**2 + dgrady**2
            p[ix] -= (d[ix]*np.vstack((dgradx, dgrady))/dgrad2).T # Project

        # 7. Termination criterion: All interior nodes move less than dptol (scaled)
        if (np.sqrt((delta_t*Ftot[d<-geps]**2).sum(1))/h0).max() < dptol:
            break

    # Clean up and plot final mesh
    p, t = dmutils.fixmesh(p, t)

    inpoly.free_context(inpoly_context)
    polygons.free_context(polygons_context)

    return p, t
