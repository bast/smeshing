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
import polygons
import math

# Local imports
import mlcompat as ml

#-----------------------------------------------------------------------------
# Functions
#-----------------------------------------------------------------------------

def distmesh2d(pv, fh, h0, bbox, pfix=None, max_num_iterations=None):
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

    polygons_context = polygons.new_context()

    polygons.add_polygon(polygons_context, pv)

    dptol = 0.001
    ttol = .1
    Fscale = 1.2
    delta_t = 0.2
    geps = .001*h0
    smallest_representable_double_float = np.finfo(np.double).eps
    deps = math.sqrt(smallest_representable_double_float)*h0
    densityctrlfreq = 30

    # Extract bounding box
    xmin, ymin, xmax, ymax = bbox
    if pfix is not None:
        pfix = np.array(pfix, dtype='d')

    # 1. Create initial distribution in bounding box (equilateral triangles)
    xs = [xmin]
    while xs[-1] <= xmax:
        xs.append(xs[-1] + h0)
    ys = [ymin]
    while ys[-1] <= ymax:
        ys.append(ys[-1] + h0*math.sqrt(3.0)/2.0)

    _points = []
    for x in xs:
        for row, y in enumerate(ys):
            if row%2 != 0:
                # shift every second row to the right
                _points.append([x + h0/2.0, y])
            else:
                _points.append([x, y])
    _points = np.array(_points)

    # 2. Remove points outside the region, apply the rejection method

    contains = polygons.contains_points(polygons_context, _points)
    p = []
    for i, point in enumerate(_points):
        if contains[i]:
            p.append(point)
    p = np.array(p)

    r0 = 1/fh(p)**2.0                                # Probability to keep point
    p = p[np.random.random(p.shape[0])<r0/r0.max()]  # Rejection method
    if pfix is not None:
        p = ml.setdiff_rows(p, pfix)                 # Remove duplicated nodes

        _pfix = pfix.tolist()
        _pfix = [tuple(x) for x in _pfix]
        _pfix = set(_pfix)
        _pfix = list(_pfix)
        pfix = np.array(_pfix)

        nfix = pfix.shape[0]
        p = np.vstack((pfix, p))                     # Prepend fix points
    else:
        nfix = 0

    count = 0
    pold = float('inf')                              # For first iteration

    while True:
        count += 1

        if max_num_iterations is not None:
            if count > max_num_iterations:
                break

        # 3. Retriangulation by the Delaunay algorithm
        dist = lambda p1, p2: np.sqrt(((p1-p2)**2).sum(1))
        if (dist(p, pold)/h0).max() > ttol:          # Any large movement?
            pold = p.copy()                          # Save current positions
            _triangles = spspatial.Delaunay(p).vertices       # List of triangles
            pmid = p[_triangles].sum(1)/3                     # Compute centroids

            # Keep interior triangles
            # in contrast to original implementation we do not use a tolerance at boundary
            # to avoid a distance computation
            contains = polygons.contains_points(polygons_context, pmid)
            t = []
            for i, triangle in enumerate(_triangles):
                if contains[i]:
                    t.append(triangle)
            t = np.array(t)

            # 4. Describe each bar by a unique pair of nodes
            bars = []
            for triangle in t:
                bars.append((triangle[0], triangle[1]))
                bars.append((triangle[1], triangle[2]))
                bars.append((triangle[2], triangle[0]))

            _bars = []
            for bar in bars:
                if bar[1] < bar[0]:
                    _bars.append((bar[1], bar[0]))
                else:
                    _bars.append(bar)
            _bars = set(_bars)
            _bars = list(_bars)
            _bars = sorted(_bars)
            bars = np.array(_bars)

        # 5. Move mesh points based on bar lengths L and forces F
        barvec = p[bars[:,0]] - p[bars[:,1]]         # List of bar vectors
        _barvec = []
        for bar in bars:
            vx = p[bar[0]][0] - p[bar[1]][0]
            vy = p[bar[0]][1] - p[bar[1]][1]
            _barvec.append([vx, vy])
        _L = []
        for bar in _barvec:
            _L.append(math.sqrt(bar[0]**2.0 + bar[1]**2.0))
        L = np.array(_L)

        bar_midpoints = []
        for bar in bars:
            _px = p[bar[0]][0] + p[bar[1]][0]
            _py = p[bar[0]][1] + p[bar[1]][1]
            bar_midpoints.append([_px/2.0, _py/2.0])
        bar_midpoints = np.array(bar_midpoints)
        hbars = fh(bar_midpoints)

        _L0 = []
        l2sum = 0.0
        hbars2sum = 0.0
        for i in range(len(_L)):
            l2sum += _L[i]**2.0
            hbars2sum += hbars[i]**2.0
        for i in range(len(_L)):
            _L0.append(hbars[i]*Fscale*math.sqrt(l2sum/hbars2sum))
        L0 = np.array(_L0)

        # Density control - remove points that are too close
        apply_density_control = False
        points_to_remove = []
        if count % densityctrlfreq == 0:
            for i in range(len(L0)):
                if L0[i] > 2.0*L[i]:
                    apply_density_control = True
                    for k in [0, 1]:
                        ip = bars[i][k]
                        if ip > nfix:
                            points_to_remove.append(ip)
        if apply_density_control:
            points_to_remove = list(set(points_to_remove))
            p = np.array([p[ip] for ip in range(len(p)) if ip not in points_to_remove])
            pold = float('inf')
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
        contains = polygons.contains_points(polygons_context, p)
        for i in range(len(p)):
            if not contains[i]:
                px = p[i] + [deps, 0]
                py = p[i] + [0, deps]
                d0 = polygons.get_distances(polygons_context, [p[i]])[0]
                dx = polygons.get_distances(polygons_context, [px])[0]
                dy = polygons.get_distances(polygons_context, [py])[0]
                dgradx = (dx - d0)/deps
                dgrady = (dy - d0)/deps
                dgrad2 = dgradx**2.0 + dgrady**2.0
                p[i][0] -= d0*dgradx/dgrad2
                p[i][1] -= d0*dgrady/dgrad2

        # 7. Termination criterion: All interior nodes move less than dptol (scaled)
        s = []
        for i in range(len(p)):
            if contains[i]:
                tmp = delta_t*Ftot[i]**2.0
                s.append(tmp[0] + tmp[1])
        if max(s) < (dptol*h0)**2.0:
            break

    polygons.free_context(polygons_context)

    return p, t
