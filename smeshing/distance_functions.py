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
import polygons

#-----------------------------------------------------------------------------
# Signed distance functions
#-----------------------------------------------------------------------------

def dcircle(p,xc,yc,r):
    """Signed distance to circle centered at xc, yc with radius r."""
    return np.sqrt(((p-np.array([xc,yc]))**2).sum(-1))-r


def huniform(p):
    """Implements the trivial uniform mesh size function h=1."""
    return np.ones(p.shape[0])
