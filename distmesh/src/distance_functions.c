// Copyright (C) 2004-2012 Per-Olof Persson. See COPYING.TXT for details.

#include <math.h>
#include "distance_functions.h"


// Quick routines
static inline double length(double x, double y)
{
    return sqrt(x*x + y*y);
}


double dsegment(double x0,double y0,double p1x,double p1y,double p2x,double p2y)
{
  double v[2] = {p2x - p1x, p2y - p1y};
  double w[2] = {x0 - p1x, y0 - p1y};

  double c1 = v[0]*w[0] + v[1]*w[1];

  if (c1 <= 0.0)
  {
    return length(x0 - p1x, y0 - p1y);
  }

  double c2 = v[0]*v[0] + v[1]*v[1];

  if (c1 >= c2)
  {
    return length(x0 - p2x, y0 - p2y);
  }
  else
  {
    return length(x0 - (p1x + c1/c2*v[0]),
                  y0 - (p1y + c1/c2*v[1]));
  }
}
