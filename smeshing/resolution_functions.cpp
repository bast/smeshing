#include <algorithm>
#include <limits>
#include <math.h>

#include "smeshing.h"

double linear_function(const double nearest_distance_at_coastline_point,
                       const double distace_to_coastline_point)
{
    int num_points_across_bay = 5;
    double resolution_at_coastline_point =
        nearest_distance_at_coastline_point / (num_points_across_bay + 1);
    double slope = 0.995792;

    return resolution_at_coastline_point + slope * distace_to_coastline_point;
}

double tanh_function(const double nearest_distance_at_coastline_point,
                     const double distace_to_coastline_point)
{
    int num_points_across_bay = 5;
    double resolution_at_coastline_point =
        nearest_distance_at_coastline_point / (num_points_across_bay + 1);
    double rcoast = resolution_at_coastline_point;
    double rfact = 2.0; // factor determining the near coastal length scales
    double dfact =
        3.0; // factor determining the middle resolution from rcoast and rmax.
    double Ld = 200.0;  // typical length from coast to open boundary
    double dev1 = 10.0; // deviding factor - near coast (default 4) higher
                        // number = steeper curve
    double dev2 = 6.0;
    double rmax = 100.0;
    double r2 = rcoast + (rmax - rcoast) / rfact;
    double x1 = 3 * rcoast;
    double x2 = dfact * r2 + rcoast;
    double a1 = (x2 - x1) / dev1;
    double a2 = (Ld - x2) / dev2;
    double xm = x1 + 2 * a1;
    double xm2 = xm + 2 * a2;
    double _x = distace_to_coastline_point;
    return rcoast + (r2 - rcoast) * (2 - (1 - tanh((_x - xm) / a1))) / 2 +
           (rmax - r2) * (2 - (1 - tanh((_x - xm2) / a2))) / 2;
}

double get_distance(const double p1x,
                    const double p1y,
                    const double p2x,
                    const double p2y)
{
    return sqrt(pow((p2x - p1x), 2.0) + pow((p2y - p1y), 2.0));
}

double get_resolution(const double x,
                      const double y,
                      const bool use_tanh,
                      const int num_points,
                      const double points_x[],
                      const double points_y[],
                      const int flanders_indices[])
{
    double r = std::numeric_limits<float>::max();
    for (int i = 0; i < num_points; i++)
    {
        // FIXME precompute this
        double nearest_distance_at_coastline_point =
            get_distance(points_x[i],
                         points_y[i],
                         points_x[flanders_indices[i]],
                         points_y[flanders_indices[i]]);

        double distace_to_coastline_point =
            get_distance(x, y, points_x[i], points_y[i]);

        if (use_tanh)
        {
            r = std::min(r,
                         tanh_function(nearest_distance_at_coastline_point,
                                       distace_to_coastline_point));
        }
        else
        {
            r = std::min(r,
                         linear_function(nearest_distance_at_coastline_point,
                                         distace_to_coastline_point));
        }
    }
    return r;
}
