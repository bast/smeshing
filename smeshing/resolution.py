import math


def linear_function(nearest_distance_at_coastline_point, distace_to_coastline_point):
    num_points_across_bay = 5
    resolution_at_coastline_point = nearest_distance_at_coastline_point/(num_points_across_bay + 1)
    slope = 0.995792
    return resolution_at_coastline_point + slope*distace_to_coastline_point


def tanh_function(nearest_distance_at_coastline_point, distace_to_coastline_point):
    num_points_across_bay = 5
    resolution_at_coastline_point = nearest_distance_at_coastline_point/(num_points_across_bay + 1)
    rcoast = resolution_at_coastline_point
    rfact = 2.0  # factor determining the near coastal length scales
    dfact = 3.0  # factor determining the middle resolution from rcoast and rmax.
    Ld    = 200.0  # typical length from coast to open boundary
    dev1 = 10.0  # deviding factor - near coast (default 4) higher number = steeper curve
    dev2 = 6.0
    rmax = 100.0
    r2 = rcoast + (rmax - rcoast) / rfact
    x1 = 3 * rcoast
    x2 = dfact * r2 + rcoast
    a1 = (x2 - x1) / dev1
    a2 = (Ld - x2) / dev2
    xm = x1 + 2 * a1
    xm2 = xm + 2 * a2
    _x = distace_to_coastline_point
    return rcoast + (r2 - rcoast)*(2 - (1 - math.tanh((_x-xm)/a1)))/2 + (rmax - r2)*(2 - (1 - math.tanh((_x - xm2)/a2)))/2
