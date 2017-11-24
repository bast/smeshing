def resolution_function(distance_to_nearest_vertex,
                        nearest_distance_at_nearest_vertex,
                        r):
    s = 0.9
    resolution = s*distance_to_nearest_vertex + nearest_distance_at_nearest_vertex/r
    if resolution < 1.0:
        resolution = 1.0
    return resolution
