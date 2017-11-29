import delaunay


def solve_delaunay(p, within_bounds_function):
    """
    Retriangulation by the Delaunay algorithm.
    """
    _triangles = delaunay.solve(p)

    triangle_centroids = []
    for triangle in _triangles:
        x = 0.0
        y = 0.0
        for ip in triangle:
            x += p[ip][0]
            y += p[ip][1]
        triangle_centroids.append((x / 3.0, y / 3.0))

    # keep interior triangles
    # in contrast to the original implementation we do not use
    # a tolerance at boundary to avoid a distance computation
    within_bounds = within_bounds_function(triangle_centroids)
    t = []
    for i, triangle in enumerate(_triangles):
        if within_bounds[i]:
            t.append(triangle)

    bars = form_bars(t)

    return bars, t


def form_bars(triangles):
    bars = []
    for triangle in triangles:
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
    return _bars
