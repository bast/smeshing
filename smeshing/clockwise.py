def edges_sum(points):
    """
    Returns positive float if polygon is clockwise.
    """
    s = 0.0
    for i in range(len(points)):
        j = (i + 1) % len(points)
        s += (points[j][0] - points[i][0]) * (points[j][1] + points[i][1])
    return s


def test_edges_sum():
    points = [(5.0, 0.0),
              (6.0, 4.0),
              (4.0, 5.0),
              (1.0, 5.0),
              (1.0, 0.0)]
    assert edges_sum(points) == -44.0


def polygon_is_clockwise(points):
    return (edges_sum(points) > 0.0)
