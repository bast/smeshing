import math
from .clockwise import polygon_is_clockwise


def compute_view_vectors(points, scale):
    """
    If scale is negative, then view vectors are towards inside.
    """

    if polygon_is_clockwise(points):
        s = -1.0 * scale
    else:
        s = +1.0 * scale

    # we remove the last point since it repeats the first
    # we assume clock-wise closed loops
    _points = [points[i] for i in range(len(points) - 1)]

    if len(_points) > 2:
        vectors = get_normal_vectors(_points, s)
    else:
        # this is a "linear" island consisting of two points
        vectors = []
        vector = (_points[1][0] - _points[0][0], _points[1][1] - _points[0][1])
        vectors.append(normalize(vector, -s))
        vectors.append(normalize(vector, s))

    # we add the removed point again
    return vectors + [vectors[0]]


def get_normal_vectors(points, s):
    num_points = len(points)
    vectors = []
    for i in range(num_points):
        i_before = i - 1
        i_after = (i + 1) % num_points
        vector = (points[i_after][1] - points[i_before][1], -(points[i_after][0] - points[i_before][0]))
        vector = normalize(vector, s)
        vectors.append(vector)
    return vectors


def normalize(vector, s):
    norm = math.sqrt(vector[0]**2.0 + vector[1]**2.0)
    return (s * vector[0] / norm, s * vector[1] / norm)
