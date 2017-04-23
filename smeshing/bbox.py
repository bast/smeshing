import sys


def get_bbox(points):

    xmin = sys.float_info.max
    xmax = -xmin
    ymin = xmin
    ymax = -ymin

    for point in points:
        xmin = min(xmin, point[0])
        xmax = max(xmax, point[0])
        ymin = min(ymin, point[1])
        ymax = max(ymax, point[1])

    return xmin, xmax, ymin, ymax
