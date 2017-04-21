from scipy.spatial import Delaunay
import random
import time
import matplotlib.pyplot as plt


def random_points(num_points):
    random.seed(1)
    points = [(random.uniform(0.0, 1.0),
               random.uniform(0.0, 1.0)) for _ in range(num_points)]
    return points


def plot(points, triangles):
    (x, y) = zip(*points)
    plt.figure()
    plt.gca().set_aspect('equal')
    plt.triplot(x, y, triangles, 'go-', markersize=1.0, linewidth=0.8)
    plt.savefig('triangles.png')


def main():

    for e in range(1, 7):
        num_points = 10**e
        points = random_points(num_points)
        t0 = time.time()
        triangles = Delaunay(points).vertices
        print('num points: {0} time spent: {1}'.format(num_points, time.time() - t0))

#   plot(points, triangles)


if __name__ == '__main__':
    main()
