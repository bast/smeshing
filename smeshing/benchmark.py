from test import polygon
from file_io import write_data

# import cProfile
# cProfile.run("p, t = polygon('test/benchmark.txt', benchmark=True)")

p, t = polygon('test/benchmark.txt', benchmark=True)
print('num points: {0}, num triangles: {1}'.format(len(p), len(t)))
write_data(p, t, 'test/data.txt')
