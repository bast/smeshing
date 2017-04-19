from test import polygon
from file_io import write_data

# import cProfile
# cProfile.run("p, t = polygon('test/benchmark.txt', benchmark=True)")

p, t = polygon('test/benchmark.txt', benchmark=True)
write_data(p, t, 'test/data.txt')
