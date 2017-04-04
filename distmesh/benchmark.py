from test import polygon, fstats
from file_io import write_data

p, t = polygon('test/benchmark.txt', uniform=True)
fstats(p, t)
write_data(p, t, 'test/data.txt')
