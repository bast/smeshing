import random


def generate_field(file_name, xs, ys, fun):
    with open(file_name, 'w') as f:
        f.write(f'{num_points}\n')
        for (x, y) in zip(xs, ys):
            z = fun(x, y)
            f.write(f'{x} {y} {z}\n')


if __name__ == '__main__':
    random.seed(0)
    num_points = 1000
    x_min, x_max = -10.0, 50.0
    y_min, y_max = 30.0, 80.0
    for (file_name, fun) in [('resolution-field-1.txt', lambda x, y: abs(y - 50.0)),
                             ('resolution-field-2.txt', lambda x, y: abs(x - 10.0))]:
        xs = [random.uniform(x_min, x_max) for _ in range(num_points)]
        ys = [random.uniform(y_min, y_max) for _ in range(num_points)]
        generate_field(file_name, xs, ys, fun)
