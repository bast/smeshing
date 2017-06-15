def read_data(file_name):
    points = []
    triangles = []

    with open(file_name, 'r') as f:
        lines = f.readlines()

        num_points = int(lines[0])
        for i in range(1, num_points + 1):
            x = float(lines[i].split()[0])
            y = float(lines[i].split()[1])
            points.append([x, y])

        num_triangles = int(lines[num_points + 1])
        for i in range(num_points + 2, num_points + num_triangles + 2):
            t1 = int(lines[i].split()[0])
            t2 = int(lines[i].split()[1])
            t3 = int(lines[i].split()[2])
            triangles.append([t1, t2, t3])

    return points, triangles


def write_data(points, triangles, file_name):
    with open(file_name, 'w') as f:
        f.write('{}\n'.format(len(points)))
        for x, y in points:
            f.write('{} {}\n'.format(x, y))
        f.write('{}\n'.format(len(triangles)))
        for t1, t2, t3 in triangles:
            f.write('{} {} {}\n'.format(t1, t2, t3))
