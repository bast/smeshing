def read_data(file_name):
    x = []
    y = []
    triangles = []

    with open(file_name, 'r') as f:
        lines = f.readlines()

        num_points = int(lines[0])
        for i in range(1, num_points+1):
            x.append(float(lines[i].split()[0]))
            y.append(float(lines[i].split()[1]))

        num_triangles = int(lines[num_points+1])
        for i in range(num_points+2, num_points+num_triangles+2):
            t1 = int(lines[i].split()[0])
            t2 = int(lines[i].split()[1])
            t3 = int(lines[i].split()[2])
            triangles.append([t1, t2, t3])

    return x, y, triangles
