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


def parse_command_line():
    import sys

    if len(sys.argv) != 3:
        sys.stderr.write('Usage: {0} [in-txt-file] [out-png-file]\n'.format(sys.argv[0]))
        sys.exit(1)

    in_file_name = sys.argv[1]
    out_file_name = sys.argv[2]

    return in_file_name, out_file_name


def generate_plot(x, y, triangles, out_file_name):
    import matplotlib.pyplot as plt

    plt.figure()
    plt.gca().set_aspect('equal')
    plt.triplot(x, y, triangles, 'go-', markersize=1.0, linewidth=0.8)
    plt.savefig(out_file_name)


in_file_name, out_file_name = parse_command_line()

x, y, triangles = read_data(in_file_name)

generate_plot(x, y, triangles, out_file_name)
