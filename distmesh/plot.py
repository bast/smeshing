import file_io


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

x, y, triangles = file_io.read_data(in_file_name)

generate_plot(x, y, triangles, out_file_name)
