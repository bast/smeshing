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

    assert len(x) == len(y)

    # remove triangles outside bounds
    # not sure how they can appear
    _triangles = []
    for triple in triangles:
        keep_triangle = True
        for t in triple:
            if not -1 < t < len(x):
                keep_triangle = False
        if keep_triangle:
                _triangles.append(triple)

    plt.figure()
    plt.gca().set_aspect('equal')
    plt.triplot(x, y, _triangles, 'g-', markersize=0.2, linewidth=0.2)
    plt.savefig(out_file_name, dpi=300)


in_file_name, out_file_name = parse_command_line()

x, y, triangles = file_io.read_data(in_file_name)

generate_plot(x, y, triangles, out_file_name)
