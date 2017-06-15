import click
import glob
from smeshing.main import run


@click.command()
@click.option('--boundary', help='File containing boundary data.')
@click.option('--islands', help='Island file names (OK to use wildcards).')
@click.option('--config', help='Configuration file.')
def main(boundary, islands, config):
    points, triangles = run(boundary_file_name=boundary,
                            island_file_names=glob.glob(islands),
                            config_file_name=config)


if __name__ == '__main__':
    main()
