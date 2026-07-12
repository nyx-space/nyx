import click
import polars as pl

from nyx_space.plots.md import orbital_elements


@click.command
@click.option("-p", "--path", type=str)
def plot_orbit_elements(path: str):
    """
    Plots the orbital elements: SMA, ECC, INC, RAAN, AOP, TA, True Longitude, AOL
    """

    df = pl.read_parquet(path)
    orbital_elements(df).show()


if __name__ == "__main__":
    plot_orbit_elements()
