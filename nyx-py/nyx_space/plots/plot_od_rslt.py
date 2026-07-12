import click
import polars as pl

from nyx_space.plots.od import (
    cr_cd,
    kalman_gains,
    filter_smoother_ratios,
    residuals,
    uncertainty,
    od_dashboard,
    orbital_element_uncertainty,
)
from nyx_space.plots.md import ric_diff

optional_est_params = ["cr", "cd"]


@click.command
@click.option("-p", "--path", type=str)
@click.option("-s", "--wstats", type=bool, default=False)
@click.option("-e", "--error_ric", type=str, default=None)
def main(path: str, wstats: bool, error_ric: str):
    df = pl.read_parquet(path)

    if error_ric:
        ricdf = pl.read_parquet(error_ric)
        ric_diff(ricdf).show()

    residuals(df, path).show()
    uncertainty(df, 3.0, path).show()
    for dash in od_dashboard(df, path):
        dash.show()

    cr_cd_plot = cr_cd(df, path)
    if cr_cd_plot is not None:
        cr_cd_plot.show()

    for plot in orbital_element_uncertainty(df, 3.0, path):
        plot.show()

    if wstats:
        gains_plot = kalman_gains(df, path)
        if gains_plot is not None:
            gains_plot.show()

        fs_ratios_plot = filter_smoother_ratios(df, path)
        if fs_ratios_plot is not None:
            fs_ratios_plot.show()


if __name__ == "__main__":
    main()
