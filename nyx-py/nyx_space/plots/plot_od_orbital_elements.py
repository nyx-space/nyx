import click
import plotly.graph_objs as go
import polars as pl
from plotly.subplots import make_subplots

from nyx_space.plots import TEMPLATE, watermark


@click.command
@click.option("-p", "--path", type=str)
@click.option("-m", "--multiplier", type=float, default=3.0)
def plot_orbit_elements(path: str, multiplier: int):
    """
    Plots the orbital elements: SMA, ECC, INC, RAAN, AOP, TA, True Longitude, AOL
    """

    df = pl.read_parquet(path)

    df = df.with_columns(
        pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
    ).sort("Epoch (UTC)", descending=False)

    columns = [
        "SemiMajorAxis (km)",
        "Eccentricity (unitless)",
        "Inclination (deg)",
        "RAAN (deg)",
        "AoP (deg)",
        "TrueAnomaly (deg)",
        "AoL (deg)",
        "TrueLongitude (deg)",
    ]

    for sigma in [False, True]:
        if sigma:
            subplot_titles = [f"{multiplier}-Sigma {col}" for col in columns]
        else:
            subplot_titles = columns
        fig = make_subplots(
            rows=4,
            cols=2,
            subplot_titles=subplot_titles,
            shared_xaxes=True,
            vertical_spacing=0.1,
        )

        row_i = 0
        col_i = 0

        for col in columns:
            if sigma:
                try:
                    y = df[f"Sigma {col}"] * multiplier
                except pl.exceptions.ColumnNotFoundError:
                    # No sigmas in this dataframe
                    return
                name = f"{multiplier}-{col}"
            else:
                y = df[col]
                name = col

            fig.add_trace(
                go.Scattergl(x=df["Epoch (UTC)"], y=y, name=name),
                row=row_i + 1,
                col=col_i + 1,
            )

            col_i = (col_i + 1) % 2
            if col_i == 0:
                row_i = (row_i + 1) % 4

        fig.update_layout(template=TEMPLATE)
        watermark(fig, path).show()


if __name__ == "__main__":
    plot_orbit_elements()
