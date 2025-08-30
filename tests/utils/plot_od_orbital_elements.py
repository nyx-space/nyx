import polars as pl
import plotly.graph_objs as go
from plotly.subplots import make_subplots
import click


@click.command
@click.option("-p", "--path", type=str)
@click.option("-m", "--multiplier", type=float, default=1.0)
def plot_orbit_elements(path: str, multiplier: int):
    """
    Plots the orbital elements: SMA, ECC, INC, RAAN, AOP, TA, True Longitude, AOL
    """

    df = pl.read_parquet(path)

    df = df.with_columns(
        pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
    ).sort("Epoch (UTC)", descending=False)

    columns = [
        "Sigma sma (km)",
        "Sigma ecc",
        "Sigma inc (deg)",
        "Sigma raan (deg)",
        "Sigma aop (deg)",
        "Sigma ta (deg)",
        "Sigma aol (deg)",
        "Sigma tlong (deg)",
    ]

    fig = make_subplots(
        rows=4,
        cols=2,
        subplot_titles=[f"{multiplier}-{col}" for col in columns],
        shared_xaxes=True,
        vertical_spacing=0.1,
    )

    row_i = 0
    col_i = 0

    for col in columns:
        name = f"{multiplier}-{col}"

        fig.add_trace(
            go.Scattergl(x=df["Epoch (UTC)"], y=df[col] * multiplier, name=name),
            row=row_i + 1,
            col=col_i + 1,
        )

        col_i = (col_i + 1) % 2
        if col_i == 0:
            row_i = (row_i + 1) % 4

    fig.show()


if __name__ == "__main__":
    plot_orbit_elements()
