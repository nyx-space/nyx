from datetime import datetime
import polars as pl
import plotly.graph_objs as go
from plotly.subplots import make_subplots


def plot_orbit_elements(
    df: pl.DataFrame,
):
    """
    Plots the orbital elements: SMA, ECC, INC, RAAN, AOP, TA, True Longitude, AOL
    """

    df = df.with_columns(pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")).sort(
        "Epoch (UTC)", descending=False
    )

    columns = [
        "sma (km)",
        "ecc",
        "inc (deg)",
        "raan (deg)",
        "aop (deg)",
        "ta (deg)",
        "aol (deg)",
        "tlong (deg)",
    ]

    fig = make_subplots(
        rows=4,
        cols=2,
        subplot_titles=columns,
        shared_xaxes=True,
        vertical_spacing=0.1,
    )

    row_i = 0
    col_i = 0

    for col in columns:
        name = col

        # # Build the color for this data frame
        # color = color_values[k % len(color_values)]
        # color = f"rgb({int(color[0])}, {int(color[1])}, {int(color[2])})"

        fig.add_trace(
            go.Scattergl(x=df["Epoch (UTC)"], y=df[col], name=name),  # , marker_color=color),
            row=row_i + 1,
            col=col_i + 1,
        )
        col_i = (col_i + 1) % 2
        if col_i == 0:
            row_i = (row_i + 1) % 4

    fig.show()


if __name__ == "__main__":
    plot_orbit_elements(pl.read_parquet("03_geo_raise.parquet"))
