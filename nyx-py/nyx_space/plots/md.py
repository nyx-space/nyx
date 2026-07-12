import plotly.express as px
import plotly.graph_objects as go
import polars as pl
from plotly.subplots import make_subplots
from nyx_space.plots import TEMPLATE, convert_units, watermark


def orbital_elements(df: pl.DataFrame, path: str | None = None) -> go.Figure:
    """
    Plot of the semi major axis, eccentricity, inclination, RAAN, Arg of Periapsis, True Anomaly,
    Arg of Latitude, and True Longitude
    """

    df = df.with_columns(
        pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
    ).sort("Epoch (UTC)", descending=False)

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

        fig.add_trace(
            go.Scattergl(x=df["Epoch (UTC)"], y=df[col], name=name),
            row=row_i + 1,
            col=col_i + 1,
        )
        col_i = (col_i + 1) % 2
        if col_i == 0:
            row_i = (row_i + 1) % 4

    return fig


def ric_diff(
    ricdf: pl.DataFrame, prefix="Delta ", path: str | None = None
) -> go.Figure:
    ricdf = ricdf.with_columns(
        pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
    ).sort("Epoch (UTC)", descending=False)
    ricdf = convert_units(ricdf)

    ric_fig = make_subplots(
        rows=2,
        cols=1,
        subplot_titles=["RIC position error (m)", "RIC velocity error (m/s)"],
        vertical_spacing=0.1,
    )
    this_ric_fig = px.line(
        ricdf,
        x="Epoch (UTC)",
        y=[f"{prefix}X (RIC) (m)", f"{prefix}Y (RIC) (m)", f"{prefix}Z (RIC) (m)"],
        template=TEMPLATE,
    )
    for trace in this_ric_fig.data:
        ric_fig.add_trace(trace, row=1, col=1)

    this_ric_fig = px.line(
        ricdf,
        x="Epoch (UTC)",
        y=[
            f"{prefix}Vx (RIC) (m/s)",
            f"{prefix}Vy (RIC) (m/s)",
            f"{prefix}Vz (RIC) (m/s)",
        ],
        template=TEMPLATE,
    )
    for trace in this_ric_fig.data:
        ric_fig.add_trace(trace, row=2, col=1)
    return watermark(ric_fig, path)
