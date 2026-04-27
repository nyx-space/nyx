import polars as pl
import plotly.graph_objects as go
from plotly.subplots import make_subplots


if __name__ == "__main__":
    df = pl.read_parquet("03_geo_raise.parquet")

    df = df.with_columns(
        pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
    ).sort("Epoch (UTC)", descending=False)

    # We create two heatmaps to compare in-plane and out-of-plane evolution
    # Using SemiMajorAxis as the Y-axis "unrolls" the spiral orbit raise
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.1,
        subplot_titles=(
            "In-Plane Thrust Angle Evolution",
            "Out-of-Plane Thrust Angle Evolution",
        ),
    )

    # In-plane thrust (typically used for SMA and Eccentricity control)
    fig.add_trace(
        go.Histogram2d(
            x=df["AoL (deg)"],
            y=df["SemiMajorAxis (km)"],
            z=df["thrust_in_plane (RCN) (deg)"],
            histfunc="avg",
            nbinsx=72,  # 5-degree bins for AoL
            nbinsy=100,  # Granularity of the SMA progression
            colorscale="Viridis",
            colorbar=dict(title="In-Plane (deg)", len=0.45, y=0.75),
        ),
        row=1,
        col=1,
    )

    # Out-of-plane thrust (typically used for Inclination and RAAN control)
    fig.add_trace(
        go.Histogram2d(
            x=df["AoL (deg)"],
            y=df["SemiMajorAxis (km)"],
            z=df["thrust_out_of_plane (RCN) (deg)"],
            histfunc="avg",
            nbinsx=72,
            nbinsy=100,
            colorscale="RdBu",  # Diverging scale often better for out-of-plane
            colorbar=dict(title="Out-of-Plane (deg)", len=0.45, y=0.25),
        ),
        row=2,
        col=1,
    )

    fig.update_layout(
        height=800,
        title_text="Steering Law Evolution vs. Orbit Growth",
        xaxis2_title="Argument of Latitude (deg)",
        yaxis_title="Semi-Major Axis (km)",
        yaxis2_title="Semi-Major Axis (km)",
    )

    fig.show()
