import polars as pl
import plotly.graph_objs as go
from plotly.subplots import make_subplots

if __name__ == "__main__":
    df_prop = pl.read_parquet("./03_geo_hf_prop.parquet")
    df_prop = df_prop.with_columns(
        pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
    ).sort("Epoch (UTC)", descending=False)
    print(df_prop.describe())

    df_lla = pl.read_parquet("./03_geo_lla.parquet")
    df_lla = df_lla.with_columns(
        pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
    ).sort("Epoch (UTC)", descending=False)
    print(df_lla.describe())

    # Create subplots
    fig = make_subplots(
        rows=2,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.1,
        subplot_titles=("Latitude vs Epoch", "Longitude vs Epoch"),
    )

    # Scattergl plot for Latitude vs Epoch
    fig.add_trace(
        go.Scattergl(
            x=df_lla["Epoch (UTC)"],
            y=df_lla["Latitude N-S (deg)"],
            mode="markers+lines",
            name="Latitude",
            marker=dict(color="blue"),
        ),
        row=1,
        col=1,
    )

    # Scattergl plot for Longitude vs Epoch
    fig.add_trace(
        go.Scattergl(
            x=df_lla["Epoch (UTC)"],
            y=df_lla["Longitude E-W (deg)"],
            mode="markers+lines",
            name="Longitude",
            marker=dict(color="blue"),
        ),
        row=2,
        col=1,
    )

    # Add vertical line to the Latitude plot
    first_latitude = df_lla["Latitude N-S (deg)"][0]
    lower_lat_deg = first_latitude - 0.05
    higher_lat_deg = first_latitude + 0.05

    fig.add_hline(
        lower_lat_deg,
        row=1,
        col=1,
        annotation_text=f"{lower_lat_deg:.3} - 0.05 deg",
        annotation_position="top right",
        line=dict(color="Red", width=1, dash="dash"),
    )
    fig.add_hline(
        higher_lat_deg,
        row=1,
        col=1,
        annotation_text=f"{lower_lat_deg:.3} + 0.05 deg",
        annotation_position="top right",
        line=dict(color="Red", width=1, dash="dash"),
    )

    # Add vertical line to the longitude plot
    first_longitude = df_lla["Longitude E-W (deg)"][0]
    lower_long_deg = first_longitude - 0.1
    higher_long_deg = first_longitude + 0.1

    fig.add_hline(
        lower_long_deg,
        row=2,
        col=1,
        annotation_text=f"{lower_long_deg:.3} - 0.1 deg",
        annotation_position="top right",
        line=dict(color="Red", width=1, dash="dash"),
    )
    fig.add_hline(
        higher_long_deg,
        row=2,
        col=1,
        annotation_text=f"{lower_lat_deg:.3} + 0.1 deg",
        annotation_position="top right",
        line=dict(color="Red", width=1, dash="dash"),
    )

    # Update layout
    fig.update_layout(
        title_text="Latitude and Longitude vs Epoch",
        xaxis_title="Epoch",
        yaxis_title="Latitude N-S (deg)",
        xaxis2_title="Epoch",
        yaxis2_title="Longitude E-W (deg)",
    )

    # Show the plot
    fig.show()
