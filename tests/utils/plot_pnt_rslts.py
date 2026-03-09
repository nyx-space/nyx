import polars as pl
import numpy as np
import plotly.graph_objects as go
from plotly.colors import sample_colorscale
import plotly.express as px
import click

color_scale = "phase"


@click.command
@click.option("-p", "--path", type=str)
def main(path: str):
    df = (
        pl.read_parquet(path)
        .with_columns(pl.col("Epoch (UTC)").cast(pl.Datetime))
        .sort("Epoch (UTC)", descending=False)
    )

    # Plot the evolution of PNT results with ellipses
    df_plot = df.filter(pl.col("Tracker").is_not_null()).with_columns(
        minutes=(pl.col("Epoch (UTC)") - pl.col("Epoch (UTC)").min()).dt.total_minutes()
    )

    fig = go.Figure()

    # Plot the trajectory
    fig.add_trace(
        go.Scatter(
            x=df_plot["Longitude (deg)"],
            y=df_plot["Latitude (deg)"],
            mode="lines+markers",
            marker=dict(
                size=6,
                color=df_plot["minutes"],
                colorscale=color_scale,
                showscale=True,
                # colorbar=dict(title="Time Index"),
            ),
            line=dict(color="lightgray", width=1),
            name="Trajectory (TRK minutes)",
            text=df_plot["Epoch (UTC)"],
            hovertemplate="Lon: %{x}<br>Lat: %{y}<br>Time: %{text}<extra></extra>",
        )
    )

    # Add uncertainty ellipses at selected points (e.g., every 10th point)
    sample_indices = range(0, len(df_plot), 10)
    max_index = df_plot["minutes"].max()

    for i in sample_indices:
        row = df_plot[i]

        # Create circle approximation for uncertainty
        theta = np.linspace(0, 2 * np.pi, 30)

        # Scale factor for visualization (e.g., 1-sigma ellipse)
        sigma_lat = row["Sigma Latitude (deg)"][0] * 3
        sigma_lon = row["Sigma Longitude (deg)"][0] * 3

        ellipse_x = row["Longitude (deg)"][0] + sigma_lon * np.cos(theta)
        ellipse_y = row["Latitude (deg)"][0] + sigma_lat * np.sin(theta)

        # Get color from colorscale
        normalized_index = row["minutes"][0] / max_index
        ellipse_color = sample_colorscale(color_scale, [normalized_index])[0]

        fig.add_trace(
            go.Scatter(
                x=ellipse_x,
                y=ellipse_y,
                mode="lines",
                line=dict(color=ellipse_color, width=2, dash="dash"),
                showlegend=False,
                hoverinfo="skip",
            )
        )

    fig.update_layout(
        title="Trajectory with Uncertainty Ellipses (3σ)",
        xaxis_title="Longitude (deg)",
        yaxis_title="Latitude (deg)",
        width=1200,
        height=1200,
    )

    fig.show()

    all_msr_types = ["Range (km)", "Doppler (km/s)", "Azimuth (deg)", "Elevation (deg)"]
    msr_type_count = 0
    msr_types = []

    for msr_type in all_msr_types:
        if f"Measurement noise: {msr_type}" in df.columns:
            print(f"Found data for {msr_type}")
            msr_type_count += 1
            msr_types += [msr_type]
            # Add the +/- 3 sigmas on measurement noise
            df = df.with_columns(
                [
                    (3.0 * pl.col(f"Measurement noise: {msr_type}")).alias(
                        f"Measurement noise 3-Sigma: {msr_type}"
                    ),
                    (-3.0 * pl.col(f"Measurement noise: {msr_type}")).alias(
                        f"Measurement noise -3-Sigma: {msr_type}"
                    ),
                ]
            )

    # Plot the measurement residuals and their noises.
    for msr in msr_types:
        y_cols = [
            f"{col}: {msr}"
            for col in [
                "Prefit residual",
                "Postfit residual",
                "Measurement noise 3-Sigma",
                "Measurement noise -3-Sigma",
            ]
        ]
        fig = px.scatter(df, x="Epoch (UTC)", y=y_cols)
        fig.update_traces(
            mode="lines",
            selector=dict(name=f"Measurement noise 3-Sigma: {msr}"),
            connectgaps=True,
            line=dict(dash="dash", color="black"),
        )
        fig.update_traces(
            mode="lines",
            selector=dict(name=f"Measurement noise -3-Sigma: {msr}"),
            connectgaps=True,
            line=dict(dash="dash", color="black"),
        )
        unit = msr.split()[-1][1:-1]
        fig.update_layout(yaxis_title=unit)
        fig.show()

        y_msr_cols = [
            f"{col}: {msr}"
            for col in [
                "Real observation",
                "Computed observation",
            ]
        ]
        px.scatter(df, x="Epoch (UTC)", y=y_msr_cols).show()


if __name__ == "__main__":
    main()
