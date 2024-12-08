import polars as pl
from scipy import stats
from scipy.stats import chi2
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

if __name__ == "__main__":
    df = pl.read_parquet("output_data/ekf_rng_dpl_az_el_odp.parquet")

    df = (
        df.with_columns(pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f"))
        .sort("Epoch (UTC)", descending=False)
        .with_columns(df["Residual ratio"].fill_null(0.0).alias("Residual ratio"))
    )

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

    # Convert the Polars column to a NumPy array for compatibility with scipy and Plotly
    residual_ratio = df["Residual ratio"].to_numpy()

    # Create QQ plot
    qq = stats.probplot(residual_ratio)
    x_qq = np.array([qq[0][0][0], qq[0][0][-1]])
    y_qq = np.array([qq[0][1][0], qq[0][1][-1]])

    # Create the QQ plot figure
    fig_qq = go.Figure()

    # Add scatter points
    fig_qq.add_trace(
        go.Scatter(
            x=qq[0][0],
            y=qq[0][1],
            mode='markers',
            name='Sample Data',
            marker=dict(color='blue')
        )
    )

    # Add the theoretical line
    fig_qq.add_trace(
        go.Scatter(
            x=x_qq,
            y=y_qq,
            mode='lines',
            name='Theoretical Normal',
            line=dict(color='red')
        )
    )

    # Update layout
    fig_qq.update_layout(
        title='Normal Q-Q Plot',
        xaxis_title='Theoretical Quantiles',
        yaxis_title='Sample Quantiles',
    )

    # Show QQ plot
    fig_qq.show()

    # Create histogram with normal distribution overlay
    hist_fig = px.histogram(
        df,
        x="Residual ratio",
        color="Tracker",
        marginal="rug",
        hover_data=df.columns,
    )

    # Calculate normal distribution parameters
    mean = residual_ratio.mean()
    std = residual_ratio.std()
    x_range = np.linspace(residual_ratio.min(), residual_ratio.max(), 100)
    y_normal = stats.norm.pdf(x_range, mean, std)

    # Scale the normal distribution to match histogram height
    max_hist_height = 100
    scaling_factor = max_hist_height / max(y_normal)
    y_normal_scaled = y_normal * scaling_factor

    # Add normal distribution curve
    hist_fig.add_trace(
        go.Scatter(
            x=x_range,
            y=y_normal_scaled,
            name='Normal Distribution',
            line=dict(color='red', width=2),
        )
    )

    # Show histogram with normal overlay
    hist_fig.show()

    px.histogram(
        df,
        x="Residual ratio",
        color="Tracker",
        marginal="rug",  # can be `box`, `violin`
        hover_data=df.columns,
    ).show()

    # Plot the residual ratios and whether they were accepted.
    px.scatter(df, x="Epoch (UTC)", y="Residual ratio", color="Residual Rejected").show()

    df_resid_ok = df.filter(df["Residual Rejected"] == False)

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
        fig = px.scatter(df_resid_ok, x="Epoch (UTC)", y=y_cols)
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

    # Plot the RIC uncertainty
    px.line(
        df, x="Epoch (UTC)", y=["Sigma X (RIC) (km)", "Sigma Y (RIC) (km)", "Sigma Z (RIC) (km)"]
    ).show()

    px.line(
        df,
        x="Epoch (UTC)",
        y=["Sigma Vx (RIC) (km/s)", "Sigma Vy (RIC) (km/s)", "Sigma Vz (RIC) (km/s)"],
    ).show()

    raise AssertionError("stop")
    # Load the RIC diff.
    for fname, errname in [
        ("04_lro_od_truth_error", "OD vs Flown"),
        ("04_lro_sim_truth_error", "Sim vs Flown (model matching)"),
    ]:
        df_ric = pl.read_parquet(f"./{fname}.parquet")
        df_ric = df_ric.with_columns(
            pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
        ).sort("Epoch (UTC)", descending=False)
        # Compute the range and range rate columns
        df_ric = df_ric.with_columns(
            [
                (
                    (
                        pl.col("Delta X (RIC) (km)") ** 2
                        + pl.col("Delta Y (RIC) (km)") ** 2
                        + pl.col("Delta Z (RIC) (km)") ** 2
                    )
                    ** 0.5
                ).alias("RIC Range (km)"),
                (
                    (
                        pl.col("Delta Vx (RIC) (km/s)") ** 2
                        + pl.col("Delta Vy (RIC) (km/s)") ** 2
                        + pl.col("Delta Vz (RIC) (km/s)") ** 2
                    )
                    ** 0.5
                ).alias("RIC Range Rate (km/s)"),
            ]
        )

        print(f"== {errname} ({fname}) ==")
        print("RIC Range (km)")
        print(df_ric["RIC Range (km)"].describe())
        print("RIC Range Rate (km/s)")
        print(df_ric["RIC Range Rate (km/s)"].describe())

        # Plot the RIC difference
        px.line(
            df_ric,
            x="Epoch (UTC)",
            y=["Delta X (RIC) (km)", "Delta Y (RIC) (km)", "Delta Z (RIC) (km)", "RIC Range (km)"],
            title=f"Position error with {errname} ({fname})",
        ).show()
        px.line(
            df_ric,
            x="Epoch (UTC)",
            y=[
                "Delta Vx (RIC) (km/s)",
                "Delta Vy (RIC) (km/s)",
                "Delta Vz (RIC) (km/s)",
                "RIC Range Rate (km/s)",
            ],
            title=f"Velocity error with {errname} ({fname})",
        ).show()
