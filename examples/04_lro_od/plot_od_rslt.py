import polars as pl
from scipy.stats import chi2
import numpy as np
import plotly.graph_objects as go
import plotly.express as px

if __name__ == "__main__":
    df = pl.read_parquet("./04_lro_od_results.parquet")

    df = df.with_columns(pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")).sort(
        "Epoch (UTC)", descending=False
    )
    # Add the +/- 3 sigmas on measurement noise
    df = df.with_columns(
        [
            (3.0 * pl.col("Measurement noise: Range (km)")).alias(
                "Measurement noise 3-Sigma: Range (km)"
            ),
            (-3.0 * pl.col("Measurement noise: Range (km)")).alias(
                "Measurement noise -3-Sigma: Range (km)"
            ),
        ]
    )
    df = df.with_columns(
        [
            (3.0 * pl.col("Measurement noise: Doppler (km/s)")).alias(
                "Measurement noise 3-Sigma: Doppler (km/s)"
            ),
            (-3.0 * pl.col("Measurement noise: Doppler (km/s)")).alias(
                "Measurement noise -3-Sigma: Doppler (km/s)"
            ),
        ]
    )

    # == Residual plots ==
    # Nyx uses the Mahalanobis distance for the residual ratios, so we test the goodness using the Chi Square distribution.
    freedoms = 2 # Two degrees of freedoms for the range and the range rate.
    x_chi = np.linspace(chi2.ppf(0.01, freedoms), chi2.ppf(0.99, freedoms), 100)
    y_chi = chi2.pdf(x_chi, freedoms)

    # Compute the scaling factor
    hist = np.histogram(df["Residual ratio"].fill_null(0.0), bins=50)[0]
    max_hist = max(hist[1:]) # Ignore the bin of zeros
    max_chi2_pdf = max(y_chi)
    scale_factor = max_hist / max_chi2_pdf

    fig = go.Figure()
    fig.add_trace(go.Scatter(x=x_chi, y=y_chi * scale_factor, mode="lines", name="Scaled Chi-Squared"))
    fig.add_trace(go.Histogram(x=df["Residual ratio"], nbinsx=100, name="Residual ratios"))
    fig.show()

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
    for msr in ["Range (km)", "Doppler (km/s)"]:
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
