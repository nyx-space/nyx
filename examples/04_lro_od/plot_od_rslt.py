import polars as pl
import plotly.express as px

if __name__ == "__main__":
    df = pl.read_parquet("./04_lro_od_results.parquet")

    df = df.with_columns(pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")).sort(
        "Epoch (UTC)", descending=False
    )
    # Add the Cr + Sigma Cr column
    df = df.with_columns((pl.col("cr") + pl.col("Sigma Cr (Moon J2000) (unitless)")).alias("Cr + Sigma"))
    df = df.with_columns((pl.col("cr") - pl.col("Sigma Cr (Moon J2000) (unitless)")).alias("Cr - Sigma"))

    # Plot the measurement residuals and their noises.
    for msr in ["Range (km)", "Doppler (km/s)"]:
        y_cols = [
            f"{col}: {msr}" for col in ["Prefit residual", "Measurement noise", "Postfit residual"]
        ]
        px.scatter(df, x="Epoch (UTC)", y=y_cols).show()

    # Plot the RIC uncertainty
    px.line(
        df, x="Epoch (UTC)", y=["Sigma X (RIC) (km)", "Sigma Y (RIC) (km)", "Sigma Z (RIC) (km)"]
    ).show()

    px.line(
        df,
        x="Epoch (UTC)",
        y=["Sigma Vx (RIC) (km/s)", "Sigma Vy (RIC) (km/s)", "Sigma Vz (RIC) (km/s)"],
    ).show()

    # Plot the Cr estimation
    px.line(df, x="Epoch (UTC)", y=["cr", "Cr + Sigma", "Cr - Sigma"]).show()

    # Load the RIC diff.
    df_ric = pl.read_parquet("./04_lro_od_truth_error.parquet")
    df_ric = df_ric.with_columns(
        pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
    ).sort("Epoch (UTC)", descending=False)
    # Plot the RIC difference
    px.line(
        df_ric,
        x="Epoch (UTC)",
        y=["Delta X (RIC) (km)", "Delta Y (RIC) (km)", "Delta Z (RIC) (km)"],
    ).show()
    px.line(
        df_ric,
        x="Epoch (UTC)",
        y=["Delta Vx (RIC) (km/s)", "Delta Vy (RIC) (km/s)", "Delta Vz (RIC) (km/s)"],
    ).show()
