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
    for fname, errname in [("04_lro_od_truth_error", "OD vs Flown"), ("04_lro_od_sim_error", "OD vs Sim"), ("04_lro_sim_truth_error", "Sim vs Flown")]:
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
