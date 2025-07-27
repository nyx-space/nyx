import click
import polars as pl
import plotly.express as px

@click.command
@click.option("-p", "--path", type=str, default="./04_lro_simulated_tracking.parquet")
def main(path: str):
    df = pl.read_parquet(path)

    df = df.with_columns(pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")).sort(
        "Epoch (UTC)", descending=False
    )

    # Build a title
    station_names = ", ".join([name for name in df["Tracking device"].unique()])
    start = df["Epoch (UTC)"][0]
    end = df["Epoch (UTC)"][-1]
    arc_duration = end - start
    title = f"Measurements from {station_names} spanning {start} to {end} ({arc_duration})"

    # Plot the overall tracking
    px.strip(df, x="Epoch (UTC)", y="Tracking device", color="Tracking device", title=title).show()

    # Plot each measurement kind
    for msr_col_name in ["Range (km)", "Doppler (km/s)"]:
        px.scatter(df, x="Epoch (UTC)", y=msr_col_name, color="Tracking device").show()

if __name__ == "__main__":
    main()
