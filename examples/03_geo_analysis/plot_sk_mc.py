import polars as pl
import plotly.graph_objs as go


if __name__ == "__main__":
    df = pl.read_parquet("03_geo_sk.parquet")

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
        "fuel_mass (kg)",
    ]

    for col in columns:
        go.Figure(data=[go.Scattergl(x=df["Epoch (UTC)"], y=df[col], name=col, opacity=0.05, showlegend=True,text=df["Monte Carlo Run Index"])]).show()
