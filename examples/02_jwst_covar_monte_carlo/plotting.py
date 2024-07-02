import polars as pl
import plotly.graph_objs as go

if __name__ == "__main__":
    df_mc = pl.read_parquet("./02_jwst_monte_carlo.parquet")
    df_mc = df_mc.with_columns(pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f"))
    print(df_mc.describe())

    df_covar = pl.read_parquet("./02_jwst_covar_map.parquet")
    df_covar = df_covar.with_columns(pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f"))
    print(df_covar.describe())

    # Build the position plots
    for coord in ["X", "Y", "Z"]:
        col = coord.lower() + " (km)"
        cov = coord + coord
        fig = go.Figure(
            data=[
                go.Scattergl(
                    x=df_mc["Epoch (UTC)"],
                    y=df_mc[col],
                    mode="lines",
                    opacity=0.05,
                    showlegend=True,
                    name=f"[MC] {coord} (km)",
                    text=df_mc["Monte Carlo Run Index"],
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col],
                    mode="lines",
                    showlegend=True,
                    name=f"[Nominal] {coord} (km)",
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col] + df_covar[f"Covariance {cov} (Earth J2000)"],
                    mode="lines",
                    showlegend=True,
                    name=coord + " (km) + Cov",
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col] - df_covar[f"Covariance {cov} (Earth J2000)"],
                    mode="lines",
                    showlegend=True,
                    name=coord + " (km) - Cov",
                ),
            ]
        )
        fig.update_layout(
            title=f"JWST {coord} state prediction - {len(df_mc)} Monte Carlo Results",
            xaxis_title="Epoch",
            yaxis_title=coord + " (km) EME2000",
            legend_title="Legend",
        )
        fig.show()

    # Build the velocity plots
    for coord in ["VX", "VY", "VZ"]:
        col = coord.lower() + " (km/s)"
        cov = coord.capitalize() + coord.capitalize()
        fig = go.Figure(
            data=[
                go.Scattergl(
                    x=df_mc["Epoch (UTC)"],
                    y=df_mc[col],
                    mode="lines",
                    opacity=0.05,
                    showlegend=True,
                    name=f"[MC] {coord} (km/s)",
                    text=df_mc["Monte Carlo Run Index"],
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col],
                    mode="lines",
                    showlegend=True,
                    name=f"[Nominal] {coord} (km/s)",
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col] + df_covar[f"Covariance {cov} (Earth J2000)"],
                    mode="lines",
                    showlegend=True,
                    name=coord + " (km/s) + Cov",
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col] - df_covar[f"Covariance {cov} (Earth J2000)"],
                    mode="lines",
                    showlegend=True,
                    name=coord + " (km/s) - Cov",
                ),
            ]
        )
        fig.update_layout(
            title=f"JWST {coord} velocity prediction - {len(df_mc)} Monte Carlo Results",
            xaxis_title="Epoch",
            yaxis_title=coord + " (km/s) EME2000",
            legend_title="Legend",
        )
        fig.show()
