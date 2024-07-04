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
                    y=df_covar[col] + 3.0 * df_covar[f"Sigma {coord} (Earth J2000) (km)"],
                    mode="lines",
                    showlegend=True,
                    name=coord + " (km) + 3-Σ",
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col] - 3.0 * df_covar[f"Sigma {coord} (Earth J2000) (km)"],
                    mode="lines",
                    showlegend=True,
                    name=coord + " (km) - 3-Σ",
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
                    y=df_covar[col]
                    + 3.0 * df_covar[f"Sigma {coord.capitalize()} (Earth J2000) (km/s)"],
                    mode="lines",
                    showlegend=True,
                    name=coord + " (km/s) + 3-Σ",
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col]
                    - 3.0 * df_covar[f"Sigma {coord.capitalize()} (Earth J2000) (km/s)"],
                    mode="lines",
                    showlegend=True,
                    name=coord + " (km/s) - 3-Σ",
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

    # Build the RIC uncertainties plots
    fig = go.Figure(
        data=[
            go.Scattergl(
                x=df_covar["Epoch (UTC)"],
                y=3.0 * df_covar[col],
                mode="lines",
                showlegend=True,
                name=f"3-{col}",
            )
            for col in ["Sigma X (RIC) (km)", "Sigma Y (RIC) (km)", "Sigma Z (RIC) (km)"]
        ]
    )
    fig.update_layout(
        title=f"JWST RIC 3-Σ position prediction",
        xaxis_title="Epoch",
        legend_title="Legend",
    )
    fig.show()
    fig = go.Figure(
        data=[
            go.Scattergl(
                x=df_covar["Epoch (UTC)"],
                y=3.0 * df_covar[col],
                mode="lines",
                showlegend=True,
                name=f"3-{col}",
            )
            for col in ["Sigma Vx (RIC) (km/s)", "Sigma Vy (RIC) (km/s)", "Sigma Vz (RIC) (km/s)"]
        ]
    )
    fig.update_layout(
        title=f"JWST RIC 3-Σ position prediction",
        xaxis_title="Epoch",
        legend_title="Legend",
    )
    fig.show()

    # Build the Keplerian uncertainty plots
    keplerian_columns = [
        "sma (km)",
        "ecc",
        "inc (deg)",
        "raan (deg)",
        "aop (deg)",
        "ta (deg)",
        "energy (km^2/s^2)",
    ]
    for col in keplerian_columns:
        fig = go.Figure(
            data=[
                go.Scattergl(
                    x=df_mc["Epoch (UTC)"],
                    y=df_mc[col],
                    mode="lines",
                    opacity=0.05,
                    showlegend=True,
                    name=f"[MC] {col}",
                    text=df_mc["Monte Carlo Run Index"],
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col],
                    mode="lines",
                    showlegend=True,
                    name=f"[Nominal] {col}",
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col] + 3.0 * df_covar[f"Sigma {col}"],
                    mode="lines",
                    showlegend=True,
                    name=col + " + 3-Σ",
                ),
                go.Scattergl(
                    x=df_covar["Epoch (UTC)"],
                    y=df_covar[col] - 3.0 * df_covar[f"Sigma {col}"],
                    mode="lines",
                    showlegend=True,
                    name=col + " - 3-Σ",
                ),
            ]
        )
        fig.show()

    print("\nThe following columns are also provided as 1-sigma in the OD dataframe:")
    for col in df_covar.columns:
        skip = False
        for pre_plotted in keplerian_columns + ["X", "Y", "Z", "Vx", "Vy", "Vz"]:
            if pre_plotted in col:
                skip = True
                break
        if not skip and "Sigma" in col:
            print(f"\t- {col}")
