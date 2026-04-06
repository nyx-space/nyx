import polars as pl
from scipy import stats
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
from plotly.subplots import make_subplots
import click


def convert_units(df):
    rename_dict = {}
    exprs = []
    for col in df.columns:
        if "(km/s)" in col:
            new_col = col.replace("(km/s)", "(m/s)")
            rename_dict[col] = new_col
            exprs.append(pl.col(col) * 1000)
        elif "(km)" in col:
            new_col = col.replace("(km)", "(m)")
            rename_dict[col] = new_col
            exprs.append(pl.col(col) * 1000)

    if exprs:
        df = df.with_columns(exprs).rename(rename_dict)
    return df


@click.command
@click.option("-p", "--path", type=str)
@click.option("-s", "--wstats", type=bool, default=False)
@click.option("-e", "--error_ric", type=str, default=None)
def main(path: str, wstats: bool, error_ric: str):
    df = pl.read_parquet(path)

    df = df.with_columns(
        pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
    ).sort("Epoch (UTC)", descending=False)
    df = convert_units(df)

    if error_ric:
        ricdf = pl.read_parquet(error_ric)
        ricdf = ricdf.with_columns(
            pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
        ).sort("Epoch (UTC)", descending=False)
        ricdf = convert_units(ricdf)

        px.line(
            ricdf,
            x="Epoch (UTC)",
            y=["Delta X (RIC) (m)", "Delta Y (RIC) (m)", "Delta Z (RIC) (m)"],
        ).show()
        px.line(
            ricdf,
            x="Epoch (UTC)",
            y=["Delta Vx (RIC) (m/s)", "Delta Vy (RIC) (m/s)", "Delta Vz (RIC) (m/s)"],
        ).show()

    all_msr_types = ["Range (m)", "Doppler (m/s)", "Azimuth (deg)", "Elevation (deg)"]
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
    fig = make_subplots(
            rows=len(msr_types),
            cols=1,
            subplot_titles=[msr for msr in msr_types],
            vertical_spacing=0.1
        )
    legend_added = set()  # Track which trace names are already in legend

    for idx, msr in enumerate(msr_types, start=1):
        unit = msr.split()[-1][1:-1]
        y_cols = [
            f"{col}: {msr}"
            for col in [
                "Prefit residual",
                "Postfit residual",
                "Measurement noise 3-Sigma",
                "Measurement noise -3-Sigma",
            ]
        ]
        for y in y_cols[:2]:
            fig.add_trace(
                go.Scatter(
                    x=df["Epoch (UTC)"],
                    y=df[y],
                    mode="markers",
                    name=y,
                    legendgroup=y,
                    marker=dict(color="blue" if "Prefit" in y else "red"),
                    showlegend=True
                ),
                row=idx, col=1
            )

        # Add 3-sigma bounds
        for y in y_cols[2:]:
            trace_type = "3-Sigma bounds"
            fig.add_trace(
                go.Scatter(
                    x=df["Epoch (UTC)"],
                    y=df[y],
                    mode="lines",
                    name=trace_type,
                    line=dict(color="black"),
                    legendgroup=trace_type,
                    connectgaps=True,
                    showlegend=(trace_type not in legend_added)
                ),
                row=idx, col=1
            )
            legend_added.add(trace_type)

        fig.update_yaxes(title_text=unit, row=idx, col=1)
    fig.update_layout(title_text="Measurement Residuals")
    fig.update_xaxes(matches='x')
    fig.show()

    # Plot the residual ratios
    px.scatter(df, x="Epoch (UTC)", y="Residual ratio", color="Residual Rejected").show()
    px.scatter(df, x="Epoch (UTC)", y="Residual ratio", color="Tracker").show()

    # Plot the RIC uncertainty
    px.line(
        df,
        x="Epoch (UTC)",
        y=["Sigma X (RIC) (m)", "Sigma Y (RIC) (m)", "Sigma Z (RIC) (m)"],
    ).show()

    if wstats:
        for msr in msr_types:
            px.scatter(
                df,
                x="Epoch (UTC)",
                y=[f"Real observation: {msr}", f"Computed observation: {msr}"],
            ).show()

        # Convert the Polars column to a NumPy array for compatibility with scipy and Plotly
        residual_ratio = df["Residual ratio"].drop_nulls().to_numpy()

        gain_columns = [c for c in df.columns if "Gain" in c]
        fs_ratio_columns = [c for c in df.columns if "Filter-smoother ratio" in c]
        is_filter_run = len(df[gain_columns].drop_nulls()) > 0

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
                mode="markers",
                name="Residuals ratios (QQ)",
                marker=dict(color="blue"),
            )
        )

        # Add the theoretical line
        fig_qq.add_trace(
            go.Scatter(
                x=x_qq,
                y=y_qq,
                mode="lines",
                name="Theoretical Normal",
                line=dict(color="red"),
            )
        )

        # Update layout
        fig_qq.update_layout(
            title="Normal Q-Q Plot",
            xaxis_title="Theoretical Quantiles",
            yaxis_title="Sample Quantiles",
        )

        # Show QQ plot
        fig_qq.show()

        px.histogram(
            df,
            x="Residual ratio",
            color="Tracker",
            marginal="rug",  # can be `box`, `violin`
            hover_data=df.columns,
        ).show()

        # Plot the filter gains or filter-smoother ratios
        if is_filter_run:
            px.scatter(df, x="Epoch (UTC)", y=gain_columns).show()
        else:
            px.scatter(df, x="Epoch (UTC)", y=fs_ratio_columns).show()


if __name__ == "__main__":
    main()
