import click
import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import polars as pl
from plotly.subplots import make_subplots
from scipy import stats

TEMPLATE = "seaborn"


def convert_units(df):
    rename_dict = {}
    exprs = []
    for col in df.columns:
        if "(km/s)" in col:
            new_col = col.replace("(km/s)", "(m/s)")
            rename_dict[col] = new_col
            exprs.append(pl.col(col) * 1000)
        elif "(km^2)" in col:
            new_col = col.replace("(km^2)", "(m^2)")
            rename_dict[col] = new_col
            exprs.append(pl.col(col) * 1_000_000)
        elif "(km)" in col:
            new_col = col.replace("(km)", "(m)")
            rename_dict[col] = new_col
            exprs.append(pl.col(col) * 1000)

    if exprs:
        df = df.with_columns(exprs).rename(rename_dict)
    return df

def autocorr(x: np.ndarray, max_lag: int) -> np.ndarray:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]

    if len(x) < 3:
        return np.array([])

    x = x - np.mean(x)
    denom = np.dot(x, x)

    if denom == 0.0:
        return np.zeros(max_lag + 1)

    max_lag = min(max_lag, len(x) - 1)

    return np.array([
        np.dot(x[:-lag], x[lag:]) / denom if lag > 0 else 1.0
        for lag in range(max_lag + 1)
    ])

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

        ric_fig = make_subplots(
            rows=2,
            cols=1,
            subplot_titles=["RIC position error (m)", "RIC velocity error (m/s)"],
            vertical_spacing=0.1,
        )
        this_ric_fig = px.line(
            ricdf,
            x="Epoch (UTC)",
            y=["Delta X (RIC) (m)", "Delta Y (RIC) (m)", "Delta Z (RIC) (m)"],
            template="seaborn",
        )
        for trace in this_ric_fig.data:
            ric_fig.add_trace(trace, row=1, col=1)

        this_ric_fig = px.line(
            ricdf,
            x="Epoch (UTC)",
            y=[
                "Delta Vx (RIC) (m/s)",
                "Delta Vy (RIC) (m/s)",
                "Delta Vz (RIC) (m/s)",
            ],
            template="seaborn",
        )
        for trace in this_ric_fig.data:
            ric_fig.add_trace(trace, row=2, col=1)
        ric_fig.show()

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
        vertical_spacing=0.1,
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
                    showlegend=True,
                ),
                row=idx,
                col=1,
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
                    showlegend=(trace_type not in legend_added),
                ),
                row=idx,
                col=1,
            )
            legend_added.add(trace_type)

        fig.update_yaxes(title_text=unit, row=idx, col=1)
    fig.update_layout(title_text="Measurement Residuals", template="seaborn")
    fig.update_xaxes(matches="x")
    fig.show()

    # Plot the RIC uncertainty for position and velocity
    sigma_fig = make_subplots(
        rows=2,
        cols=1,
        subplot_titles=[
            "RIC 1-sigma position uncertainty (m)",
            "RIC 1-sigma velocity uncertainty (m/s)",
        ],
        vertical_spacing=0.1,
    )
    this_fig = px.line(
        df,
        x="Epoch (UTC)",
        y=["Sigma X (RIC) (m)", "Sigma Y (RIC) (m)", "Sigma Z (RIC) (m)"],
        template=TEMPLATE,
    )
    for trace in this_fig.data:
        sigma_fig.add_trace(trace, row=1, col=1)

    this_fig = px.line(
        df,
        x="Epoch (UTC)",
        y=["Sigma Vx (RIC) (m/s)", "Sigma Vy (RIC) (m/s)", "Sigma Vz (RIC) (m/s)"],
        template=TEMPLATE,
    )
    for trace in this_fig.data:
        sigma_fig.add_trace(trace, row=2, col=1)
    sigma_fig.show()

    # Create one QQ plot per whitened residual (there is only one if scalar processing of measurements)
    whitened_resids = [c for c in df.columns if "Whitened residual" in c]

    # Plot the residual ratios
    for whitened_resid in whitened_resids:
        # 1. Create the subplot figure: 2 rows, 1 column
        # shared_xaxes=True links zooming/panning across both rows
        fig = make_subplots(
            rows=4,
            cols=2,
            shared_xaxes=False,
            vertical_spacing=0.08,
            horizontal_spacing=0.07,
            specs=[
                [{"colspan": 2}, None],
                [{"colspan": 2}, None],
                [{"colspan": 2}, None],
                [{}, {}],
            ],
            subplot_titles=(
                "Rejected Status",
                "By Tracker",
                "Autocorrelation by Tracker",
                "Accepted Residuals Histogram",
                "Accepted Residuals QQ Plot",
            ),
        )
        # Generate the "Rejected" plot
        # We use px to get the traces, then add them to the subplot
        fig_rej = px.scatter(
            df, x="Epoch (UTC)", y=whitened_resid, color="Residual Rejected"
        )
        for trace in fig_rej.data:
            fig.add_trace(trace, row=1, col=1)

        # Generate the "Tracker" plot
        fig_track = px.scatter(df, x="Epoch (UTC)", y=whitened_resid, color="Tracker")
        for trace in fig_track.data:
            # We set showlegend=True to avoid duplicate legend entries if needed
            fig.add_trace(trace, row=2, col=1)

        df_accepted = df.filter(~pl.col("Residual Rejected"))
        fig_hist = px.histogram(
            df_accepted, x=whitened_resid, color="Tracker", barmode="overlay"
        )
        for trace in fig_hist.data:
            trace.showlegend = False
            fig.add_trace(trace, row=4, col=1)

        # Generate the QQ plot
        sample = df_accepted[whitened_resid].drop_nulls().to_numpy()

        # QQ data
        (osm, osr), (slope, intercept, r) = stats.probplot(sample, dist="norm", fit=True)

        x_line = np.linspace(np.min(osm), np.max(osm), 200)

        # QQ scatter
        fig.add_trace(
            go.Scatter(
                x=osm,
                y=osr,
                mode="markers",
                name=f"{whitened_resid} QQ",
                marker=dict(color="blue"),
            ),
            row=4,
            col=2,
        )

        # Standard-normal expected line: this is the important one for whitened residuals
        fig.add_trace(
            go.Scatter(
                x=x_line,
                y=x_line,
                mode="lines",
                name="Expected N(0,1)",
                line=dict(color="red", dash="dash"),
            ),
            row=4,
            col=2,
        )

        # Optional fitted line, useful to see empirical bias/scale
        fig.add_trace(
            go.Scatter(
                x=x_line,
                y=slope * x_line + intercept,
                mode="lines",
                name=f"Fitted normal: μ={intercept:.3f}, σ={slope:.3f}, R={r:.3f}",
                line=dict(color="gray"),
            ),
            row=4,
            col=2,
        )

        for tracker in df_accepted["Tracker"].unique().to_list():
            tracker_df = (
                df_accepted
                .filter(pl.col("Tracker") == tracker)
                .sort("Epoch (UTC)")
            )

            x = tracker_df[whitened_resid].drop_nulls().to_numpy()

            if len(x) < 5:
                continue

            max_lag = 30
            rho = autocorr(x, max_lag=max_lag)
            lags = np.arange(len(rho))

            fig.add_trace(
                go.Bar(
                    x=lags,
                    y=rho,
                    name=f"{tracker} ACF",
                    showlegend=False,
                    opacity=0.7,
                ),
                row=3,
                col=1,
            )

            # Approximate 95% white-noise bounds
            bound = 1.96 / np.sqrt(len(x))

            fig.add_trace(
                go.Scatter(
                    x=[0, len(rho) - 1],
                    y=[bound, bound],
                    mode="lines",
                    line=dict(dash="dash"),
                    name=f"{tracker} +95%",
                    showlegend=False,
                ),
                row=3,
                col=1,
            )

            fig.add_trace(
                go.Scatter(
                    x=[0, len(rho) - 1],
                    y=[-bound, -bound],
                    mode="lines",
                    line=dict(dash="dash"),
                    name=f"{tracker} -95%",
                    showlegend=False,
                ),
                row=3,
                col=1,
            )

        # Global Layout Updates
        fig.update_layout(
            title_text=f"Comprehensive Residual Analysis: {whitened_resid}",
            template="seaborn",
        )
        # Update axes titles for clarity
        fig.update_yaxes(title_text=whitened_resid, row=1, col=1)
        fig.update_yaxes(title_text=whitened_resid, row=2, col=1)
        fig.update_yaxes(title_text="Count", row=3, col=1)
        fig.update_xaxes(title_text="Theoretical N(0,1) Quantiles", row=3, col=2)
        fig.update_yaxes(title_text="Sample Quantiles", row=3, col=2)

        fig.show()

    if wstats:
        for msr in msr_types:
            px.scatter(
                df,
                x="Epoch (UTC)",
                y=[f"Real observation: {msr}", f"Computed observation: {msr}"],
            ).show()

        gain_columns = [c for c in df.columns if "Gain" in c]
        fs_ratio_columns = [c for c in df.columns if "Filter-smoother ratio" in c]
        is_filter_run = len(df[gain_columns].drop_nulls()) > 0

        # Plot the filter gains or filter-smoother ratios
        if is_filter_run:
            px.scatter(df, x="Epoch (UTC)", y=gain_columns).show()
        else:
            px.scatter(df, x="Epoch (UTC)", y=fs_ratio_columns).show()


if __name__ == "__main__":
    main()
