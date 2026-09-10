import numpy as np
import plotly.express as px
import plotly.graph_objects as go
import polars as pl
from plotly.subplots import make_subplots
from scipy import stats

from nyx_space.plots import TEMPLATE, convert_units, watermark

optional_est_params = ["cr", "cd"]

all_msr_types = ["Range (m)", "Doppler (m/s)", "Azimuth (deg)", "Elevation (deg)"]

# Canonical Seaborn deep colorway
SEABORN_COLORS = [
    "#4C72B0",
    "#55A868",
    "#C44E52",
    "#8172B2",
    "#CCB974",
    "#64B5CD",
    "#8C8C8C",
    "#E377C2",
    "#BCBD22",
    "#17BECF",
]


def _get_tracker_color_map(
    df: pl.DataFrame,
) -> tuple[list[str], dict[str, str]]:
    if "Tracker" not in df.columns:
        return [], {}
    trackers = sorted(df["Tracker"].drop_nulls().unique().to_list())
    return trackers, {
        trk: SEABORN_COLORS[i % len(SEABORN_COLORS)] for i, trk in enumerate(trackers)
    }


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

    return np.array(
        [
            np.dot(x[:-lag], x[lag:]) / denom if lag > 0 else 1.0
            for lag in range(max_lag + 1)
        ]
    )


def residuals(df: pl.DataFrame, path: str | None = None) -> go.Figure:
    """
    Plots physical residuals per measurement type decomposed by ground station tracker,
    retaining postfit as primary points, prefit as secondary markers, and pass-segmented
    3-sigma noise envelopes without interpolating across inter-pass gaps.
    """
    df = convert_units(df)
    msr_types = [
        msr for msr in all_msr_types if f"Measurement noise: {msr}" in df.columns
    ]
    trackers, tracker_colors = _get_tracker_color_map(df)

    fig = make_subplots(
        rows=len(msr_types),
        cols=1,
        subplot_titles=[f"{msr} Residuals by Station" for msr in msr_types],
        vertical_spacing=0.1,
    )

    for idx, msr in enumerate(msr_types, start=1):
        unit = msr.split()[-1][1:-1]
        postfit_col = f"Postfit residual: {msr}"
        prefit_col = f"Prefit residual: {msr}"
        noise_col = f"Measurement noise: {msr}"

        df_msr = df.filter(pl.col(postfit_col).is_not_null())
        if df_msr.height == 0:
            continue

        for trk in trackers:
            trk_df = df_msr.filter(pl.col("Tracker") == trk)
            if trk_df.height == 0:
                continue

            color = tracker_colors.get(trk, "#4C72B0")

            # Postfit residuals (solid markers)
            fig.add_trace(
                go.Scatter(
                    x=trk_df["Epoch (UTC)"],
                    y=trk_df[postfit_col],
                    mode="markers",
                    name=f"{trk} Postfit",
                    legendgroup=trk,
                    marker=dict(color=color, symbol="circle", size=6),
                    hovertemplate=(
                        f"<b>{trk}</b> Postfit<br>"
                        "Epoch: %{x}<br>"
                        f"Residual: %{{y:.4f}} {unit}<br>"
                        "<extra></extra>"
                    ),
                    showlegend=(idx == 1),
                ),
                row=idx,
                col=1,
            )

            # Prefit residuals (open markers, toggled off by default for visual clarity)
            if prefit_col in trk_df.columns:
                fig.add_trace(
                    go.Scatter(
                        x=trk_df["Epoch (UTC)"],
                        y=trk_df[prefit_col],
                        mode="markers",
                        name=f"{trk} Prefit",
                        legendgroup=trk,
                        marker=dict(color=color, symbol="circle-open", size=5),
                        visible="legendonly",
                        hovertemplate=(
                            f"<b>{trk}</b> Prefit<br>"
                            "Epoch: %{x}<br>"
                            f"Residual: %{{y:.4f}} {unit}<br>"
                            "<extra></extra>"
                        ),
                        showlegend=(idx == 1),
                    ),
                    row=idx,
                    col=1,
                )

        # 3-Sigma measurement noise envelope (do not connect gaps across tracking passes)
        if noise_col in df_msr.columns:
            for mult, name in [(3.0, "+3σ Noise"), (-3.0, "-3σ Noise")]:
                fig.add_trace(
                    go.Scatter(
                        x=df_msr["Epoch (UTC)"],
                        y=mult * df_msr[noise_col],
                        mode="lines",
                        name=name,
                        line=dict(
                            color="rgba(80, 80, 80, 0.6)", dash="dash", width=1.2
                        ),
                        legendgroup="Noise Bounds",
                        connectgaps=False,
                        showlegend=(idx == 1 and mult > 0),
                        hoverinfo="skip",
                    ),
                    row=idx,
                    col=1,
                )

        fig.update_yaxes(title_text=f"Residual ({unit})", row=idx, col=1)

    fig.update_layout(
        title_text="Orbit Determination Measurement Residuals (Prefit / Postfit)",
        template=TEMPLATE,
    )
    fig.update_xaxes(matches="x")
    return watermark(fig, path)


def uncertainty(
    df: pl.DataFrame, sigmas: float = 3.0, path: str | None = None
) -> go.Figure:
    df = convert_units(df)
    sigma_fig = make_subplots(
        rows=2,
        cols=1,
        subplot_titles=[
            f"RIC {sigmas}-sigma position uncertainty (m)",
            f"RIC {sigmas}-sigma velocity uncertainty (m/s)",
        ],
        vertical_spacing=0.1,
    )

    df = df.with_columns(
        [
            (sigmas * pl.col(f"Sigma {axis} (RIC) (m)")).alias(
                f"{sigmas} Sigma {axis} (RIC) (m)"
            )
            for axis in ["X", "Y", "Z"]
        ]
        + [
            (sigmas * pl.col(f"Sigma V{axis.lower()} (RIC) (m/s)")).alias(
                f"{sigmas} Sigma V{axis.lower()} (RIC) (m/s)"
            )
            for axis in ["x", "y", "z"]
        ]
    )

    for trace in px.line(
        df,
        x="Epoch (UTC)",
        y=[f"{sigmas} Sigma {axis} (RIC) (m)" for axis in ["X", "Y", "Z"]],
        template=TEMPLATE,
    ).data:
        sigma_fig.add_trace(trace, row=1, col=1)

    for trace in px.line(
        df,
        x="Epoch (UTC)",
        y=[f"{sigmas} Sigma V{axis.lower()} (RIC) (m/s)" for axis in ["x", "y", "z"]],
        template=TEMPLATE,
    ).data:
        sigma_fig.add_trace(trace, row=2, col=1)

    return watermark(sigma_fig, path)


def od_dashboard(df: pl.DataFrame, path: str | None = None) -> list[go.Figure]:
    """
    Orbit determination dashboard separating scalar measurement streams by observable type,
    ground station attribution, non-occluding station ACF curves, normalized density histograms
    with theoretical N(0,1) envelopes, and station-stratified Q-Q distributions.
    """
    df = convert_units(df)
    msr_types = [
        msr for msr in all_msr_types if f"Measurement noise: {msr}" in df.columns
    ]
    whitened_cols = [c for c in df.columns if "Whitened residual" in c]
    trackers, tracker_colors = _get_tracker_color_map(df)

    # Determine pairing between observable types and whitened columns.
    # Sequential scalar processing uses 'Whitened residual #0' for all observables,
    # requiring dataset partitioning by observable presence (I assume scalar row-wise updates
    # here based on standard Nyx filter output; vectorized runs provide 1:1 column cardinality).
    job_configs = []
    if len(whitened_cols) == 1 and len(msr_types) > 1:
        for msr in msr_types:
            job_configs.append((msr, whitened_cols[0], msr))
    elif len(whitened_cols) == len(msr_types):
        for msr, w_col in zip(msr_types, whitened_cols):
            job_configs.append((msr, w_col, msr))
    else:
        for w_col in whitened_cols:
            job_configs.append((w_col, w_col, None))

    plots = []

    for display_title, w_col, msr_type in job_configs:
        # Slice DataFrame if observable-specific filtering is needed
        if msr_type is not None:
            active_df = df.filter(pl.col(f"Prefit residual: {msr_type}").is_not_null())
        else:
            active_df = df

        if active_df.height == 0:
            continue

        fig = make_subplots(
            rows=4,
            cols=2,
            shared_xaxes=False,
            vertical_spacing=0.08,
            horizontal_spacing=0.08,
            specs=[
                [{"colspan": 2}, None],
                [{"colspan": 2}, None],
                [{"colspan": 2}, None],
                [{}, {}],
            ],
            subplot_titles=(
                f"Residual Rejection Status & 3σ Gate ({display_title})",
                f"Whitened Residuals by Ground Station ({display_title})",
                f"Autocorrelation by Ground Station ({display_title})",
                f"Accepted Residuals Density vs N(0,1) ({display_title})",
                f"Normal Q-Q Distribution by Station ({display_title})",
            ),
        )

        df_acc = active_df.filter(~pl.col("Residual Rejected"))
        df_rej = active_df.filter(pl.col("Residual Rejected"))

        # Row 1: Timeline with explicit Rejection demarcation
        for trk in trackers:
            trk_acc = df_acc.filter(pl.col("Tracker") == trk)
            color = tracker_colors.get(trk, "#4C72B0")
            if trk_acc.height > 0:
                fig.add_trace(
                    go.Scatter(
                        x=trk_acc["Epoch (UTC)"],
                        y=trk_acc[w_col],
                        mode="markers",
                        name=f"{trk} (Accepted)",
                        legendgroup=trk,
                        marker=dict(color=color, symbol="circle", size=5, opacity=0.75),
                        hovertemplate=f"<b>{trk}</b> (Acc)<br>Epoch: %{{x}}<br>Norm Resid: %{{y:.3f}}σ<extra></extra>",
                        showlegend=False,
                    ),
                    row=1,
                    col=1,
                )

        if df_rej.height > 0:
            fig.add_trace(
                go.Scatter(
                    x=df_rej["Epoch (UTC)"],
                    y=df_rej[w_col],
                    mode="markers",
                    name="Rejected",
                    legendgroup="Rejected",
                    marker=dict(
                        color="#C44E52", symbol="x", size=7, line=dict(width=1.5)
                    ),
                    hovertemplate="<b>REJECTED</b><br>Station: %{text}<br>Epoch: %{x}<br>Norm Resid: %{y:.3f}σ<extra></extra>",
                    text=df_rej["Tracker"].to_list(),
                    showlegend=True,
                ),
                row=1,
                col=1,
            )

        # Draw 3-sigma innovation edit boundaries
        for bound_val, bound_name in [(3.0, "+3σ Limit"), (-3.0, "-3σ Limit")]:
            fig.add_trace(
                go.Scatter(
                    x=[active_df["Epoch (UTC)"].min(), active_df["Epoch (UTC)"].max()],
                    y=[bound_val, bound_val],
                    mode="lines",
                    name=bound_name,
                    line=dict(color="rgba(196, 78, 82, 0.7)", dash="dot", width=1.5),
                    hoverinfo="skip",
                    showlegend=False,
                ),
                row=1,
                col=1,
            )

        # Row 2: Station timeline with physical observable context in hover
        for trk in trackers:
            trk_df = active_df.filter(pl.col("Tracker") == trk)
            if trk_df.height == 0:
                continue

            color = tracker_colors.get(trk, "#4C72B0")
            hover_text = []
            for row in trk_df.iter_rows(named=True):
                txt = f"<b>{trk}</b><br>Epoch: {row['Epoch (UTC)']}<br>Norm: {row[w_col]:.3f}σ"
                if msr_type and f"Postfit residual: {msr_type}" in row:
                    unit_str = msr_type.split()[-1][1:-1]
                    txt += f"<br>Postfit: {row[f'Postfit residual: {msr_type}']:.4f} {unit_str}"
                txt += f"<br>Status: {'REJECTED' if row['Residual Rejected'] else 'Accepted'}"
                hover_text.append(txt)

            fig.add_trace(
                go.Scatter(
                    x=trk_df["Epoch (UTC)"],
                    y=trk_df[w_col],
                    mode="markers",
                    name=trk,
                    legendgroup=trk,
                    marker=dict(color=color, symbol="circle", size=5.5),
                    hoverinfo="text",
                    text=hover_text,
                    showlegend=True,
                ),
                row=2,
                col=1,
            )

        # Row 3: Autocorrelation by Tracker (Lines + Markers to eliminate bar occlusion)
        for trk in trackers:
            trk_acc = df_acc.filter(pl.col("Tracker") == trk).sort("Epoch (UTC)")
            x_series = trk_acc[w_col].drop_nulls().to_numpy()
            if len(x_series) < 5:
                continue

            color = tracker_colors.get(trk, "#4C72B0")
            max_lag = min(30, len(x_series) - 2)
            rho = autocorr(x_series, max_lag=max_lag)
            lags = np.arange(len(rho))

            fig.add_trace(
                go.Scatter(
                    x=lags,
                    y=rho,
                    mode="lines+markers",
                    name=f"{trk} ACF",
                    legendgroup=trk,
                    line=dict(color=color, width=1.5),
                    marker=dict(color=color, size=4),
                    hovertemplate=f"<b>{trk}</b> Lag %{{x}}: ρ = %{{y:.3f}}<extra></extra>",
                    showlegend=False,
                ),
                row=3,
                col=1,
            )

            # 95% Bartlett confidence limits
            bound = 1.96 / np.sqrt(len(x_series))
            for b_sign in [1.0, -1.0]:
                fig.add_trace(
                    go.Scatter(
                        x=[0, max_lag],
                        y=[b_sign * bound, b_sign * bound],
                        mode="lines",
                        line=dict(color=color, dash="dot", width=0.8),
                        hoverinfo="skip",
                        showlegend=False,
                    ),
                    row=3,
                    col=1,
                )

        # Row 4, Col 1: Histogram normalized to Probability Density with N(0,1) PDF
        for trk in trackers:
            trk_acc = df_acc.filter(pl.col("Tracker") == trk)
            sample = trk_acc[w_col].drop_nulls().to_numpy()
            if len(sample) < 3:
                continue

            color = tracker_colors.get(trk, "#4C72B0")
            mu, std = np.mean(sample), np.std(sample)

            fig.add_trace(
                go.Histogram(
                    x=sample,
                    histnorm="probability density",
                    name=f"{trk} (μ={mu:.2f}, σ={std:.2f})",
                    legendgroup=trk,
                    marker=dict(color=color),
                    opacity=0.45,
                    showlegend=False,
                ),
                row=4,
                col=1,
            )

        # Standard normal reference PDF curve
        z_grid = np.linspace(-4.0, 4.0, 200)
        norm_pdf = stats.norm.pdf(z_grid, loc=0.0, scale=1.0)
        fig.add_trace(
            go.Scatter(
                x=z_grid,
                y=norm_pdf,
                mode="lines",
                name="N(0,1) Ideal",
                line=dict(color="#1A1A1A", dash="dash", width=1.8),
                hoverinfo="skip",
                showlegend=True,
            ),
            row=4,
            col=1,
        )

        # Row 4, Col 2: Q-Q Plot Stratified by Station
        for trk in trackers:
            trk_acc = df_acc.filter(pl.col("Tracker") == trk)
            sample = trk_acc[w_col].drop_nulls().to_numpy()
            if len(sample) < 3:
                continue

            color = tracker_colors.get(trk, "#4C72B0")
            (osm, osr), (slope, intercept, _) = stats.probplot(
                sample, dist="norm", fit=True
            )

            fig.add_trace(
                go.Scatter(
                    x=osm,
                    y=osr,
                    mode="markers",
                    name=f"{trk} QQ",
                    legendgroup=trk,
                    marker=dict(color=color, size=4.5),
                    hovertemplate=f"<b>{trk}</b> Quantiles<br>Theoretical: %{{x:.2f}}<br>Sample: %{{y:.2f}}<extra></extra>",
                    showlegend=False,
                ),
                row=4,
                col=2,
            )

        # Expected N(0,1) 45-degree reference line
        fig.add_trace(
            go.Scatter(
                x=[-3.5, 3.5],
                y=[-3.5, 3.5],
                mode="lines",
                name="Expected N(0,1)",
                line=dict(color="#C44E52", dash="dash", width=1.5),
                hoverinfo="skip",
                showlegend=False,
            ),
            row=4,
            col=2,
        )

        # Axes updates
        fig.update_yaxes(title_text="Normalized (σ)", row=1, col=1)
        fig.update_yaxes(title_text="Normalized (σ)", row=2, col=1)
        fig.update_yaxes(title_text="Autocorrelation (ρ)", row=3, col=1)
        fig.update_xaxes(title_text="Lag (samples)", row=3, col=1)
        fig.update_yaxes(title_text="Probability Density", row=4, col=1)
        fig.update_xaxes(title_text="Normalized Innovation (σ)", row=4, col=1)
        fig.update_yaxes(title_text="Sample Quantiles (σ)", row=4, col=2)
        fig.update_xaxes(title_text="Theoretical Quantiles", row=4, col=2)

        fig.update_layout(
            title_text=f"Orbit Determination Residual Analysis: {display_title}",
            barmode="overlay",
            template=TEMPLATE,
        )
        plots.append(watermark(fig, path))

    return plots


def cr_cd(df: pl.DataFrame, path: str | None = None) -> go.Figure | None:
    plots_to_make = []
    for col in optional_est_params:
        if df[col].max() != df[col].min():
            sigma_col = next(
                c for c in df.columns if col in c.lower() and "sigma" in c.lower()
            )
            plots_to_make.append((col.capitalize(), col, sigma_col))

    if not plots_to_make:
        return None

    fig = make_subplots(
        rows=len(plots_to_make),
        cols=1,
        subplot_titles=[p[0] for p in plots_to_make],
        vertical_spacing=0.1,
    )

    legend_added = False
    for idx, (title, val_col, sigma_col) in enumerate(plots_to_make, start=1):
        fig.add_trace(
            go.Scatter(
                x=df["Epoch (UTC)"],
                y=df[val_col],
                mode="lines+markers",
                name=title,
                legendgroup=title,
                marker=dict(color="#4C72B0" if "cr" in title.lower() else "#55A868"),
                showlegend=True,
            ),
            row=idx,
            col=1,
        )

        df = df.with_columns(
            [
                (pl.col(val_col) + 3.0 * pl.col(sigma_col)).alias(f"{title} +3-Sigma"),
                (pl.col(val_col) - 3.0 * pl.col(sigma_col)).alias(f"{title} -3-Sigma"),
            ]
        )
        for bound in [f"{title} +3-Sigma", f"{title} -3-Sigma"]:
            fig.add_trace(
                go.Scatter(
                    x=df["Epoch (UTC)"],
                    y=df[bound],
                    mode="lines",
                    name="3-Sigma bounds",
                    line=dict(color="black", dash="dash"),
                    legendgroup="3-Sigma bounds",
                    connectgaps=True,
                    showlegend=(not legend_added),
                ),
                row=idx,
                col=1,
            )
            legend_added = True

        fig.update_yaxes(title_text="Value (unitless)", row=idx, col=1)

    fig.update_layout(
        title_text=" ".join([x[0] for x in plots_to_make]),
        template=TEMPLATE,
    )
    fig.update_xaxes(matches="x")
    return watermark(fig, path)


def kalman_gains(df: pl.DataFrame, path: str | None = None) -> go.Figure | None:
    df = convert_units(df)
    gain_columns = [c for c in df.columns if "Gain" in c]
    if len(df[gain_columns].drop_nulls()) > 0:
        fig = px.scatter(df, x="Epoch (UTC)", y=gain_columns, template=TEMPLATE)
        return watermark(fig, path)
    return None


def filter_smoother_ratios(
    df: pl.DataFrame, path: str | None = None
) -> go.Figure | None:
    df = convert_units(df)
    gain_columns = [c for c in df.columns if "Gain" in c]
    fs_ratio_columns = [c for c in df.columns if "Filter-smoother ratio" in c]
    if (
        len(df[gain_columns].drop_nulls()) == 0
        and len(df[fs_ratio_columns].drop_nulls()) > 0
    ):
        fig = px.scatter(df, x="Epoch (UTC)", y=fs_ratio_columns, template=TEMPLATE)
        return watermark(fig, path)
    return None


def orbital_element_uncertainty(
    df: pl.DataFrame, sigmas: float = 3.0, path: str | None = None
) -> list[go.Figure]:
    columns = [
        "SemiMajorAxis (km)",
        "Eccentricity (unitless)",
        "Inclination (deg)",
        "RAAN (deg)",
        "AoP (deg)",
        "TrueAnomaly (deg)",
        "AoL (deg)",
        "TrueLongitude (deg)",
    ]

    plots = []
    for sigma in [False, True]:
        subplot_titles = (
            [f"{sigmas}-Sigma {col}" for col in columns] if sigma else columns
        )
        fig = make_subplots(
            rows=4,
            cols=2,
            subplot_titles=subplot_titles,
            shared_xaxes=True,
            vertical_spacing=0.1,
        )

        row_i, col_i = 0, 0
        for col in columns:
            if sigma:
                col_name = f"Sigma {col}"
                if col_name not in df.columns:
                    raise ValueError(
                        "Provided dataframe does not contain covariance sigma columns."
                    )
                y = df[col_name] * sigmas
                name = f"{sigmas}-{col}"
            else:
                y = df[col]
                name = col

            fig.add_trace(
                go.Scattergl(x=df["Epoch (UTC)"], y=y, name=name),
                row=row_i + 1,
                col=col_i + 1,
            )

            col_i = (col_i + 1) % 2
            if col_i == 0:
                row_i = (row_i + 1) % 4

        fig.update_layout(template=TEMPLATE)
        plots.append(watermark(fig, path))
    return plots
