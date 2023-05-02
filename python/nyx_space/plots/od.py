"""
Nyx, blazing fast astrodynamics
Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Affero General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Affero General Public License for more details.

You should have received a copy of the GNU Affero General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

import plotly.graph_objects as go

from .utils import plot_with_error, plot_line, finalize_plot, colors

import pandas as pd

import numpy as np
from scipy.stats import norm
from scipy.special import erfcinv


def plot_estimates(
    dfs,
    title,
    msr_df=None,
    time_col_name="Epoch:GregorianUtc",
    html_out=None,
    copyright=None,
    pos_fig=None,
    vel_fig=None,
    show=True,
):
    """
    Plots the estimates from an orbit determination solution

    Args:
        dfs (pandas.DataFrame): The data frame containing the orbit determination solution (or a list thereof)
        title (str): The title of the plot
        msr_df (pandas.DataFrame): The data frame containing the measurements
        time_col_name (str): The name of the time column
        html_out (str): The name of the HTML file to save the plot to
        copyright (str): The copyright to display on the plot
        pos_fig (plotly.graph_objects.Figure): The figure to plot the position estimates on
        vel_fig (plotly.graph_objects.Figure): The figure to plot the velocity estimates on
        show (bool): Whether to show the plot. If set to false, the figure will be returned.
    """

    if not isinstance(dfs, list):
        dfs = [dfs]

    if pos_fig is None:
        pos_fig = go.Figure()

    if vel_fig is None:
        vel_fig = go.Figure()

    for df in dfs:
        try:
            orig_tim_col = df[time_col_name]
        except KeyError:
            # Find the time column
            try:
                col_name = [x for x in df.columns if x.startswith("Epoch")][0]
            except IndexError:
                raise KeyError("Could not find any time column")
            print(f"Could not find time column {time_col_name}, using `{col_name}`")
            orig_tim_col = df[col_name]

        # Build a Python datetime column
        time_col = pd.to_datetime(orig_tim_col)
        x_title = "Epoch {}".format(time_col_name[-3:])

        covar = {}

        # Reference the covariance frames
        for covar_var in [
            "cx_x",
            "cy_y",
            "cz_z",
            "cx_dot_x_dot",
            "cy_dot_y_dot",
            "cz_dot_z_dot",
        ]:
            possible_frames = ["", ":ric", ":rcn", ":vnc"]
            for frame in possible_frames:
                if f"{covar_var}{frame}" in df:
                    covar[covar_var] = f"{covar_var}{frame}"

            if covar_var not in covar:
                raise KeyError(f"Cannot find {covar_var} covariance column")

        plot_with_error(
            pos_fig,
            df,
            time_col,
            "Estimate:X (km)",
            covar["cx_x"],
            "blue",
            x_title,
            "Position estimate (km)",
            f"Position estimates: {title}",
            copyright,
        )
        plot_with_error(
            pos_fig,
            df,
            time_col,
            "Estimate:Y (km)",
            covar["cy_y"],
            "green",
            x_title,
            "Position estimate (km)",
            f"Position estimates: {title}",
            copyright,
        )
        plot_with_error(
            pos_fig,
            df,
            time_col,
            "Estimate:Z (km)",
            covar["cz_z"],
            "orange",
            x_title,
            "Position estimate (km)",
            f"Position estimates: {title}",
            copyright,
        )

        plot_with_error(
            vel_fig,
            df,
            time_col,
            "Estimate:VX (km/s)",
            covar["cx_dot_x_dot"],
            "blue",
            x_title,
            "Velocity estimate (km)",
            f"Velocity estimates: {title}",
            copyright,
        )
        plot_with_error(
            vel_fig,
            df,
            time_col,
            "Estimate:VY (km/s)",
            covar["cy_dot_y_dot"],
            "green",
            x_title,
            "Velocity estimate (km)",
            f"Velocity estimates: {title}",
            copyright,
        )
        plot_with_error(
            vel_fig,
            df,
            time_col,
            "Estimate:VZ (km/s)",
            covar["cz_dot_z_dot"],
            "orange",
            x_title,
            "Velocity estimate (km)",
            f"Velocity estimates: {title}",
            copyright,
        )

    if msr_df is not None:
        # Plot the measurements on both plots
        pos_fig = plot_measurements(
            msr_df, title, time_col_name, fig=pos_fig, show=False
        )

        vel_fig = plot_measurements(
            msr_df, title, time_col_name, fig=vel_fig, show=False
        )

    if html_out:
        html_out = html_out.replace(".html", "_{}.html")
        this_output = html_out.format("pos_est")
        with open(this_output, "w") as f:
            f.write(pos_fig.to_html())
        print(f"Saved HTML to {this_output}")

        this_output = html_out.format("vel_est")
        with open(this_output, "w") as f:
            f.write(vel_fig.to_html())
        print(f"Saved HTML to {this_output}")

    if show:
        pos_fig.show()
        vel_fig.show()
    else:
        return pos_fig, vel_fig


def plot_covar(
    dfs,
    title,
    msr_df=None,
    time_col_name="Epoch:GregorianUtc",
    html_out=None,
    copyright=None,
    pos_fig=None,
    vel_fig=None,
    show=True,
):
    """
    Plot only the covariance from an orbit determination solution

    Args:
        dfs (pandas.DataFrame): The data frame containing the orbit determination solution (or a list thereof)
        title (str): The title of the plot
        msr_df (pandas.DataFrame): The data frame containing the measurements
        time_col_name (str): The name of the time column
        html_out (str): The name of the HTML file to save the plot to
        copyright (str): The copyright to display on the plot
        pos_fig (plotly.graph_objects.Figure): The figure to plot the position estimates on
        vel_fig (plotly.graph_objects.Figure): The figure to plot the velocity estimates on
        show (bool): Whether to show the plot. If set to false, the figure will be returned.
    """

    if not isinstance(dfs, list):
        dfs = [dfs]

    if pos_fig is None:
        pos_fig = go.Figure()

    if vel_fig is None:
        vel_fig = go.Figure()

    for df in dfs:

        try:
            orig_tim_col = df[time_col_name]
        except KeyError:
            # Find the time column
            try:
                col_name = [x for x in df.columns if x.startswith("Epoch")][0]
            except IndexError:
                raise KeyError("Could not find any time column")
            print(f"Could not find time column {time_col_name}, using `{col_name}`")
            orig_tim_col = df[col_name]

        # Build a Python datetime column
        time_col = pd.to_datetime(orig_tim_col)
        x_title = "Epoch {}".format(time_col_name[-3:])

        covar = {}

        # Reference the covariance frames
        for covar_var in [
            "cx_x",
            "cy_y",
            "cz_z",
            "cx_dot_x_dot",
            "cy_dot_y_dot",
            "cz_dot_z_dot",
        ]:
            possible_frames = ["", ":ric", ":rcn", ":vnc"]
            for frame in possible_frames:
                if f"{covar_var}{frame}" in df:
                    covar[covar_var] = f"{covar_var}{frame}"

            if covar_var not in covar:
                raise KeyError(f"Cannot find {covar_var} covariance column")

        plot_line(
            pos_fig,
            df,
            time_col,
            covar["cx_x"],
            "blue",
            x_title,
            "Position covariance (km)",
            title,
            copyright,
        )
        plot_line(
            pos_fig,
            df,
            time_col,
            covar["cy_y"],
            "green",
            x_title,
            "Position covariance (km)",
            title,
            copyright,
        )
        plot_line(
            pos_fig,
            df,
            time_col,
            covar["cz_z"],
            "orange",
            x_title,
            "Position covariance (km)",
            title,
            copyright,
        )

        # Autoscale
        hwpt = int(len(df) / 2)
        max_cov_y = max(
            max(df[covar["cx_x"]][hwpt:]),
            max(df[covar["cy_y"]][hwpt:]),
            max(df[covar["cz_z"]][hwpt:]),
        )
        pos_fig.update_layout(yaxis_range=[-0.1 * max_cov_y, 1.5 * max_cov_y])

        plot_line(
            vel_fig,
            df,
            time_col,
            covar["cx_dot_x_dot"],
            "blue",
            x_title,
            "Velocity covariance (km/s)",
            title,
            copyright,
        )
        plot_line(
            vel_fig,
            df,
            time_col,
            covar["cy_dot_y_dot"],
            "green",
            x_title,
            "Velocity covariance (km/s)",
            title,
            copyright,
        )
        plot_line(
            vel_fig,
            df,
            time_col,
            covar["cz_dot_z_dot"],
            "orange",
            x_title,
            "Velocity covariance (km/s)",
            title,
            copyright,
        )
        # Autoscale
        hwpt = int(len(df) / 2)
        max_cov_y = max(
            max(df[covar["cx_dot_x_dot"]][hwpt:]),
            max(df[covar["cy_dot_y_dot"]][hwpt:]),
            max(df[covar["cz_dot_z_dot"]][hwpt:]),
        )

        vel_fig.update_layout(yaxis_range=[-0.1 * max_cov_y, 1.5 * max_cov_y])

    if msr_df is not None:
        # Plot the measurements on both plots
        pos_fig = plot_measurements(
            msr_df, title, time_col_name, fig=pos_fig, show=False
        )

        vel_fig = plot_measurements(
            msr_df, title, time_col_name, fig=vel_fig, show=False
        )

    if html_out:

        html_out = html_out.replace(".html", "_{}.html")
        this_output = html_out.format("pos_cov")
        with open(this_output, "w") as f:
            f.write(pos_fig.to_html())
        print(f"Saved HTML to {this_output}")

        this_output = html_out.format("vel_cov")
        with open(this_output, "w") as f:
            f.write(vel_fig.to_html())
        print(f"Saved HTML to {this_output}")

    if show:
        pos_fig.show()
        vel_fig.show()
    else:
        return pos_fig, vel_fig


def plot_state_deviation(
    dfs,
    title,
    msr_df=None,
    traj_err=None,
    time_col_name="Epoch:GregorianUtc",
    html_out=None,
    copyright=None,
    pos_fig=None,
    vel_fig=None,
    show=True,
):
    """
    Plot the state deviation from an orbit determination solution

    Args:
        dfs (pandas.DataFrame): The data frame containing the orbit determination solution (or list thereof)
        title (str): The title of the plot
        msr_df (pandas.DataFrame): The data frame containing the measurements
        traj_err (pandas.DataFrame): The data frame containing the trajectory error compared to truth (or list thereof)
        time_col_name (str): The name of the time column
        html_out (str): The name of the HTML file to save the plot to
        copyright (str): The copyright to display on the plot
        pos_fig (plotly.graph_objects.Figure): The figure to plot the position estimates on
        vel_fig (plotly.graph_objects.Figure): The figure to plot the velocity estimates on
        show (bool): Whether to show the plot. If set to false, the figure will be returned.
    """

    if not isinstance(dfs, list):
        dfs = [dfs]

    if pos_fig is None:
        pos_fig = go.Figure()

    if vel_fig is None:
        vel_fig = go.Figure()

    x_title = "" # Forward declaration, will store the last time column name

    for df in dfs:

        try:
            orig_tim_col = df[time_col_name]
        except KeyError:
            # Find the time column
            try:
                col_name = [x for x in df.columns if x.startswith("Epoch")][0]
            except IndexError:
                raise KeyError("Could not find any time column")
            print(f"Could not find time column {time_col_name}, using `{col_name}`")
            orig_tim_col = df[col_name]

        # Build a Python datetime column
        time_col = pd.to_datetime(orig_tim_col)
        x_title = "Epoch {}".format(time_col_name[-3:])

        # If there is a trajectory error, rename the state deviations with the actual trajectory errors
        if traj_err is not None:
            df["delta_x"] = traj_err["x_err_km"]
            df["delta_y"] = traj_err["y_err_km"]
            df["delta_z"] = traj_err["z_err_km"]

            df["delta_vx"] = traj_err["vx_err_km_s"]
            df["delta_vy"] = traj_err["vy_err_km_s"]
            df["delta_vz"] = traj_err["vz_err_km_s"]
            print("Overwrote the state deviations with the trajectory error")

        covar = {}

        # Reference the covariance frames
        for covar_var in [
            "cx_x",
            "cy_y",
            "cz_z",
            "cx_dot_x_dot",
            "cy_dot_y_dot",
            "cz_dot_z_dot",
        ]:
            possible_frames = ["", ":ric", ":rcn", ":vnc"]
            for frame in possible_frames:
                if f"{covar_var}{frame}" in df:
                    covar[covar_var] = f"{covar_var}{frame}"

            if covar_var not in covar:
                raise KeyError(f"Cannot find {covar_var} covariance column")

        plot_with_error(
            pos_fig,
            df,
            time_col,
            "delta_x",
            covar["cx_x"],
            "blue",
            x_title,
            "Position deviation (km)",
            title,
            copyright,
        )
        plot_with_error(
            pos_fig,
            df,
            time_col,
            "delta_y",
            covar["cy_y"],
            "green",
            x_title,
            "Position deviation (km)",
            title,
            copyright,
        )
        plot_with_error(
            pos_fig,
            df,
            time_col,
            "delta_z",
            covar["cz_z"],
            "orange",
            x_title,
            "Position deviation (km)",
            title,
            copyright,
        )
        # Autoscale
        hwpt = int(len(df) / 2)
        max_cov_y = max(
            max(df[covar["cx_x"]][hwpt:]),
            max(df[covar["cy_y"]][hwpt:]),
            max(df[covar["cz_z"]][hwpt:]),
        )


        pos_fig.update_layout(yaxis_range=[-1.5 * max_cov_y, 1.5 * max_cov_y])

        plot_with_error(
            vel_fig,
            df,
            time_col,
            "delta_vx",
            covar["cx_dot_x_dot"],
            "blue",
            x_title,
            "Velocity deviation (km/s)",
            title,
            copyright,
        )
        plot_with_error(
            vel_fig,
            df,
            time_col,
            "delta_vy",
            covar["cy_dot_y_dot"],
            "green",
            x_title,
            "Velocity deviation (km/s)",
            title,
            copyright,
        )
        plot_with_error(
            vel_fig,
            df,
            time_col,
            "delta_vz",
            covar["cz_dot_z_dot"],
            "orange",
            x_title,
            "Velocity deviation (km/s)",
            title,
            copyright,
        )
        # Autoscale
        hwpt = int(len(df) / 2)
        max_cov_y = max(
            max(df[covar["cx_dot_x_dot"]][hwpt:]),
            max(df[covar["cy_dot_y_dot"]][hwpt:]),
            max(df[covar["cz_dot_z_dot"]][hwpt:]),
        )

        vel_fig.update_layout(yaxis_range=[-1.5 * max_cov_y, 1.5 * max_cov_y])

    if msr_df is not None:
        # Plot the measurements on both plots
        pos_fig = plot_measurements(
            msr_df, title, time_col_name, fig=pos_fig, show=False
        )

        vel_fig = plot_measurements(
            msr_df, title, time_col_name, fig=vel_fig, show=False
        )
    
    finalize_plot(pos_fig, title, x_title, "Position deviation (km)")
    finalize_plot(vel_fig, title, x_title, "Velocity deviation (km/s)")

    if html_out:
        html_out = html_out.replace(".html", "_{}.html")

        this_output = html_out.format("pos_dev")
        with open(this_output, "w") as f:
            f.write(pos_fig.to_html())
        print(f"Saved HTML to {this_output}")

        this_output = html_out.format("vel_dev")
        with open(this_output, "w") as f:
            f.write(vel_fig.to_html())
        print(f"Saved HTML to {this_output}")

    if show:
        pos_fig.show()
        vel_fig.show()
    else:
        return pos_fig, vel_fig


def plot_measurements(
    dfs,
    title,
    time_col_name="Epoch:GregorianUtc",
    html_out=None,
    copyright=None,
    fig=None,
    show=True,
):
    if not isinstance(dfs, list):
        dfs = [dfs]

    if fig is None:
        fig = go.Figure()

    color_values = list(colors.values())

    station_colors = {}

    for df in dfs:

        try:
            orig_tim_col = df[time_col_name]
        except KeyError:
            # Find the time column
            try:
                col_name = [x for x in df.columns if x.startswith("Epoch")][0]
            except IndexError:
                raise KeyError("Could not find any time column")
            print(f"Could not find time column {time_col_name}, using `{col_name}`")
            orig_tim_col = df[col_name]

        # Build a Python datetime column
        time_col = pd.to_datetime(orig_tim_col)
        x_title = "Epoch {}".format(time_col_name[-3:])

        # Diff the epochs of the measurements to find when there is a start and end.

        # Get Epoch TAI diffs of diffs
        diff2 = df["Epoch:TAI (s)"].diff().diff()

        # Lists to hold segments
        station_names = []
        start_indices = []
        end_indices = []

        # Loop through diff2 values
        for i, val in enumerate(diff2):
            # Positive spike indicates end of station segment at the prior index and the start of a new segment
            if val > 0:
                end_idx = i - 1
                end_indices.append(end_idx)
                start_idx = i
                start_indices.append(start_idx)
                station_names.append(df.loc[start_idx, "Tracking device"])
        # Add the start and end data
        start_indices = [0] + start_indices
        end_indices += [len(df) - 1]
        station_names = [df["Tracking device"][0]] + station_names

        # Zip into segments
        station_segments = list(zip(station_names, start_indices, end_indices))

        # Add tracking passes to the plot
        for ii, (name, start_idx, end_idx) in enumerate(station_segments):
            start_epoch = time_col[start_idx]
            end_epoch = time_col[end_idx]

            try:
                color = station_colors[name]
            except KeyError:
                color = color_values[(len(fig.data) + ii) % len(color_values)]
                color = f"rgb({int(color[0])}, {int(color[1])}, {int(color[2])})"
                station_colors[name] = color

            fig.add_vrect(
                x0=start_epoch,
                x1=end_epoch,
                annotation_text=name,
                annotation_position="top left",
                fillcolor=color,
                opacity=0.25,
                line_width=0,
            )

        finalize_plot(fig, title, x_title, copyright, show)

    if html_out:
        with open(html_out, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {html_out}")

    if show:
        fig.show()
    else:
        return fig


def plot_residuals(
    df,
    title,
    kind="prefit",
    time_col_name="Epoch:GregorianUtc",
    msr_df=None,
    copyright=None,
):
    """
    Plot of residuals, with 3-σ lines
    """

    try:
        orig_tim_col = df[time_col_name]
    except KeyError:
        # Find the time column
        try:
            col_name = [x for x in df.columns if x.startswith("Epoch")][0]
        except IndexError:
            raise KeyError("Could not find any time column")
        print(f"Could not find time column {time_col_name}, using `{col_name}`")
        orig_tim_col = df[col_name]

    # Build a Python datetime column
    time_col = pd.to_datetime(orig_tim_col)
    x_title = "Epoch {}".format(time_col_name[-3:])

    plt_any = False

    for col in df.columns:
        if col.endswith(kind):
            fig = go.Figure()
            residuals = df[col]

            # Add scatter trace of residual vs normal quantiles
            fig.add_trace(
                go.Scatter(
                    x=time_col,
                    y=residuals,
                    mode="markers",
                    name=f"Residuals {col}",
                )
            )

            mean = np.mean(residuals)
            std = np.std(residuals)
            # Add the 3-σ lines
            three_sig_color = colors["red_ish"]
            three_sig_color = f"rgb({int(three_sig_color[0])}, {int(three_sig_color[1])}, {int(three_sig_color[2])})"
            fig.add_hline(
                y=mean + 3 * std, line_dash="dash", line_color=three_sig_color
            )
            fig.add_hline(
                y=mean - 3 * std, line_dash="dash", line_color=three_sig_color
            )
            # Add the 1-σ lines
            one_sig_color = colors["bright_green"]
            one_sig_color = f"rgb({int(one_sig_color[0])}, {int(one_sig_color[1])}, {int(one_sig_color[2])})"
            fig.add_hline(y=mean + std, line_dash="dot", line_color=one_sig_color)
            fig.add_hline(y=mean - std, line_dash="dot", line_color=one_sig_color)

            if msr_df is not None:
                # Plot the measurements on both plots
                fig = plot_measurements(
                    msr_df, title, time_col_name, fig=fig, show=False
                )

            finalize_plot(
                fig, title=f"{title} {col}", xtitle=x_title, copyright=copyright
            )
            fig.show()
            plt_any = True

    if not plt_any:
        raise ValueError(f"No columns ending with {kind} found -- nothing plotted")


def plot_residual_histogram(df, title, kind="prefit", copyright=None):
    """
    Histogram of residuals
    """

    for col in df.columns:
        if col.endswith(kind):
            residuals = df[col]
            fig = go.Figure()
            fig.add_trace(
                go.Histogram(
                    x=residuals,
                    name=f"Histogram of Residuals for {col}",
                )
            )

            # Add normal distribution with same mean and std
            mean = np.mean(residuals)
            std = np.std(residuals)
            x = np.linspace(mean - 3 * std, mean + 3 * std, 100)
            pdf = norm.pdf(x, mean, std)
            fig.add_trace(
                go.Scatter(
                    x=x,
                    y=pdf,
                    name=f"Normal Distribution for {col}",
                    line_dash="dash",
                    line_color="red",
                )
            )

            finalize_plot(fig, title, xtitle=None, copyright=copyright)

            fig.show()


def plot_residual_qq(
    df,
    title,
    col="residual_ratio",
    copyright=None,
):
    """
    Plot of residuals, with 3-σ lines
    """

    residuals = df[col]

    # Determine the number of residuals
    # n = len(residuals)
    quantiles = np.arange(0, 1.01, 0.01)

    # Calculate the residuals quantiles using np.quantile:

    residual_quantiles = np.quantile(residuals, quantiles)

    # Calculate the theoretical quantiles for a normal distribution with mean=0 and std=1:

    normal_quantiles = np.sqrt(2) * erfcinv(2 * quantiles)

    # Create figure
    fig = go.Figure()

    # Add scatter trace of residual vs normal quantiles
    fig.add_trace(
        go.Scatter(
            x=normal_quantiles, y=residual_quantiles, mode="markers", name="Quantiles"
        )
    )

    # Add line representing normal distribution
    fig.add_trace(
        go.Scatter(
            x=[min(normal_quantiles), max(normal_quantiles)],
            y=[min(normal_quantiles), max(normal_quantiles)],
            mode="lines",
            name="Normal",
            line=dict(color="red"),
        )
    )

    # Update layout
    fig.update_layout(
        title="Normal Q-Q Plot",
        xaxis_title="Theoretical Quantiles",
        yaxis_title="Residual Quantiles",
    )

    # Display interactive plot
    finalize_plot(fig, title, copyright=copyright)
    fig.show()


def plot_qq(df, title, col="residual_ratio", copyright=None):

    # Your list of residuals
    residuals = df[col]

    # Standardize the residuals
    standardized_residuals = (residuals - np.mean(residuals)) / np.std(residuals)

    # Compute the theoretical normal quantiles (percentiles)
    residuals_sorted = np.sort(standardized_residuals)
    theoretical_quantiles = norm.ppf(
        np.linspace(0.5 / len(residuals), 1 - 0.5 / len(residuals), len(residuals))
    )

    # Create the Q-Q plot using Plotly
    fig = go.Figure()

    # Add the scatter plot of standardized residuals
    fig.add_trace(
        go.Scatter(
            x=theoretical_quantiles,
            y=residuals_sorted,
            mode="markers",
            name="Residuals",
        )
    )

    # Add the 45-degree line representing the normal distribution
    min_value = np.min([theoretical_quantiles.min(), residuals_sorted.min()])
    max_value = np.max([theoretical_quantiles.max(), residuals_sorted.max()])
    fig.add_shape(
        type="line",
        x0=min_value,
        x1=max_value,
        y0=min_value,
        y1=max_value,
        yref="y",
        xref="x",
        line=dict(color="red"),
        name="Normal Distribution",
    )

    fig.update_layout(
        title="Q-Q Plot",
        xaxis_title="Theoretical Normal Quantiles",
        yaxis_title="Sorted Standardized Residuals",
    )
    fig.show()
