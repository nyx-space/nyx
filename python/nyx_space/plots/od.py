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
import plotly.express as px
from datetime import datetime

from .utils import plot_with_error, plot_line, finalize_plot, colors

import re
from slugify import slugify

import pandas as pd

import numpy as np
from scipy.stats import norm

from nyx_space.time import Epoch

def plot_estimates(
    dfs,
    title,
    cov_frame="RIC",
    cov_fmt="sqrt",
    cov_sigma=3.0,
    msr_df=None,
    time_col_name="Epoch:Gregorian UTC",
    ref_traj=None,
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
        cov_frame (str): Frame in which to plot the covariance, defaults to RIC frame. File also contains the integration frame.
        cov_fmt (str): Defines how to plot the covariance, defaults to "sqrt" (applies np.sqrt to the column), which will convert the covariance from an expectation to an uncertainty.
        cov_sigma (float): Defines the number fo sigmas of uncertainty or covariance to plot, uses [68, 95, 99.7% rule](https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule) for 1.0, 2.0, and 3.0 of sigma.
        msr_df (pandas.DataFrame): The data frame containing the measurements
        time_col_name (str): The name of the time column
        ref_trajs (pandas.DataFrame): If provided, the difference between the reference trajectory and the estimates will be plotted instead of the estimates alone (or list thereof).
        html_out (str): The name of the HTML file to save the plot to
        copyright (str): The copyright to display on the plot
        pos_fig (plotly.graph_objects.Figure): The figure to plot the position estimates on
        vel_fig (plotly.graph_objects.Figure): The figure to plot the velocity estimates on
        show (bool): Whether to show the plot. If set to false, the figure will be returned.
    """

    if not isinstance(dfs, list):
        dfs = [dfs]

    if ref_traj is not None and not isinstance(ref_traj, list):
        ref_traj = [ref_traj]

    if pos_fig is None:
        pos_fig = go.Figure()

    if vel_fig is None:
        vel_fig = go.Figure()

    for num, df in enumerate(dfs):
        try:
            orig_tim_col = df[time_col_name]
        except KeyError:
            # Find the time column
            try:
                col_name = [x for x in df.columns if x.startswith("Epoch")][0]
            except IndexError:
                raise KeyError("Could not find any Epoch column")
            print(f"Could not find time column {time_col_name}, using `{col_name}`")
            orig_tim_col = df[col_name]

        # Build a Python datetime column
        pd_ok_epochs = []
        for epoch in orig_tim_col:
            epoch = epoch.replace("UTC", "").strip()
            if "." not in epoch:
                epoch += ".0"
            pd_ok_epochs += [datetime.fromisoformat(str(epoch).replace("UTC", "").strip())]
        time_col = pd.Series(pd_ok_epochs)
        x_title = "Epoch {}".format(time_col_name[-3:])

        # Check that the requested covariance frame exists
        frames = set(
            [
                re.search(r"\((.*?)\)", c).group(1)
                for c in df.columns
                if c.startswith("Covariance")
            ]
        )
        if cov_frame not in frames:
            raise ValueError(
                f"Covariance frame `{cov_frame}` not in one of the available frames from dataframe: {frames}"
            )
        covar = {}

        # Reference the covariance frames
        for covar_var, covar_col in {
            "cx_x": "Covariance XX",
            "cy_y": "Covariance YY",
            "cz_z": "Covariance ZZ",
            "cx_dot_x_dot": "Covariance VxVx",
            "cy_dot_y_dot": "Covariance VyVy",
            "cz_dot_z_dot": "Covariance VzVz",
        }.items():
            # Create a new column with the transformed covariance. (e.g. "Covariance VzVz (RIC) 1.0-sigma sqrt")
            if cov_fmt is None:
                cov_col_name = f"{covar_col} ({cov_frame}) {cov_sigma}-sigma"
                # No transformation here
                df[cov_col_name] = df[f"{covar_col} ({cov_frame})"]
            else:
                cov_col_name = f"{covar_col} ({cov_frame}) {cov_sigma}-sigma {cov_fmt}"
                # Transform the current column
                df[cov_col_name] = eval(f"np.{cov_fmt}")(
                    df[f"{covar_col} ({cov_frame})"]
                )
            covar[f"{covar_var}"] = cov_col_name

        plt_df = df
        plt_columns = [
            "x (km)",
            "y (km)",
            "z (km)",
            "vx (km/s)",
            "vy (km/s)",
            "vz (km/s)",
        ]
        plt_modifier = ""

        if ref_traj is not None:
            # Merge both dataframes, keeping all of the values
            merged_df = pd.merge(
                df,
                ref_traj[num],
                on=time_col_name,
                how="outer",
                suffixes=("_od", "_ref_traj"),
            )
            if merged_df.shape[0] > max(df.shape[0], ref_traj[num].shape[0]):
                print(
                    "Reference epochs differ between dataframes, performing a merge by closest key on TAI epoch instead"
                )
                merged_df = pd.merge_asof(
                    df, ref_traj[num], on="Epoch:TAI (s)", suffixes=("_od", "_ref_traj")
                )

            # Add the difference to reference to the dataframe
            for coord in plt_columns:
                merged_df[f"delta {coord}"] = (
                    merged_df[f"{coord}_od"] - merged_df[f"{coord}_ref_traj"]
                )
            # Set the plotted dataframe to this one
            plt_df = merged_df
            plt_columns = [f"delta {c}" for c in plt_columns]
            plt_modifier = " error"

        plot_with_error(
            pos_fig,
            plt_df,
            time_col,
            plt_columns[0],
            covar["cx_x"],
            "blue",
            x_title,
            f"Position estimate {plt_modifier}(km)",
            f"Position estimates: {title}",
            copyright,
        )
        plot_with_error(
            pos_fig,
            plt_df,
            time_col,
            plt_columns[1],
            covar["cy_y"],
            "green",
            x_title,
            f"Position estimate {plt_modifier}(km)",
            f"Position estimates: {title}",
            copyright,
        )
        plot_with_error(
            pos_fig,
            plt_df,
            time_col,
            plt_columns[2],
            covar["cz_z"],
            "orange",
            x_title,
            f"Position estimate {plt_modifier}(km)",
            f"Position estimates: {title}",
            copyright,
        )

        plot_with_error(
            vel_fig,
            plt_df,
            time_col,
            plt_columns[3],
            covar["cx_dot_x_dot"],
            "blue",
            x_title,
            f"Velocity estimate {plt_modifier}(km/s)",
            f"Velocity estimates: {title}",
            copyright,
        )
        plot_with_error(
            vel_fig,
            plt_df,
            time_col,
            plt_columns[4],
            covar["cy_dot_y_dot"],
            "green",
            x_title,
            f"Velocity estimate {plt_modifier}(km/s)",
            f"Velocity estimates: {title}",
            copyright,
        )
        plot_with_error(
            vel_fig,
            plt_df,
            time_col,
            plt_columns[5],
            covar["cz_dot_z_dot"],
            "orange",
            x_title,
            f"Velocity estimate {plt_modifier}(km/s)",
            f"Velocity estimates: {title}",
            copyright,
        )

    if msr_df is not None:
        # Plot the measurements on both plots
        pos_fig = overlay_measurements(
            pos_fig, msr_df, title, time_col_name, show=False
        )

        vel_fig = overlay_measurements(
            vel_fig, msr_df, title, time_col_name, show=False
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
    cov_frame="RIC",
    cov_fmt="sqrt",
    cov_sigma=3.0,
    msr_df=None,
    time_col_name="Epoch:Gregorian UTC",
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
        cov_frame (str): Frame in which to plot the covariance, defaults to RIC frame. File also contains the integration frame.
        cov_fmt (str): Defines how to plot the covariance, defaults to "sqrt" (applies np.sqrt to the column), which will convert the covariance from an expectation to an uncertainty.
        cov_sigma (float): Defines the number fo sigmas of uncertainty or covariance to plot, uses [68, 95, 99.7% rule](https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule) for 1.0, 2.0, and 3.0 of sigma.
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
                raise KeyError("Could not find any Epoch column")
            print(f"Could not find time column {time_col_name}, using `{col_name}`")
            orig_tim_col = df[col_name]

        # Build a Python datetime column
        pd_ok_epochs = []
        for epoch in orig_tim_col:
            epoch = epoch.replace("UTC", "").strip()
            if "." not in epoch:
                epoch += ".0"
            pd_ok_epochs += [datetime.fromisoformat(str(epoch).replace("UTC", "").strip())]
        time_col = pd.Series(pd_ok_epochs)
        x_title = "Epoch {}".format(time_col_name[-3:])

        # Check that the requested covariance frame exists
        frames = set(
            [
                re.search(r"\((.*?)\)", c).group(1)
                for c in df.columns
                if c.startswith("Covariance")
            ]
        )
        if cov_frame not in frames:
            raise ValueError(
                f"Covariance frame `{cov_frame}` not in one of the available frames from dataframe: {frames}"
            )
        covar = {}

        # Reference the covariance frames
        for covar_var, covar_col in {
            "cx_x": "Covariance XX",
            "cy_y": "Covariance YY",
            "cz_z": "Covariance ZZ",
            "cx_dot_x_dot": "Covariance VxVx",
            "cy_dot_y_dot": "Covariance VyVy",
            "cz_dot_z_dot": "Covariance VzVz",
        }.items():
            # Create a new column with the transformed covariance. (e.g. "Covariance VzVz (RIC) 1.0-sigma sqrt")
            cov_col_name = f"{covar_col} ({cov_frame}) {cov_sigma}-sigma {cov_fmt}"
            # Transform the current column
            df[cov_col_name] = eval(f"np.{cov_fmt}")(df[f"{covar_col} ({cov_frame})"])
            covar[f"{covar_var}"] = cov_col_name

        plot_line(
            pos_fig,
            df,
            time_col,
            covar["cx_x"],
            "blue",
            x_title,
            f"Position covariance {cov_sigma}-sigma (km)",
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
            f"Position covariance {cov_sigma}-sigma (km)",
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
            f"Position covariance {cov_sigma}-sigma (km)",
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
            f"Velocity covariance {cov_sigma}-sigma (km/s)",
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
            f"Velocity covariance {cov_sigma}-sigma (km/s)",
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
            f"Velocity covariance {cov_sigma}-sigma (km/s)",
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
        pos_fig = overlay_measurements(
            pos_fig, msr_df, title, time_col_name, show=False
        )

        vel_fig = overlay_measurements(
            vel_fig, msr_df, title, time_col_name, show=False
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


def overlay_measurements(
    fig,
    dfs,
    title,
    time_col_name="Epoch:Gregorian UTC",
    html_out=None,
    copyright=None,
    show=True,
):
    """
    Given a plotly figure, overlay the measurements as shaded regions on top of the existing plot.
    For a plot of measurements only, use `plot_measurements`.
    """
    if not isinstance(dfs, list):
        dfs = [dfs]

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
                raise KeyError("Could not find any Epoch column")
            print(f"Could not find time column {time_col_name}, using `{col_name}`")
            orig_tim_col = df[col_name]

        # Build a Python datetime column
        pd_ok_epochs = []
        for epoch in orig_tim_col:
            epoch = epoch.replace("UTC", "").strip()
            if "." not in epoch:
                epoch += ".0"
            pd_ok_epochs += [datetime.fromisoformat(str(epoch).replace("UTC", "").strip())]
        time_col = pd.Series(pd_ok_epochs)
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

        finalize_plot(fig, title, x_title, None, copyright)

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
    kind="Prefit",
    time_col_name="Epoch:Gregorian UTC",
    msr_df=None,
    copyright=None,
    html_out=None,
    show=True,
):
    """
    Plot of residuals, with 3-σ lines. Returns a tuple of the plots if show=False.
    """

    try:
        orig_tim_col = df[time_col_name]
    except KeyError:
        # Find the time column
        try:
            col_name = [x for x in df.columns if x.startswith("Epoch")][0]
        except IndexError:
            raise KeyError("Could not find any Epoch column")
        print(f"Could not find time column {time_col_name}, using `{col_name}`")
        orig_tim_col = df[col_name]

    # Build a Python datetime column
    pd_ok_epochs = []
    for epoch in orig_tim_col:
        epoch = epoch.replace("UTC", "").strip()
        if "." not in epoch:
            epoch += ".0"
        pd_ok_epochs += [datetime.fromisoformat(str(epoch).replace("UTC", "").strip())]
        time_col = pd.Series(pd_ok_epochs)
    x_title = "Epoch {}".format(time_col_name[-3:])

    plt_any = False

    rtn_plots = []

    for col in df.columns:
        if col.startswith(kind):
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
                y=mean + 3 * std,
                line_dash="dash",
                line_color=three_sig_color,
                annotation_text="+ 3σ",
            )
            fig.add_hline(
                y=mean - 3 * std,
                line_dash="dash",
                line_color=three_sig_color,
                annotation_text="- 3σ",
            )
            # Add the 1-σ lines
            one_sig_color = colors["bright_green"]
            one_sig_color = f"rgb({int(one_sig_color[0])}, {int(one_sig_color[1])}, {int(one_sig_color[2])})"
            fig.add_hline(
                y=mean + std,
                line_dash="dot",
                line_color=one_sig_color,
                annotation_text="+ 1σ",
            )
            fig.add_hline(
                y=mean - std,
                line_dash="dot",
                line_color=one_sig_color,
                annotation_text="- 1σ",
            )

            if msr_df is not None:
                # Plot the measurements on both plots
                fig = overlay_measurements(
                    fig, msr_df, title, time_col_name, show=False
                )

            finalize_plot(
                fig, title=f"{title} {col}", xtitle=x_title, copyright=copyright
            )

            plt_any = True

            if html_out:
                this_output = html_out.replace(".html", f"_{slugify(col)}.html")
                with open(this_output, "w") as f:
                    f.write(fig.to_html())
                print(f"Saved HTML to {this_output}")

            if show:
                fig.show()
            else:
                rtn_plots += [fig]

    if not plt_any:
        raise ValueError(f"No columns ending with {kind} found -- nothing plotted")

    if not show:
        return rtn_plots


def plot_residual_histogram(
    df, title, kind="Prefit", copyright=None, html_out=None, show=True
):
    """
    Histogram of residuals
    """

    for col in df.columns:
        if col.startswith(kind):
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

            if html_out:
                this_output = html_out.replace(".html", f"_{slugify(col)}.html")
                with open(this_output, "w") as f:
                    f.write(fig.to_html())
                print(f"Saved HTML to {this_output}")

            if show:
                fig.show()

def plot_measurements(
    df,
    msr_type=None,
    title=None,
    time_col_name="Epoch:Gregorian UTC",
    html_out=None,
    copyright=None,
    show=True,
):
    """
    Plot the provided measurement type, fuzzy matching of the column name, or plot all as a strip
    """

    if title is None:
        # Build a title
        station_names = ", ".join([name for name in df["Tracking device"].unique()])
        start = Epoch(df["Epoch:Gregorian UTC"].iloc[0])
        end = Epoch(df["Epoch:Gregorian UTC"].iloc[-1])
        arc_duration = end.timedelta(start)
        title = f"Measurements from {station_names} spanning {start} to {end} ({arc_duration})"

    try:
        orig_tim_col = df[time_col_name]
    except KeyError:
        # Find the time column
        try:
            col_name = [x for x in df.columns if x.startswith("Epoch")][0]
        except IndexError:
            raise KeyError("Could not find any Epoch column")
        print(f"Could not find time column {time_col_name}, using `{col_name}`")
        orig_tim_col = df[col_name]

    # Build a Python datetime column
    pd_ok_epochs = []
    for epoch in orig_tim_col:
        epoch = epoch.replace("UTC", "").strip()
        if "." not in epoch:
            epoch += ".0"
        pd_ok_epochs += [datetime.fromisoformat(str(epoch).replace("UTC", "").strip())]
    df["time_col"] = pd.Series(pd_ok_epochs)
    x_title = "Epoch {}".format(time_col_name[-3:])

    if msr_type is None:
        fig = px.strip(df, x="time_col", y="Tracking device", color="Tracking device")
        finalize_plot(fig, title, x_title, "All tracking data", copyright)
    else:
        msr_col_name = [col for col in df.columns if msr_type in col.lower()]

        fig = px.scatter(df, x="time_col", y=msr_col_name, color="Tracking device")
        finalize_plot(fig, title, x_title, msr_col_name[0], copyright)

    if html_out:
        with open(html_out, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {html_out}")

    if show:
        fig.show()
    else:
        return fig
