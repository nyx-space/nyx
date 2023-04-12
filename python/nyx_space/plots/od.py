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

from .utils import plot_with_error, plot_line

import pandas as pd


def plot_estimates(
    dfs,
    title,
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
            raise ValueError(f"Could not find time column {time_col_name}")
        else:
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
            raise ValueError(f"Could not find time column {time_col_name}")
        else:
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
            raise ValueError(f"Could not find time column {time_col_name}")
        else:
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


# TODO:
# - Add a function to plot the residuals
# - Add a function to plot the measurements on top of any other plotted data


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
