import plotly.graph_objects as go

from .utils import plot_with_error, plot_line

import pandas as pd


def plot_estimates(
    df, title, time_col_name="Epoch:GregorianUtc", html_out=None, copyright=None
):
    """
    Plots the estimates from an orbit determination solution
    """

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

    fig = go.Figure()
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:X (km)",
        covar["cx_x"],
        "blue",
        x_title,
        "Position estimate (km)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:Y (km)",
        covar["cy_y"],
        "green",
        x_title,
        "Position estimate (km)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:Z (km)",
        covar["cz_z"],
        "orange",
        x_title,
        "Position estimate (km)",
        title,
        copyright,
    )
    fig.show()
    if html_out:
        html_out = html_out.replace(".html", "_{}.html")
        this_output = html_out.format("pos_est")
        with open(this_output, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {this_output}")

    fig = go.Figure()
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:VX (km/s)",
        covar["cx_dot_x_dot"],
        "blue",
        x_title,
        "Velocity estimate (km)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:VY (km/s)",
        covar["cy_dot_y_dot"],
        "green",
        x_title,
        "Velocity estimate (km)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:VZ (km/s)",
        covar["cz_dot_z_dot"],
        "orange",
        x_title,
        "Velocity estimate (km)",
        title,
        copyright,
    )

    fig.show()

    if html_out:
        html_out = html_out.replace(".html", "_{}.html")
        this_output = html_out.format("vel_est")
        with open(this_output, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {this_output}")


def plot_covar(
    df, title, time_col_name="Epoch:GregorianUtc", html_out=None, copyright=None
):
    """
    Plot only the covariance from an orbit determination solution
    """

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

    fig = go.Figure()
    plot_line(
        fig,
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
        fig,
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
        fig,
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
    fig.update_layout(yaxis_range=[-0.1 * max_cov_y, 1.5 * max_cov_y])
    fig.show()
    if html_out:
        html_out = html_out.replace(".html", "_{}.html")
        this_output = html_out.format("pos_cov")
        with open(this_output, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {this_output}")

    fig = go.Figure()
    plot_line(
        fig,
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
        fig,
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
        fig,
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

    fig.update_layout(yaxis_range=[-0.1 * max_cov_y, 1.5 * max_cov_y])
    fig.show()
    if html_out:
        html_out = html_out.replace(".html", "_{}.html")
        this_output = html_out.format("vel_cov")
        with open(this_output, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {this_output}")


def plot_state_deviation(
    df, title, time_col_name="Epoch:GregorianUtc", html_out=None, copyright=None
):
    """
    Plot the state deviation from an orbit determination solution
    """

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

    fig = go.Figure()
    plot_with_error(
        fig,
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
        fig,
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
        fig,
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
    fig.update_layout(yaxis_range=[-1.5 * max_cov_y, 1.5 * max_cov_y])
    fig.show()
    if html_out:
        html_out = html_out.replace(".html", "_{}.html")
        this_output = html_out.format("pos_dev")
        with open(this_output, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {this_output}")

    fig = go.Figure()
    plot_with_error(
        fig,
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
        fig,
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
        fig,
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

    fig.update_layout(yaxis_range=[-1.5 * max_cov_y, 1.5 * max_cov_y])
    fig.show()
    if html_out:
        html_out = html_out.replace(".html", "_{}.html")
        this_output = html_out.format("vel_dev")
        with open(this_output, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {this_output}")
