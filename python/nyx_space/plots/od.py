import plotly.graph_objects as go

from .utils import plot_with_error, plot_line, plot_raw_line


def plot_estimates(
    df, title, time_col_name="Epoch:GregorianUtc", html_out=None, copyright=None
):
    """
    Plots the estimates from an orbit determination solution
    """

    try:
        time_col = df[time_col_name]
    except KeyError:
        raise ValueError(f"Could not find time column {time_col_name}")

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
                locals()[covar_var] = f"{covar_var}{frame}"

        if covar_var not in locals():
            raise KeyError(f"Cannot find {covar_var} covariance column")

    fig = go.Figure()
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:X (km)",
        cx_x,
        "blue",
        "Position estimate (km)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:Y (km)",
        cy_y,
        "green",
        "Position estimate (km)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:Z (km)",
        cz_z,
        "orange",
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
        cx_dot_x_dot,
        "blue",
        "Velocity estimate (km)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:VY (km/s)",
        cy_dot_y_dot,
        "green",
        "Velocity estimate (km)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "Estimate:VZ (km/s)",
        cz_dot_z_dot,
        "orange",
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
        time_col = df[time_col_name]
    except KeyError:
        raise ValueError(f"Could not find time column {time_col_name}")

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
                locals()[covar_var] = f"{covar_var}{frame}"

        if covar_var not in locals():
            raise KeyError(f"Cannot find {covar_var} covariance column")

    fig = go.Figure()
    plot_line(
        fig, df, time_col, cx_x, "blue", "Position covariance (km)", title, copyright
    )
    plot_line(
        fig, df, time_col, cy_y, "green", "Position covariance (km)", title, copyright
    )
    plot_line(
        fig, df, time_col, cz_z, "orange", "Position covariance (km)", title, copyright
    )

    # Autoscale
    hwpt = int(len(df) / 2)
    max_cov_y = max(max(df[cx_x][hwpt:]), max(df[cy_y][hwpt:]), max(df[cz_z][hwpt:]))
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
        cx_dot_x_dot,
        "blue",
        "Velocity covariance (km/s)",
        title,
        copyright,
    )
    plot_line(
        fig,
        df,
        time_col,
        cy_dot_y_dot,
        "green",
        "Velocity covariance (km/s)",
        title,
        copyright,
    )
    plot_line(
        fig,
        df,
        time_col,
        cz_dot_z_dot,
        "orange",
        "Velocity covariance (km/s)",
        title,
        copyright,
    )
    # Autoscale
    hwpt = int(len(df) / 2)
    max_cov_y = max(
        max(df[cx_dot_x_dot][hwpt:]),
        max(df[cy_dot_y_dot][hwpt:]),
        max(df[cz_dot_z_dot][hwpt:]),
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
        time_col = df[time_col_name]
    except KeyError:
        raise ValueError(f"Could not find time column {time_col_name}")

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
                locals()[covar_var] = f"{covar_var}{frame}"

        if covar_var not in locals():
            raise KeyError(f"Cannot find {covar_var} covariance column")

    fig = go.Figure()
    plot_with_error(
        fig,
        df,
        time_col,
        "delta_x",
        cx_x,
        "blue",
        "Position deviation (km)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "delta_y",
        cy_y,
        "green",
        "Position deviation (km)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "delta_z",
        cz_z,
        "orange",
        "Position deviation (km)",
        title,
        copyright,
    )
    # Autoscale
    hwpt = int(len(df) / 2)
    max_cov_y = max(max(df[cx_x][hwpt:]), max(df[cy_y][hwpt:]), max(df[cz_z][hwpt:]))
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
        cx_dot_x_dot,
        "blue",
        "Velocity deviation (km/s)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "delta_vy",
        cy_dot_y_dot,
        "green",
        "Velocity deviation (km/s)",
        title,
        copyright,
    )
    plot_with_error(
        fig,
        df,
        time_col,
        "delta_vz",
        cz_dot_z_dot,
        "orange",
        "Velocity deviation (km/s)",
        title,
        copyright,
    )
    # Autoscale
    hwpt = int(len(df) / 2)
    max_cov_y = max(
        max(df[cx_dot_x_dot][hwpt:]),
        max(df[cy_dot_y_dot][hwpt:]),
        max(df[cz_dot_z_dot][hwpt:]),
    )

    fig.update_layout(yaxis_range=[-1.5 * max_cov_y, 1.5 * max_cov_y])
    fig.show()
    if html_out:
        html_out = html_out.replace(".html", "_{}.html")
        this_output = html_out.format("vel_dev")
        with open(this_output, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {this_output}")
