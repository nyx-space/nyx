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
from datetime import datetime
import numpy as np

# Color blindness friendly colors, from https://jfly.uni-koeln.de/color/#pallet
# And from https://davidmathlogic.com/colorblind/#%23E1BE6A-%2340B0A6-%237bcb12-%233d9087
colors = {
    "sky_blue": (86, 180, 233),  # #56b4e9
    "purple": (204, 121, 167),  # #cc79a7
    "green": (0, 158, 115),  # #009e73
    "red_ish": (213, 94, 0),  # #d55e00
    "orange": (230, 159, 0),  # #e69f00
    "blue": (0, 114, 178),  # #0072b2
    "yellow": (240, 228, 66),  # #f0e442
    "yellow2": (255, 194, 10),  # #ffc20a
    "blue2": (12, 123, 220),  # #0c7bdc
    "bright_green": (26, 255, 26),  # #1aff1a
    "purple2": (75, 0, 146),  # #4b0092
    "brown": (153, 79, 0),  # #994f00
    "blue3": (0, 108, 209),  # #006cd1
    "light_yellow": (254, 254, 98),  # #fefe62
    "pink": (211, 95, 183),  # #d35fb7
    "black": (0, 0, 0),
    "gray": (231, 231, 231),  # #e7e7e7
    "extra-1": (216, 27, 96),  # #D81B60
    "extra-2": (30, 136, 229),  # #1E88E5
    "extra-4": (0, 77, 64),  # #004D40
    "extra-5": (215, 106, 161),  # #d76aa1
    "extra-6": (46, 18, 240),  # #2e12f0
    "extra-7": (190, 222, 214),  # #beded6
    "extra-8": (27, 179, 109),  # #1bb36d
    "extra-9": (232, 17, 141),  # #e8118d
    "extra-10": (27, 33, 247),  # #1b21f7
    "extra-3": (255, 193, 7),  # #FFC107
}

# Radii of different celestial bodies in km
radii = {
    "Sun": 696342.000000,
    "Mercury": 2439.700000,
    "Venus": 6051.900000,
    "Earth": 6378.136300,
    "Luna": 1738.100000,
    "Mars": 3397.000000,
    "Jupiter": 71492.000000,
    "Saturn": 60268.000000,
    "Uranus": 25559.000000,
    "Neptune": 25269.000000,
}

body_color = {
    "Sun": colors["yellow"],
    "Mercury": colors["extra-7"],
    "Venus": colors["orange"],
    "Earth": colors["sky_blue"],
    "Luna": colors["gray"],
    "Mars": colors["red_ish"],
    "Jupiter": colors["yellow2"],
    "Saturn": colors["extra-3"],
    "Uranus": colors["extra-2"],
    "Neptune": colors["extra-10"],
}


def _add_watermark(who):
    nyx_tpl = go.layout.Template()
    year = datetime.now().year
    nyx_tpl.layout.annotations = [
        dict(
            name="watermark",
            text=f"{who} © {year}",
            textangle=0,
            opacity=0.1,
            font=dict(color="black", size=100),
            xref="paper",
            yref="paper",
            x=0.1,
            y=0.1,
            showarrow=False,
        )
    ]


year = datetime.now().year
nyx_tpl = go.layout.Template()
nyx_tpl.layout.annotations = [
    dict(
        name="watermark",
        text=f"Powered by Nyx Space © {year}",
        opacity=0.75,
        font=dict(color="#3d84e8", size=12),
        xref="paper",
        yref="paper",
        x=0.95,
        y=0.08,
        showarrow=False,
    )
]


def plot_with_error(
    fig, df, x, y_name, err_name, color_name, xtitle, ytitle, title, copyright=None
):
    """
    Adds the traces for the provided X vs Y with ERR around Y.
    """
    color = colors[color_name]

    fig.add_trace(
        go.Scatter(
            name=y_name,
            x=x,
            y=df[y_name],
            mode="lines",
            line=dict(color=f"rgb({color[0]}, {color[1]}, {color[2]})"),
        )
    )

    fig.add_trace(
        go.Scatter(
            name=f"+{y_name}",
            x=x,
            y=df[y_name] + df[err_name],
            mode="lines",
            marker=dict(
                color=f"rgb({int(color[0]*0.85)}, {int(color[1]*0.85)}, {int(color[2]*0.85)})"
            ),
            line=dict(width=0),
            showlegend=False,
        )
    )
    fig.add_trace(
        go.Scatter(
            name=f"-{y_name}",
            x=x,
            y=df[y_name] - df[err_name],
            marker=dict(
                color=f"rgb({int(color[0]*0.85)}, {int(color[1]*0.85)}, {int(color[2]*0.85)})"
            ),
            line=dict(width=0),
            mode="lines",
            fillcolor=f"rgba({color[0]}, {color[1]}, {color[1]}, 0.3)",
            fill="tonexty",
            showlegend=False,
        )
    )

    finalize_plot(fig, title, xtitle, ytitle, copyright)


def plot_line(fig, df, x, y_name, color_name, xtitle, ytitle, title, copyright):
    """
    Adds the trace for the provided X vs Y.
    """

    plot_raw_line(
        fig, x, df[y_name], y_name, color_name, xtitle, ytitle, title, copyright
    )


def plot_raw_line(fig, x, y, y_name, color_name, xtitle, ytitle, title, copyright=None):
    """
    Adds the trace for the provided X vs Y.
    """
    color = colors[color_name]

    fig.add_trace(
        go.Scatter(
            name=y_name,
            x=x,
            y=y,
            mode="lines",
            line=dict(color=f"rgb({color[0]}, {color[1]}, {color[2]})"),
        )
    )

    finalize_plot(fig, title, xtitle, ytitle, copyright)


def finalize_plot(fig, title, xtitle=None, ytitle=None, copyright=None):
    """
    Add titles, copyright and watermark
    """

    annotations = [dict(templateitemname="watermark")]
    if copyright is not None:
        annotations += [
            dict(
                templateitemname="watermark",
                text=f"{copyright} © {year}",
                font=dict(color="black", size=16),
                y=0.1,
            )
        ]

    if xtitle:
        fig.update_layout(xaxis_title=xtitle)
    if ytitle:
        fig.update_layout(yaxis_title=ytitle)

    fig.update_layout(
        title=title,
        hovermode="x",
        template=nyx_tpl,
        annotations=annotations,
    )


def build_sphere(size, color, num_points=500, opacity=1.0):
    """
    Source: https://python.plainenglish.io/how-to-create-a-3d-model-of-the-solar-system-with-plotly-in-python-2a74e672b771
    """

    theta = np.linspace(0, 2 * np.pi, num_points)
    phi = np.linspace(0, np.pi, num_points)

    # Set up coordinates for points on the sphere
    x0 = size * np.outer(np.cos(theta), np.sin(phi))
    y0 = size * np.outer(np.sin(theta), np.sin(phi))
    z0 = size * np.outer(np.ones(num_points), np.cos(phi))

    # Set up trace
    trace = go.Surface(
        x=x0,
        y=y0,
        z=z0,
        colorscale=[
            [0, f"rgba({int(color[0])}, {int(color[1])}, {int(color[2])}, {opacity})"],
            [1, f"rgba({int(color[0])}, {int(color[1])}, {int(color[2])}, {opacity})"],
        ],
        hoverinfo="none",
    )
    trace.update(showscale=False)

    return trace


def build_3dline(x, y, z, text, color=colors["sky_blue"], width=2, name=None):
    """
    Super simple, but allows consistency
    """

    trace = go.Scatter3d(
        x=x,
        y=y,
        z=z,
        marker=dict(size=0.1),
        line=dict(
            color=f"rgb({int(color[0])}, {int(color[1])}, {int(color[2])})",
            width=width,
        ),
        text=text,
        name=name,
    )
    return trace
