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

from .utils import (
    radii,
    build_sphere,
    build_3dline,
    colors,
    finalize_plot,
    body_color,
)


def plot_traj(
    dfs,
    title,
    html_out=None,
    copyright=None,
    fig=None,
    center="Earth",
    show=True,
):
    """
    Plot a trajectory in 3D

    Args:
        dfs (pandas.DataFrame): The data frame containing the trajectory (or a list thereof)
        title (str): The title of the plot
        html_out (str, optional): The path to save the HTML to. Defaults to None.
        copyright (str, optional): The copyright to display on the plot. Defaults to None.
        fig (plotly.graph_objects.Figure, optional): The figure to add the trajectory to. Defaults to None.
        center (str, optional): The name of the center object, e.g. `Luna` (Moon). Defaults to "Earth".
        show (bool, optional): Whether to show the plot. Defaults to True. If set to false, the figure will be returned.
    """

    if fig is None:
        fig = go.Figure()

    # Only plot the sphere if there currently is nothing on the plot
    if len(fig.data) == 0:
        try:
            traces = [build_sphere(radii[center], body_color[center], opacity=0.8)]
        except KeyError:
            print("Body must be one of:" + ", ".join(radii.keys()))
            raise
    else:
        traces = []

    if not isinstance(dfs, list):
        dfs = [dfs]

    for df in dfs:
        print(df.describe())
        print(df.columns)
        color_values = list(colors.values())
        traces += [
            build_3dline(
                df["x (km)"],
                df["y (km)"],
                df["z (km)"],
                df["Epoch:Gregorian UTC"],
                color=color_values[len(fig.data) % len(color_values)],
                name=title,
            )
        ]

    # Now that we have the data, let's plot it
    fig.add_traces(traces)
    finalize_plot(fig, title, copyright=copyright)
    fig.update_layout(scene_aspectmode="data")

    if html_out:
        with open(html_out, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {html_out}")

    if show:
        fig.show()
    else:
        return fig


def plot_ground_track(
    dfs,
    title,
    names,
    projection="orthographic",
    landmarks={},
    not_earth=False,
    html_out=None,
    copyright=None,
    fig=None,
    show=True,
):
    """
    Plot a ground track

    Args:
        dfs (pandas.DataFrame): The data frame containing the trajectory (or a list thereof)
        title (str): The title of the plot
        names (str): The name of the trajectory (or a list thereof)
        projection (str, optional): The projection to use. Defaults to "orthographic".
        landmarks (dict, optional): A dictionary of landmarks to plot. Defaults to {}. Example: {"Eiffel Tower": (48.8584, 2.2945)} (latitude, longitude)
        not_earth (bool, optional): Whether to plot the Earth. Defaults to False.
        html_out (str, optional): The path to save the HTML to. Defaults to None.
        copyright (str, optional): The copyright to display on the plot. Defaults to None.
        fig (plotly.graph_objects.Figure, optional): The figure to add the trajectory to. Defaults to None.
        show (bool, optional): Whether to show the plot. Defaults to True. If set to false, the figure will be returned.
    """
    if fig is None:
        fig = go.Figure()

    color_values = list(colors.values())

    if not isinstance(dfs, list):
        dfs = [dfs]

    if not isinstance(names, list):
        names = [names]

    assert len(names) == len(dfs), "Number of names must match number of dataframes"

    for (name, df) in zip(names, dfs):
        i = len(fig.data)
        try:
            color = color_values[i % len(color_values)]
            fig.add_trace(
                go.Scattergeo(
                    lon=df["geodetic_longitude (deg)"],
                    lat=df["geodetic_latitude (deg)"],
                    mode="lines",
                    name=name,
                    text=df["Epoch:Gregorian UTC"],
                    line=dict(
                        color=f"rgb({int(color[0])}, {int(color[1])}, {int(color[2])})",
                    ),
                )
            )
        except KeyError as e:
            print(
                "\n\nLatitude and longitude not found in the data.\n"
                "Generate the output with traj.to_csv_with_groundtrack(...), \n"
                'Or specify `"geodetic_latitude", "geodetic_longitude"` in the output headers\n\n'
            )
            raise e

    for name, (lat, long) in landmarks.items():
        print(f"Adding landmark {name} at lat={lat}, long={long}")
        fig.add_trace(
            go.Scattergeo(lon=[float(long)], lat=[float(lat)], name=name.strip())
        )

    fig.update_geos(
        lataxis_showgrid=True,
        lonaxis_showgrid=True,
        projection_type=projection,
        visible=not not_earth,
    )
    fig.update_layout(margin={"r": 0, "t": 0, "l": 0, "b": 0})

    finalize_plot(fig, f"Ground track: {title}", copyright=copyright)

    if html_out:
        with open(html_out, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {html_out}")

    if show:
        fig.show()
    else:
        return fig
