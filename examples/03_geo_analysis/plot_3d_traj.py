import argparse
import polars as pl
import plotly.graph_objs as go
import numpy as np


def build_sphere(size, num_points=500, opacity=1.0):
    """
    Source: https://python.plainenglish.io/how-to-create-a-3d-model-of-the-solar-system-with-plotly-in-python-2a74e672b771
    """

    color = (86, 180, 233)  # #56b4e9

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


def plot_traj(
    df: pl.DataFrame,
    colored_by="fuel_mass (kg)",
    color_descr="fuel mass (kg)",
    scale=1.0
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

    df = df.with_columns(pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")).sort(
        "Epoch (UTC)", descending=False
    )

    fig = go.Figure()

    earth_radius = 6378.136300

    traces = [build_sphere(earth_radius, opacity=0.8)]

    traces += [
        go.Scatter3d(
            x=df["x (km)"],
            y=df["y (km)"],
            z=df["z (km)"],
            mode="lines",
            line=dict(
                color=df[colored_by]*scale,
                colorscale="Viridis",
                colorbar=dict(title=color_descr),
            ),
            name="GEO",
        )
    ]

    # Now that we have the data, let's plot it
    fig.add_traces(traces)
    fig.update_layout(scene_aspectmode="data")

    fig.show()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="3D traj plotter")
    parser.add_argument("pq", type=str, help="Path to the parquet file")
    args = parser.parse_args()
    
    df = pl.read_parquet(args.pq)
    plot_traj(df)
    plot_traj(df, colored_by="penumbra event light-source: Sun J2000, shadows casted by: Earth J2000, Moon J2000", color_descr="Illumination %", scale=100.0)
