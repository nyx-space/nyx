import argparse
import pandas as pd
import plotly.graph_objects as go

from plt_utils import (
    radii,
    build_sphere,
    build_3dline,
    colors,
    finalize_plot,
    body_color,
)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot a trajectory in 3D")
    parser.add_argument(
        "pq",
        metavar="N",
        type=str,
        nargs="+",
        help="list of parquet files, each will have a different color (hopefully)",
    )
    parser.add_argument(
        "-c",
        "--center",
        type=str,
        default="Earth",
        required=False,
        help="The name of the center object, e.g. `Luna` (Moon)",
    )
    parser.add_argument(
        "--opacity",
        type=float,
        default=0.8,
        help="Opacity of the center body sphere",
    )
    parser.add_argument("--html", type=str, help="Save the HTML to a file")
    args = parser.parse_args()

    # Load all of the CSV files in data frames
    fname = ""
    try:
        traces = [
            build_sphere(
                radii[args.center], body_color[args.center], opacity=args.opacity
            )
        ]
    except KeyError:
        print("Body must be one of:" + ", ".join(radii.keys()))
        raise

    color_values = list(colors.values())
    for i, fpath in enumerate(args.pq):
        df = pd.read_parquet(fpath)
        print(df.describe())
        print(df.columns)
        traces += [
            build_3dline(
                df["x (km)"],
                df["y (km)"],
                df["z (km)"],
                df["Epoch:Gregorian UTC"],
                color=color_values[i % len(color_values)],
                name=fpath.split("/")[-1],
            )
        ]

    # Now that we have the data, let's plot it
    fig = go.Figure(data=traces)
    finalize_plot(fig, title="3D trajectory plot")
    fig.update_layout(scene_aspectmode="data")
    fig.show()

    if args.html:
        with open(args.html, "w") as f:
            f.write(fig.to_html())
        print(f"Saved HTML to {args.html}")
