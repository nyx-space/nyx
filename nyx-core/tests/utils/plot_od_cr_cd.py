import click
import plotly.express as px
import plotly.graph_objects as go
import polars as pl
from plotly.subplots import make_subplots

TEMPLATE = "seaborn"


@click.command
@click.option(
    "-p",
    "--path",
    type=str,
    required=True,
    help="Path to the parquet file containing OD results",
)
def main(path: str):
    df = pl.read_parquet(path)

    df = df.with_columns(
        pl.col("Epoch (UTC)").str.to_datetime("%Y-%m-%dT%H:%M:%S%.f")
    ).sort("Epoch (UTC)", descending=False)

    cols = df.columns
    cr_col = "cr" if "cr" in cols else None
    cd_col = "cd" if "cd" in cols else None

    sigma_cr_col = next(
        (col for col in cols if col.startswith("Sigma Cr") and "(unitless)" in col),
        None,
    )
    sigma_cd_col = next(
        (col for col in cols if col.startswith("Sigma Cd") and "(unitless)" in col),
        None,
    )

    plots_to_make = []
    if cr_col:
        plots_to_make.append(("Cr Estimate", cr_col, sigma_cr_col))
    if cd_col:
        plots_to_make.append(("Cd Estimate", cd_col, sigma_cd_col))

    if not plots_to_make:
        print("No Cr or Cd columns found in the data.")
        return

    fig = make_subplots(
        rows=len(plots_to_make),
        cols=1,
        subplot_titles=[p[0] for p in plots_to_make],
        vertical_spacing=0.1,
    )

    legend_added = False

    for idx, (title, val_col, sigma_col) in enumerate(plots_to_make, start=1):
        # Add the estimate trace
        fig.add_trace(
            go.Scatter(
                x=df["Epoch (UTC)"],
                y=df[val_col],
                mode="lines+markers",
                name=title,
                legendgroup=title,
                marker=dict(color="blue" if "Cr" in title else "green"),
                showlegend=True,
            ),
            row=idx,
            col=1,
        )

        # Add 3-sigma bounds if available
        if sigma_col:
            df = df.with_columns(
                [
                    (pl.col(val_col) + 3.0 * pl.col(sigma_col)).alias(
                        f"{title} +3-Sigma"
                    ),
                    (pl.col(val_col) - 3.0 * pl.col(sigma_col)).alias(
                        f"{title} -3-Sigma"
                    ),
                ]
            )
            for bound in [f"{title} +3-Sigma", f"{title} -3-Sigma"]:
                show_this_legend = not legend_added
                fig.add_trace(
                    go.Scatter(
                        x=df["Epoch (UTC)"],
                        y=df[bound],
                        mode="lines",
                        name="3-Sigma bounds",
                        line=dict(color="black", dash="dash"),
                        legendgroup="3-Sigma bounds",
                        connectgaps=True,
                        showlegend=show_this_legend,
                    ),
                    row=idx,
                    col=1,
                )
                legend_added = True

        fig.update_yaxes(title_text="Value (unitless)", row=idx, col=1)

    fig.update_layout(title_text="Cr and Cd Estimates", template=TEMPLATE)
    fig.update_xaxes(matches="x")
    fig.show()


if __name__ == "__main__":
    main()
