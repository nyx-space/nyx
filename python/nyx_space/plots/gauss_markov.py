import plotly.express as px

from .utils import finalize_plot


def plot_gauss_markov(df, title="Gauss Markov Process", tau=None):
    """
    Plots a Gauss Markov process from the provided data frame
    """

    # Grab the column
    try:
        col_name = [col for col in df.columns if "Bias" in col][0]
    except IndexError:
        raise ValueError("No bias column found in the provided data frame")

    fig = px.scatter(
        df,
        x="Delta Time (s)",
        y=col_name,
        color="Run",
        opacity=0.4,
        marginal_y="rug",
    )

    if tau:
        fig.add_vline(
            x=tau,
            line_width=2,
            line_dash="dash",
            line_color="red",
            row=1,
            col=1,
        )

    finalize_plot(fig, title=title)

    fig.show()
