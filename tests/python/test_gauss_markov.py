from nyx_space.orbit_determination import GaussMarkov
from nyx_space.time import Unit
from pathlib import Path
import pandas as pd
import plotly.express as px


def test_fogm(plot=False):
    """
    Tests a first order Gauss Markov process
    """

    # Base path
    root = Path(__file__).joinpath("../../../").resolve()
    outpath = root.joinpath("output_data/")

    gm = GaussMarkov(tau=Unit.Hour * 24, sigma=0.1, steady_state=0.5)
    print(gm)
    assert (
        str(gm)
        == "First order Gauss-Markov process with τ = 1 days, σ = 0.1, q = 0.5: bias = None"
    )
    gm.simulate(str(outpath.joinpath("fogm.parquet")))

    # Read in the file
    df = pd.read_parquet(str(outpath.joinpath("fogm.parquet")))

    # The seeds are chosen from entry each time, so the bias should be non-zero but we can't predict what it will be.
    assert df["Bias (unitless)"].mean() != 0.0
    assert df["Bias (unitless)"].max() != 0.0
    assert df["Bias (unitless)"].min() != 0.0
    # assert df["Bias (unitless)"].count() == 300100

    if plot:
        fig = px.scatter(
            df,
            x="Delta Time (s)",
            y="Bias (unitless)",
            color="Run",
            opacity=0.4,
        )
        fig.update_layout(title=f"{gm}")


def test_defaults(kinds=["Range", "RangeHP", "Doppler", "DopplerHP"], plot=False):
    """
    Tests the four default models
    """

    # Base path
    root = Path(__file__).joinpath("../../../").resolve()
    outpath = root.joinpath("output_data/")

    for kind in kinds:
        gm = GaussMarkov.from_default(kind)
        gm.simulate(str(outpath.joinpath(f"{kind}.parquet")))

        # Read in the file
        df = pd.read_parquet(str(outpath.joinpath(f"{kind}.parquet")))

        if plot:
            fig = px.scatter(
                df, x="Delta Time (s)", y="Bias (unitless)", color="Run", opacity=0.25
            )
            fig.update_layout(title=f"{kind} bias: {gm}")
            fig.add_vline(
                x=gm.tau.to_seconds(), line_width=2, line_dash="dash", line_color="red"
            )
            fig.show()


if __name__ == "__main__":
    test_fogm(True)
    test_defaults(["Range"], True)
