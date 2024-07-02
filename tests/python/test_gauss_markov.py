from nyx_space.orbit_determination import GaussMarkov
from nyx_space.time import Unit
from nyx_space.plots import plot_gauss_markov
from pathlib import Path
import pandas as pd


def test_fogm(plot=False):
    """
    Tests a first order Gauss Markov process
    """

    # Base path
    root = Path(__file__).joinpath("../../../").resolve()
    outpath = root.joinpath("output_data/")

    gm = GaussMarkov(tau=Unit.Hour * 24, sigma=5.6, steady_state=0.5)
    print(gm)
    assert str(gm) == "First order Gauss-Markov process with τ = 1 days, σ_b = 5.6, σ_q = 0.5"
    gm.simulate(str(outpath.joinpath("fogm.parquet")))

    # Read in the file
    df = pd.read_parquet(str(outpath.joinpath("fogm.parquet")))

    # The seeds are chosen from entropy each time,
    # so the bias should be non-zero but we can't predict what it will be.
    assert df["Bias (unitless)"].mean() != 0.0
    assert df["Bias (unitless)"].max() != 0.0
    assert df["Bias (unitless)"].min() != 0.0
    assert df["Bias (unitless)"].count() == 12525

    if plot:
        plot_gauss_markov(df, title=f"{gm}", tau=gm.tau.to_seconds())


def test_defaults(kinds=["Range", "RangeHP", "Doppler", "DopplerHP"], plot=False):
    """
    Tests the four default models
    """

    # Base path
    root = Path(__file__).joinpath("../../../").resolve()
    outpath = root.joinpath("output_data/")

    for kind in kinds:
        unit = "km" if "Range" in kind else "km/s"
        gm = GaussMarkov.default(kind)
        gm.simulate(str(outpath.joinpath(f"{kind}.parquet")), unit=unit)

        # Read in the file
        df = pd.read_parquet(str(outpath.joinpath(f"{kind}.parquet")))

        if plot:
            plot_gauss_markov(df, title=f"{kind} - {gm}", tau=gm.tau.to_seconds())


def test_load(plot=False):
    """
    Tests that we can load a Gauss Markov process from a file and simulate it.
    """
    # Base path
    root = Path(__file__).joinpath("../../../").resolve()
    config = root.joinpath("data/tests/config/high-prec-network.yaml")
    outpath = root.joinpath("output_data/")

    models = GaussMarkov.load_named(str(config))

    for name, gm in models.items():
        gm.simulate(str(outpath.joinpath(f"{name}.parquet")))

        # Read in the file
        df = pd.read_parquet(str(outpath.joinpath(f"{name}.parquet")))

        if plot:
            plot_gauss_markov(df, title=f"{name} - {gm}", tau=gm.tau.to_seconds())


if __name__ == "__main__":
    test_fogm(True)
    test_defaults(plot=True)
    test_load(True)
