"""
Nyx, blazing fast astrodynamics
Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

import plotly.express as px

from .utils import finalize_plot


def plot_gauss_markov(df, title="Gauss Markov Process", tau=None):
    """
    Plots a Gauss Markov process from the provided data frame
    """

    # Grab the column
    try:
        col_name = [col for col in df.columns if "Bias" in col][0]
        variance_name = [col for col in df.columns if "Variance" in col][0]
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

    fig = px.line(
        df,
        x="Delta Time (s)",
        y=variance_name,
        color="Run",
        opacity=0.5,
    )

    finalize_plot(fig, title=title)

    fig.show()
