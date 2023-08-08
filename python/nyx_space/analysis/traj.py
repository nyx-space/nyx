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

import pandas as pd
from nyx_space.time import Unit, TimeSeries
from nyx_space.mission_design import TrajectoryLoader
from nyx_space.monte_carlo import StateParameter


def diff_traj_parquet(path1: str, path2: str, step=Unit.Minute * 1) -> pd.DataFrame:
    """
    Load both trajectory files from parquet, sample them at the provided step, and return an error data frame that can be plotted with the typical tools.
    """

    oe_list_str = [
        "sma",
        "ecc",
        "inc",
        "raan",
        "aop",
        "ta",
        "tlong",
        "aol",
        "x",
        "y",
        "z",
        "vx",
        "vy",
        "vz",
    ]

    oe_list = [StateParameter(name) for name in oe_list_str]

    def loader(path: str) -> TrajectoryLoader:
        # Try to be somewhat clever
        if path.lower().endswith(".oem"):
            # Load as OEM and build the parquet file
            return TrajectoryLoader(
                str(path), "oem", path[:-3] + "parquet"
            ).to_orbit_traj()
        else:
            return TrajectoryLoader(str(path), "parquet").to_orbit_traj()

    traj1 = loader(path1)
    traj2 = loader(path2)

    err_data = {"Epoch:Gregorian UTC": []}
    for oe in oe_list:
        err_data[str(oe)] = []

    # Sample to build the data.

    for epoch in TimeSeries(
        traj1.first().epoch, traj1.last().epoch, step, inclusive=True
    ):
        try:
            state1 = traj1.at(epoch)
            state2 = traj2.at(epoch)
        except:
            # Missing data at either of these, let's ignore
            pass
        else:
            err_data["Epoch:Gregorian UTC"] += [str(epoch)]
            for oe in oe_list:
                val1 = state1.value_of(oe)
                val2 = state2.value_of(oe)
                err_data[str(oe)] += [val1 - val2]

    # Build the data frame
    err_df = pd.DataFrame(err_data, columns=list(err_data.keys()))

    return err_df
