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

from dataclasses import dataclass
from pathlib import Path

import yaml
from nyx_space.time import Epoch, Unit

_default_range_noise_km = {
    "tau": "24 h",
    "bias_sigma": 5.0e-3,  # 5 m
    "steady_state_sigma": 0.1e-3,  # 0.1 m
}

_default_doppler_noise_km_s = {
    "tau": "24 h",
    "bias_sigma": 50.0e-6,  # 5 cm/s
    "steady_state_sigma": 1.5e-6,  # 0.15 cm/s
}

f_2ghz = 2.036250e9  # Hz -- the  2 GHz frequency of DORIS
f_400mhz = 250e6  # Hz -- the 400 MHz frequency of DORIS


@dataclass
class RnxMsr:
    """
    Stores a RINEX measurement
    """

    epoch: Epoch  # Includes the clock offset information
    station_name: str
    l1_phase_cycles: float  #  2GHz phase cycle data
    l2_phase_cycles: float  # 400MHz phase cycle data
    c1_pseudo_range_km: float  # 2GHz pseudo range data
    c2_pseudo_range_km: float  # 400MHz pseudo range data
    rel_freq_offset: float  # in 1e-11


def _to_dec(in_str):
    deg, min, sec = [float(x) for x in in_str.split()]
    return deg + min / 60.0 + sec / 3600.0


def snx_to_groundstation(
    snx_path,
    output_path,
    frame="IAU Earth",
    elevation_mask_deg=0.0,
    range_noise_km=_default_range_noise_km,
    doppler_noise_km_s=_default_doppler_noise_km_s,
):
    """
    Convert a SINEX file to a ground station YAML file.

    This only converts the data between the `+SITE/ID` and `-SITE/ID` lines.
    """

    offsets = {
        "name": (0, 5),
        "PT": (6, 2),
        "__DOMES__": (9, 9),
        "T": (10, 1),
        "description": (21, 22),
        "longitude_deg": (44, 11),
        "latitude_deg": (56, 11),
        "height_km": (68, 7),
    }

    entries = []

    with open(snx_path, "r") as snx_file:
        parsing = False
        for line in snx_file.readlines():
            if not parsing and line.startswith("+SITE/ID"):
                parsing = True
                continue
            elif parsing:
                if line.startswith("-SITE/ID"):
                    print(
                        "Parsed {} ground stations from {}.".format(
                            len(entries), snx_path
                        )
                    )
                    break
                elif line.startswith("*CODE"):
                    continue
                else:
                    # Parse the data line
                    kwrds = [
                        line[idx : idx + klen].strip() for idx, klen in offsets.values()
                    ]
                    # Create the YAML
                    this_entry = dict(zip(offsets.keys(), kwrds))
                    # Fix the latitude and longitude to be decimals
                    this_entry["longitude_deg"] = _to_dec(this_entry["longitude_deg"])
                    this_entry["latitude_deg"] = _to_dec(this_entry["latitude_deg"])
                    # Convert the height to km
                    this_entry["height_km"] = float(this_entry["height_km"]) / 1000.0

                    # Append the Nyx specific data
                    this_entry["frame"] = frame
                    this_entry["elevation_mask_deg"] = elevation_mask_deg
                    this_entry["range_noise_km"] = range_noise_km
                    this_entry["doppler_noise_km_s"] = doppler_noise_km_s
                    this_entry["light_time_correction"] = True
                    # Remove the unused keys
                    del this_entry["__DOMES__"]
                    del this_entry["T"]
                    del this_entry["PT"]

                    entries += [this_entry]

    output_path = str(Path(output_path).resolve())

    # Dump the YAML without creating anchors
    # cf. https://github.com/yaml/pyyaml/issues/103
    class NoAliasDumper(yaml.Dumper):
        def ignore_aliases(self, _data):
            return True

    yaml.dump(entries, open(output_path, "w"), Dumper=NoAliasDumper)
    print("Stored {} ground stations in {}.".format(len(entries), output_path))


def basic_rinex(rnx_path, stations):
    """
    Parse a RINEX file to extract the basic information.
    """

    action = "hdr"

    rnx_stations = {}
    data = []

    with open(rnx_path) as rnx_file:
        for lno, line in enumerate(rnx_file.readlines()):
            if action == "hdr":
                if line[0] == "D" and len(line.split()[0]) == 3:
                    # We've found the first station
                    action = "stations"
            if action == "stations":
                # Parse the station data
                if line.startswith(" "):
                    # End of the stations
                    action = "new_record"
                    continue
                splt = line.split()
                # The following will raise an error if the station code is not in the stations dictionary
                # Note that we ignore the dome number and the frequency shift factor
                rnx_stations[splt[0]] = stations[splt[1]]

            elif action == "msr_start":
                if line[0] == ">":
                    action = "new_record"
                else:
                    # This is a new entry in the measurement block
                    splt = line.split()
                    # Note that we set the relative frequency offset to 0.0 but we'll modify it at the next line
                    try:
                        msr = RnxMsr(
                            station_name=rnx_stations[splt[0]].name,
                            epoch=epoch,
                            l1_phase_cycles=float(splt[1]),
                            l2_phase_cycles=float(splt[2]),
                            # Documentation says the data is in units of 0.01 km
                            c1_pseudo_range_km=float(splt[3]) / 100.0,
                            c2_pseudo_range_km=float(splt[4]) / 100.0,
                            rel_freq_offset=0.0,
                        )
                    except ValueError:
                        print(f"Error at line {lno} - stopping here")
                        break
                    data += [msr]
                    action = "msr_continued"
            elif action == "msr_continued":
                splt = line.split()
                if len(splt) != 9:
                    # Something is wrong, let's remove this measurement
                    data.pop()
                    action = "msr_start"
                    continue
                rel_freq_offset = float(splt[2])
                msr = data.pop()
                msr.rel_freq_offset = rel_freq_offset
                data += [msr]
                action = "msr_start"

            if action == "new_record" and line[0] == ">":
                # Compute the epoch of this record
                splt = line.split()
                if int(splt[-1]) != 0:
                    print(f"Skipping unusable measurement block: {line}")
                    continue
                offset = Unit.Second * float(splt[-2])
                esplt = splt[1:-4]
                epoch = (
                    Epoch(
                        f"{esplt[0]}-{esplt[1]}-{esplt[2]} {esplt[3]}:{esplt[4]}:{esplt[5]} TAI"
                    )
                    + offset
                )
                action = "msr_start"

    print("Parsed {} measurements from {}.".format(len(data), rnx_path))
    return data


def prototype(gs_yaml, rnx_path, plot=True):
    import numpy as np
    import plotly.graph_objects as go
    from nyx_space.orbit_determination import GroundStation

    stations = GroundStation.load_many(gs_yaml)
    ground_stations = {}
    for station in stations:
        ground_stations[station.name] = station
    data = basic_rinex(rnx_path, ground_stations)
    gs_to_msrs = {}
    for msr in data:
        if msr.station_name not in gs_to_msrs:
            gs_to_msrs[msr.station_name] = []
        gs_to_msrs[msr.station_name] += [msr]

    # Build the time arrays one ground station at a time
    for gs_name, msrs in gs_to_msrs.items():
        time_offsets = [
            0.0,
        ]
        c1_data = []
        ref_epoch = msrs[0].epoch
        for msr in msrs[1:]:
            if msr.c1_pseudo_range_km > 0.0:
                # This data is invalid
                continue
            time_offsets += [msr.epoch.timedelta(ref_epoch).to_seconds()]
            c1_data += [msr.c1_pseudo_range_km]

        # Some preprocessing to try to make this zero meaned
        c1_data = (np.array(c1_data) - np.mean(c1_data)) / np.std(c1_data)

        try:
            c1_corr = np.correlate(time_offsets, c1_data, mode="same")
        except ValueError:
            print(f"Skipping {gs_name} due to a ValueError")
            continue

        tau_val = np.max(c1_corr) / np.exp(1)
        tau = time_offsets[np.argmin(np.abs(c1_corr - tau_val))]
        # Compute the stddev of the data after 2*tau
        c1_stddev = np.std([x for t, x in zip(time_offsets, c1_data) if t > 2 * tau])

        dt = msrs[-1].epoch.timedelta(ref_epoch)
        title = f"For {ground_stations[gs_name]}: τ = {tau} s, σ = {c1_stddev * np.std(c1_data)} km (with {len(msrs)} measurements over {dt})"
        print(title)
        print(f"\tC1: {c1_corr}")
        if len(msrs) > 100 and plot:
            fig = go.Figure()
            fig.add_trace(
                go.Scatter(
                    x=time_offsets, y=c1_corr, mode="lines", name="C1 data (2 GHz)"
                )
            )
            # Add τ
            fig.add_vline(
                x=tau,
                line_width=2,
                line_dash="dash",
                line_color="red",
                row=1,
                col=1,
            )
            fig.update_layout(title=title)
            fig.show()
            # Plot the data itself too
            fig = go.Figure()
            fig.add_trace(
                go.Scatter(
                    x=time_offsets, y=c1_data, mode="lines", name="C1 data (2 GHz)"
                )
            )
            fig.update_layout(title=title)
            fig.show()

    """
    Iterate through different true anomaly values and check which lead to the least errors
    """
    import logging

    import numpy as np
    import pandas as pd
    from nyx_space.cosmic import Cosm, Orbit, Spacecraft
    from nyx_space.mission_design import (
        TrajectoryLoader,
        SpacecraftDynamics,
        propagate,
    )
    from nyx_space.orbit_determination import GroundTrackingArcSim, TrkConfig

    FORMAT = "%(levelname)s %(name)s %(asctime)-15s %(filename)s:%(lineno)d %(message)s"
    logging.basicConfig(format=FORMAT)
    logging.getLogger().setLevel(logging.INFO)

    cosm = Cosm.de438()
    eme2k = cosm.frame("EME2000")

    # Base path
    root = Path(__file__).joinpath("../../../../").resolve()
    config_path = root.joinpath("./data/tests/config/")
    outpath = root.joinpath("output_data/")

    # Load the dynamics
    dynamics = SpacecraftDynamics.load_named(str(config_path.joinpath("dynamics.yaml")))

    for gs_name, msrs in gs_to_msrs.items():
        if len(msrs) < 100:
            continue
        # Build the approximate orbit
        # We'll use the first measurement as the reference epoch
        ref_epoch = msrs[0].epoch

        # We'll use the first measurement as the reference station
        ref_station = stations[gs_name]
        # Setup a continuous tracking for this station
        trk_cfg = {ref_station.name: TrkConfig()}

        for ta_deg in np.linspace(0, 360, 10):
            orbit = Orbit.from_keplerian(
                sma_km=7349136.3e-3,
                ecc=0.00117,
                inc_deg=99.35,
                raan_deg=108.48,
                aop_deg=110.9766,
                ta_deg=ta_deg,
                epoch=ref_epoch,
                frame=eme2k,
            )

            # Build a spacecraft
            sc = Spacecraft(
                orbit,
                dry_mass_kg=1500.0,
                fuel_mass_kg=0.0,
            )

            # Propagate this trajectory for one day
            # An propagate for two periods (we only care about the trajectory)
            _, traj = propagate(sc, dynamics["hifi"], Unit.Day * 1)
            print(traj)

            # This is a bit of a pain but at the moment, we must serialize the trajectory to a file
            # and reload it into a Dynamic trajectory before passing it to the Python object.

            traj_file = str(outpath.joinpath("./od_val_with_arc_truth_ephem.parquet"))
            traj.to_parquet(traj_file)
            traj = TrajectoryLoader(traj_file)

            # Build a tracking arc
            arc_sim = GroundTrackingArcSim([ref_station], traj, trk_cfg, 0)
            # Generate the measurements
            path = arc_sim.generate_measurements(
                str(outpath.joinpath(f"./from_{ref_station.name}-ta_{ta_deg}.parquet")),
                timestamp=False,
                metadata={},
            )
            # Load this parquet and print the head
            df = pd.read_parquet(path)
            print(df.head())

            break
        break
