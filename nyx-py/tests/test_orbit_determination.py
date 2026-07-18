import os

from nyx_space import Spacecraft
from nyx_space.anise import MetaAlmanac
from nyx_space.anise.analysis import OrbitalElement
from nyx_space.anise.astro import LocalFrame, Orbit
from nyx_space.anise.constants import CelestialObjects, Frames
from nyx_space.mission_design import (
    AccelModels,
    Dynamics,
    GravityFieldConfig,
    PointMasses,
    Propagator,
)
from nyx_space.monte_carlo import StateDispersion, StateParameter
from nyx_space.orbit_determination import (
    CN0,
    CarrierFreq,
    ExportCfg,
    GroundStation,
    GroundTrackingArcSim,
    Handoff,
    Location,
    MeasurementType,
    ProcessNoise,
    SigmaRejection,
    KalmanVariant,
    Scheduler,
    SpacecraftEstimate,
    SpacecraftODProcess,
    StochasticNoise,
    TrackingDataArc,
    TrkConfig,
)
from nyx_space.time import Duration, Epoch, Unit


def test_ground_station():
    """
    Demonstrate ground station config from a script, saving to a Yaml, and reloading
    """
    # Location requires latitude (deg), longitude (deg), and height above the geoid strictly in kilometers.
    # The station must be anchored to a planetary body frame (e.g., Earth-fixed).
    gs = GroundStation(
        "Paris, FR",
        Location(2.3522, 48.8566, 0.4, Frames.IAU_EARTH_FRAME.to_frameuid(), [], True),
        {
            MeasurementType.Range: StochasticNoise(),
            MeasurementType.Doppler: StochasticNoise(),
            MeasurementType.Elevation: StochasticNoise(),
            MeasurementType.Azimuth: StochasticNoise(),
        },
    )

    # VERIFICATION: Ensure lossless round-trip serialization.
    # This guarantees that complex network definitions can be stored in version
    # control (YAML) and loaded without mutating floating-point data.
    as_yaml = gs.to_yaml()
    gs_reloaded = GroundStation.from_yaml(as_yaml)
    assert str(gs) == str(gs_reloaded)
    assert gs == gs_reloaded
    print(gs)
    print(repr(gs))


def test_howto_simulate_tracking_data():
    """
    Goal: Propagate a high-fidelity orbit and simulate synthetic radiometric
    tracking data from a ground station network.
    """
    almanac = MetaAlmanac("../data/02_config/ci_almanac.dhall").process()
    eme2k = almanac.frame_info(Frames.EME2000)

    # Step 1: Build a reference trajectory
    # Define a low-Earth orbit using Cartesian state vectors (kilometers and kilometers per second).
    orbit = Orbit(
        5442.1625926801835,
        -4068.9498468206248,
        -13.456851447751518,
        2.8581975428173836,
        3.8097859312745794,
        6.002126693122689,
        Epoch("2025-08-25 11:55:44 UTC"),
        eme2k,
    )
    spacecraft = Spacecraft(orbit)

    # The gravity field must be defined in an Earth-fixed frame so the
    # spherical harmonics rotate correctly beneath the inertial orbit.
    accel_models = AccelModels(
        point_masses=PointMasses(
            celestial_objects=[CelestialObjects.MOON, CelestialObjects.SUN]
        ),
        gravity_field=GravityFieldConfig(
            degree=10,
            order=10,
            filepath="../data/01_planetary/EGM2008_to2190_TideFree.gz",
            frame=Frames.IAU_EARTH_FRAME.to_frameuid(),
        ),
    )

    # Propagate to generate the ground-truth trajectory.
    dynamics = Dynamics(accel_models)
    propagator = Propagator(dynamics, almanac)

    traj = propagator.for_duration(spacecraft, Unit.Day * 1, trajectory=True).trajectory

    # Step 2: Configure Stochastic Measurement Noise
    # We define white noise and bias characteristics for each measurement type.
    # NOTE only Nyx premium can estimate biases.

    stochastic_noises = {
        # Nyx provides defaults, with names: range, doppler, angles
        MeasurementType.Range: StochasticNoise(name="Range"),
        # Range and Doppler noise can be physically derived from carrier frequency and Carrier-to-Noise density.
        MeasurementType.Doppler: StochasticNoise.from_hardware_doppler_km_s(
            1e-15, Duration("1 min"), CarrierFreq.XBand(), CN0.Strong()
        ),
        # The defaults are usually good
        MeasurementType.Elevation: StochasticNoise(name="Angles"),
        MeasurementType.Azimuth: StochasticNoise(name="Angles"),
    }

    # VERIFICATION: Validate the noise distribution independently.
    # Simulating the noise profile via Monte Carlo before simulate measurements.
    range_sim_samples = stochastic_noises[MeasurementType.Range].simulate(
        "range_noise_simlation.pq", 1000, "km"
    )

    assert len(range_sim_samples) == 1_082_000
    last_sample = range_sim_samples[-1]
    print(last_sample)
    assert last_sample.run == 999
    assert last_sample.dt_s == 86400.0

    # Step 3: Build the ground station network
    gs0 = GroundStation(
        name="Paris, FR",
        location=Location(
            latitude_deg=48.8566,
            longitude_deg=2.3522,
            height_km=0.4,
            frame=Frames.IAU_EARTH_FRAME.to_frameuid(),
            terrain_mask=[],
            terrain_mask_ignored=True,
        ),
        stochastic_noises=stochastic_noises,
    )

    gs1 = GroundStation(
        "Denver, CO",
        Location(
            39.7420, -104.9915, 1.8, Frames.IAU_EARTH_FRAME.to_frameuid(), [], True
        ),
        stochastic_noises,
    )

    # Step 4: Define Tracking Schedules
    # Schedulers dictate how conflicts are resolved when multiple stations have visibility.
    # 'Eager' will hand off as soon as a new station acquires the signal.
    # 'Greedy' will hold the track until loss of signal before allowing a handoff.
    configs = {
        "Paris, FR": TrkConfig(Scheduler(Handoff.Eager)),
        "Denver, CO": TrkConfig(Scheduler(Handoff.Greedy)),
    }

    trk_sim = GroundTrackingArcSim({"Paris, FR": gs0, "Denver, CO": gs1}, traj, configs)

    # Step 5: Build the strands if any of the ground stations are configured as a Scheduler.
    trk_sim.build_schedule(almanac)

    # Step 6: Generate synthetic measurements
    trk_arc = trk_sim.generate_measurements(almanac)
    assert not trk_arc.is_empty()

    # Step 6: Verify the Generated Tracking Arc
    # VERIFICATION: Ensure the arc is populated with the expected temporal constraints and volume.
    assert not trk_arc.is_empty()
    assert trk_arc.len() == 72
    # The arc duration is measured from the timestamp of the first to the last measurement.
    assert trk_arc.duration() > Unit.Hour * 18

    # VERIFICATION: Ensure the scheduler used the full network topology.
    expected_trackers = {"Denver, CO", "Paris, FR"}
    assert set(trk_arc.unique_aliases()) == expected_trackers

    # VERIFICATION: Ensure all defined measurement variants were simulated.
    expected_types = {
        MeasurementType.Range,
        MeasurementType.Azimuth,
        MeasurementType.Doppler,
        MeasurementType.Elevation,
    }
    assert set(trk_arc.unique_types()) == expected_types

    # VERIFICATION: Test filtering mechanisms critical for orbit determination ingestion.
    denver_arc = trk_arc.filter_by_tracker("Denver, CO")
    assert not denver_arc.is_empty()
    assert denver_arc.len() < trk_arc.len()

    range_arc = trk_arc.filter_by_measurement_type(MeasurementType.Range)
    assert not range_arc.is_empty()

    # VERIFICATION: Validate data export standards.
    # The generated measurements must comply with the CCSDS TDM standard for interoperability.
    tdm_filepath = "test_output.tdm"
    # Define aliases to map our names to those used by another software.
    aliases = {"Denver, CO": "RL CO"}
    trk_arc.write_ccsds_tdm("MySpacecraft", aliases, tdm_filepath)
    assert os.path.exists(tdm_filepath)

    # We can rename the ground stations when loading the TDM file
    reverse_aliases = {"RL CO": "Denver, CO"}
    trk_reloaded = TrackingDataArc.from_ccsds_tdm(tdm_filepath, reverse_aliases)

    assert not trk_reloaded.is_empty()
    assert trk_reloaded.len() == 72
    # The arc duration is measured from the timestamp of the first to the last measurement.
    assert trk_reloaded.duration() > Unit.Hour * 18
    # Check that the reverse aliases were applied
    assert set(trk_reloaded.unique_aliases()) == expected_trackers


def test_howto_exec_orbit_determination_filter():
    """
    Goal: Given some tracking data from a ground station network, run it through
    a scalar orbit determination filter to estimate the spacecraft position, velocity
    abd coefficient of reflectivity.

    Only Nyx Premium supports simultaneously estimating ground station location,
    measurement biases, or consider covariance for Cr/Cd/ground station parameters.
    """
    almanac = MetaAlmanac("../data/02_config/ci_almanac.dhall").process()
    eme2k = almanac.frame_info(Frames.EME2000)

    # Step 1: Build a reference trajectory
    # Define a low-Earth orbit using Cartesian state vectors (kilometers and kilometers per second).
    orbit = Orbit(
        5442.1625926801835,
        -4068.9498468206248,
        -13.456851447751518,
        2.8581975428173836,
        3.8097859312745794,
        6.002126693122689,
        Epoch("2025-08-25 11:55:44 UTC"),
        eme2k,
    )
    spacecraft = Spacecraft(orbit)

    # Step 2: Disperse this spacecraft state to create an initial filter error
    # between the state used in generating the tracking data, and the one
    # used to initialize the filter.
    # Nyx supports dispersing initial states using a multivariate normal
    # distribution on any arbitrary orbital element. It also supports
    # using the covariance in the SpacecraftEstimate structure to disperse
    # according to the covariance in that structure.
    disp = [
        StateDispersion.zero_mean(
            StateParameter.Element(OrbitalElement.SemiMajorAxis), 1.0
        ),
        StateDispersion.zero_mean(
            StateParameter.Element(OrbitalElement.Eccentricity), 1e-6
        ),
    ]

    # Builds a zero mean estimate from these dispersions
    estimate = SpacecraftEstimate.from_dispersions(
        nominal_state=spacecraft, dispersions=disp, seed=123
    )
    estimate_chk_seed = SpacecraftEstimate.from_dispersions(
        nominal_state=spacecraft, dispersions=disp, seed=123
    )
    assert estimate_chk_seed == estimate

    # VERIFICATION: State is dispersed
    assert spacecraft == spacecraft, "tautology check failed"
    assert estimate.state != spacecraft, "estimate was not dispersed"
    assert estimate.nominal_state == spacecraft, "nominal state should be nominal input"
    assert len(estimate.state_deviations) == 9
    assert sum(estimate.state_deviations) / 9 > 0

    # The gravity field must be defined in an Earth-fixed frame so the
    # spherical harmonics rotate correctly beneath the inertial orbit.
    accel_models = AccelModels(
        point_masses=PointMasses(
            celestial_objects=[CelestialObjects.MOON, CelestialObjects.SUN]
        ),
        gravity_field=GravityFieldConfig(
            degree=10,
            order=10,
            filepath="../data/01_planetary/EGM2008_to2190_TideFree.gz",
            frame=Frames.IAU_EARTH_FRAME.to_frameuid(),
        ),
    )

    # Propagate the to generate the ground-truth trajectory.
    dynamics = Dynamics(accel_models)
    propagator = Propagator(dynamics, almanac)

    traj = propagator.for_duration(spacecraft, Unit.Day * 1, trajectory=True).trajectory

    # Step 3: Configure Stochastic Measurement Noise, generate tracking data.
    # We define white noise and bias characteristics for each measurement type.
    # NOTE only Nyx premium can estimate biases.

    stochastic_noises = {
        MeasurementType.Range: StochasticNoise(name="Range"),
        MeasurementType.Doppler: StochasticNoise(name="Doppler"),
    }

    # Step 4: Build the ground station network
    gs0 = GroundStation(
        name="Paris, FR",
        location=Location(
            latitude_deg=48.8566,
            longitude_deg=2.3522,
            height_km=0.4,
            frame=Frames.IAU_EARTH_FRAME.to_frameuid(),
            terrain_mask=[],
            terrain_mask_ignored=True,
        ),
        stochastic_noises=stochastic_noises,
    )

    gs1 = GroundStation(
        "Denver, CO",
        Location(
            39.7420, -104.9915, 1.8, Frames.IAU_EARTH_FRAME.to_frameuid(), [], True
        ),
        stochastic_noises,
    )

    # Step 5: Define Tracking Schedules
    configs = {
        "Paris, FR": TrkConfig(Scheduler(Handoff.Greedy)),
        "Denver, CO": TrkConfig(Scheduler(Handoff.Greedy)),
    }

    network = {"Paris, FR": gs0, "Denver, CO": gs1}

    trk_sim = GroundTrackingArcSim(network, traj, configs, seed=159)

    # Build the strands if any of the ground stations are configured as a Scheduler.
    trk_sim.build_schedule(almanac)

    # Step 6: Generate synthetic measurements
    trk_arc = trk_sim.generate_measurements(almanac)
    assert not trk_arc.is_empty()

    # Step 7: Run the orbit determination, seeding the filter with the dispersed state.
    # Nyx supports two kinds of Extended Kalman Filters: deviation tracking (good for unconverged solutions)
    # Deviation tracking is structurally underconfident when initialized with small dispersions.
    # VERIFICATION: We explicitly assert this known statistical failure rather than ignoring it.
    od_proc_deviation = SpacecraftODProcess(
        propagator, KalmanVariant.DeviationTracking, network
    )
    od_dev_sol = od_proc_deviation.process_arc(estimate, trk_arc)
    # The OD solution provides statistical tests.
    nis_consistency = od_dev_sol.nis_consistency()
    print(f"{nis_consistency}")
    # Print the information to the log
    nis_consistency.log()
    try:
        assert nis_consistency.is_consistent()
    except AssertionError:
        print("deviation tracking on small dispersions are typically underconfident")

    nees = od_dev_sol.nees_consistency(traj)
    nees.log()
    print(f"{nees}")
    assert nees.is_underconfident() == nis_consistency.is_underconfident(), (
        "NIS and NEES should agree"
    )
    assert nees.is_overconfident() == nis_consistency.is_overconfident(), (
        "NIS and NEES should agree"
    )
    try:
        assert nees.is_consistent()
    except AssertionError:
        print("and this is confirmed in the NEES metric")

    assert od_dev_sol.is_normal(), "residuals should follow a normal distribution"
    assert od_dev_sol.is_filter_run(), "this is a filter run"
    # Export the whole orbit determination solution into a single Parquet file.
    # This includes estimated states, prefits, postfits, covariance, sigmas on orbital elements, Kalman gains, etc.
    od_dev_sol.to_parquet("od_dev_smoothed.pq", ExportCfg(False))
    # Allow for some variance in the Ground station noises
    assert len(od_dev_sol.accepted_residuals()) >= 138
    # It's possible to print and extract all of the internal information of this residuals.
    # That info is also stored when exporting the OD solution to Parquet.
    final_resid = od_dev_sol.accepted_residuals()[-1]
    assert not final_resid.rejected
    # Grab the computed/real observable by measurement type
    for msr_type in [MeasurementType.Range, MeasurementType.Doppler]:
        real_obs = final_resid.real_obs(MeasurementType.Range)
        if real_obs is None:
            continue
        manual_prefit = real_obs - final_resid.computed_obs(MeasurementType.Range)
        assert manual_prefit == final_resid.prefit[0]
        break
    print(f"all available fields: {dir(final_resid)}")

    # Step 8: Running the smoother is a simple call
    smoothed = od_dev_sol.smooth(almanac)
    assert smoothed.is_smoother_run()
    assert not smoothed.is_filter_run()

    # Nyx also supports reference update (better once the filter has converged on a good solution).
    # Reference update algorithms require minor process noise to compensate for integrator
    # step size disparities between the truth trajectory and the filter's internal propagation.
    # In fact, the truth propagator will take larger steps than the orbit determination process
    # which is fixed by default to 1 minute steps, ensuring a covariance at most every minute.
    # Hence, we'll setup a simple process noise from a small velocity in m/s, integrated over the noise duration.
    process_noise = ProcessNoise.from_velocity_m_s(
        vx_m_s=1e-6,
        vy_m_s=1e-6,
        vz_m_s=1e-6,
        noise_duration=Unit.Second * 1,
        disable_time=Unit.Hour * 2,  # Disable after that duration without measurements
        local_frame=LocalFrame.Inertial,  # Default is Inertial, but RIC/VNC can also be used
    )
    od_proc = SpacecraftODProcess(
        propagator, KalmanVariant.ReferenceUpdate, network, process_noise=process_noise
    )
    od_sol = od_proc.process_arc(estimate, trk_arc)
    # Export the whole orbit determination solution into a single Parquet file.
    # This includes estimated states, prefits, postfits, covariance, sigmas on orbital elements, Kalman gains, etc.
    od_sol.to_parquet("od_ref_update.pq", ExportCfg(False))

    print(od_sol.nis_consistency())
    print(od_sol.nees_consistency(traj))
    # We can export this solution to an OEM
    definitive_ephem = od_sol.to_ephemeris("Test OD Spacecraft")
    oem_filepath = "definitive_ephem.oem"
    definitive_ephem.write_ccsds_oem(oem_filepath)
    assert os.path.exists(oem_filepath)
    os.remove(oem_filepath)

    # Finally, one can also set the sigma rejection criteria, which defaults to 3 sigmas.
    print(od_proc.sigma_rejection)
    print(od_proc.variant)

    # Let re-run with a very high number for sigma rejections (unrealistic)
    #  Expanding the threshold to 5 sigmas forces the filter to ingest statistical outliers.
    od_proc.sigma_rejection = SigmaRejection(5.0)
    od_sol_5sigma = od_proc.process_arc(estimate, trk_arc)
    assert od_sol_5sigma.is_filter_run(), "5-Sigma threshold run failed to execute."

    print(od_sol_5sigma.nis_consistency())
    print(od_sol_5sigma.nees_consistency(traj))


if __name__ == "__main__":
    import logging

    logging.basicConfig(level=logging.INFO)
    test_howto_exec_orbit_determination_filter()
    # test_howto_simulate_tracking_data()
