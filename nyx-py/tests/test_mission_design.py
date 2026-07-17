import logging

import polars as pl

from nyx_space import DragData, ExportCfg, Mass, Spacecraft, SRPData
from nyx_space.anise import MetaAlmanac
from nyx_space.anise.analysis import Condition, Event, OrbitalElement, ScalarExpr
from nyx_space.anise.astro import Orbit
from nyx_space.anise.constants import CelestialObjects, Frames
from nyx_space.mission_design import (
    AccelModels,
    AtmDensity,
    Drag,
    Dynamics,
    ForceModels,
    GravityFieldConfig,
    IntegratorMethod,
    IntegratorOptions,
    PointMasses,
    Propagator,
    SolarPressure,
    SolidTides,
    TidalPerturber,
)
from nyx_space.monte_carlo import MvnSpacecraft, StateDispersion, StateParameter
from nyx_space.time import Duration, Epoch, Unit


def test_howto_propagate_with_perturbations():
    """
    Goal: Propagate a spacecraft's orbit over time while accounting for
    Solar Radiation Pressure (SRP) and atmospheric drag.

    IMPORTANT: This demonstrates SIMPLE propagation segments. It is NOT representative of the sequencing
    available in the premium package of Nyx. Refer to the website for these cases.
    """

    # Step 1: Define Integrator Options
    # We default to a variable step integrator. This ensures numerical accuracy
    # without wasting CPU cycles on regions of the orbit where dynamics change slowly.
    opts = IntegratorOptions()
    logging.info(f"Integrator configured with variable steps: {opts}")
    # TEST Printing the options confirm its variable-step nature
    assert str(opts) == "min_step: 1 ms, max_step: 45 min, tol: 1e-12, attempts: 50"
    # TEST Fixing the min and max step to the same value will force this to be a fixed step integrator
    assert str(IntegratorOptions(Unit.Second * 1, Unit.Second * 1)) == "fixed step: 1 s"

    # Step 2: Define the Environment
    # Define which celestial bodies exert a point-mass gravitational pull.
    accel_models = AccelModels(
        point_masses=PointMasses(
            celestial_objects=[
                CelestialObjects.EARTH,
                CelestialObjects.MOON,
                CelestialObjects.SUN,
            ],
        ),
    )

    # Load planetary ephemerides.
    almanac = (
        MetaAlmanac("../data/02_config/ci_almanac.dhall")
        .process()
        .load("../data/01_planetary/earth_2025_250826_2125_predict.bpc")
    )

    # Step 3: Define Non-Gravitational Forces
    # Configure Solar Radiation Pressure targeting Earth and Moon frames.
    srp = SolarPressure(
        [Frames.EARTH_J2000, Frames.MOON_J2000], almanac, estimate=False
    )
    # Configure atmospheric drag using a standard exponential atmospheric model.
    drag = Drag(
        AtmDensity.earth_exponential(),
        almanac.frame_info(Frames.IAU_EARTH_FRAME),
        estimate=False,
    )

    # Combine accelerations and forces into a unified Dynamics model.
    dynamics = Dynamics(accel_models, force_models=ForceModels(srp, drag))

    # Step 4: Initialize the Propagator
    # The Propagator encapsulates the physics, the math (RungeKutta89), and the universe (almanac).
    # It remains stateless regarding the spacecraft, allowing it to be reused.
    propagator = Propagator(dynamics, almanac, IntegratorMethod.RungeKutta89, opts)

    # Step 5: Define the Spacecraft Initial State
    # We define the starting Keplerian orbital elements.
    eme2k = almanac.frame_info(Frames.EME2000)
    orbit = Orbit.from_keplerian(
        sma_km=6800.0,
        ecc=1e-4,
        inc_deg=45.0,
        raan_deg=60.0,
        aop_deg=75.0,
        ta_deg=90.0,
        epoch=Epoch("2030-12-29 01:02:03 TDB"),
        frame=eme2k,
    )

    # Define the physical characteristics of the spacecraft necessary to compute drag and SRP.
    # Surfaces are in square meters, masses are in kilograms.
    my_sc = Spacecraft(
        orbit,
        Mass(dry_mass_kg=123.0),
        SRPData(area_m2=10.0, coeff_reflectivity=1.2),
        DragData(area_m2=10.0, coeff_drag=2.0),
    )

    # Compute the acceleration imparted on a spacecraft provided the dynamics, in the reference frame
    # of spacecraft's orbit
    accel_km_s2 = propagator.accel_km_s2(my_sc)
    expect_m_s2 = 3.4576  # Purely a concidental number!
    assert abs(sum(accel_km_s2) * 1e3 / 3 - expect_m_s2) < 1e-3, (
        "wrong computed acceleration"
    )

    # Step 6: Execute Propagations
    # Scenario A: Propagate until a specific temporal epoch.
    target_epoch = Epoch("2030-12-30 01:02:03 UTC")
    result_epoch = propagator.until_epoch(my_sc, target_epoch, trajectory=True)

    # The trajectory may be converted to an ANISE Ephemeris, e.g. for serialization to a BSP/OEM file.
    ephem = result_epoch.trajectory.to_ephemeris("MySpacecraftTraj", ExportCfg())
    logging.info(
        f"Generated {ephem}, ending at {result_epoch.state.orbit.ta_deg():.3f} deg True Anomaly."
    )
    true_anomalies = [(record.orbit.epoch, record.orbit.ta_deg()) for record in ephem]
    # TEST true anomalies progress in time because the ephem iterator steps through the states chronologically.
    prev_ta_deg = true_anomalies[0][1]
    for epoch, ta_deg in true_anomalies[1:10]:
        # Quick test -- all of the validation is done in Rust.
        assert ta_deg > prev_ta_deg
        prev_ta_deg = ta_deg

    # Scenario B: Propagate until a specific orbital event (apoapsis).
    # This demonstrates the root-finding capabilities of Nyx to stop exactly at an event,
    # bounded by a maximum duration to prevent infinite searching.
    result_event = propagator.until_event(
        my_sc, Event.apoapsis(), max_duration=Duration("1 day")
    )
    logging.info(
        f"Apoapsis found {result_event.state.orbit.epoch - orbit.epoch} after initial state."
    )


def test_howto_configure_spherical_harmonics():
    """
    Goal: Configure a high-fidelity spherical harmonic gravity field model
    for precise near-Earth propagation.
    """

    almanac = (
        MetaAlmanac("../data/02_config/ci_almanac.dhall")
        .process()
        .load("../data/01_planetary/earth_2025_250826_2125_predict.bpc")
    )
    eme2k = almanac.frame_info(Frames.EME2000)

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

    # Step 1: Define the Gravity Field Configuration
    # We apply the EGM2008 gravity model. The degree and order dictate the resolution
    # of the spherical harmonics. A 10x10 field captures primary planetary oblateness (J2)
    # and significant localized anomalies. The frame must be Earth-fixed (e.g., ITRF93) for the
    # field to rotate with the planet.
    accel_models = AccelModels(
        point_masses=PointMasses(
            celestial_objects=[CelestialObjects.MOON, CelestialObjects.SUN]
        ),
        gravity_field=GravityFieldConfig(
            degree=10,
            order=10,
            filepath="../data/01_planetary/EGM2008_to2190_TideFree.gz",
            frame=Frames.EARTH_ITRF93.to_frameuid(),
        ),
    )

    # Step 2: Propagate
    # For high-precision gravity fields, Dormand-Prince 7-8 is a robust choice.
    dynamics = Dynamics(accel_models)
    propagator = Propagator(dynamics, almanac, IntegratorMethod.DormandPrince78)

    final_state = propagator.for_duration(spacecraft, Unit.Day * 1)
    logging.info(
        f"Final position magnitude after 1 day: {final_state.state.orbit.rmag_km():.3f} km"
    )

    # TEST Case built in collaboration with IOAstro the .NET astrodynamics tool.
    expected = Orbit(
        -5276.159136,
        4263.301774,
        -404.522560,
        -2.724561,
        -3.933753,
        -5.983828,
        Epoch("2025-08-26 11:55:44 UTC"),
        eme2k,
    )
    ric_diff = expected.ric_difference(final_state.state.orbit)

    pos_err_m = ric_diff.rmag_km() * 1e3
    vel_err_m_s = ric_diff.vmag_km_s() * 1e3

    print(
        f"== RIC diff ==\nPosition: {pos_err_m:.1f} m\nVelocity: {vel_err_m_s:.1e} m/s"
    )
    assert pos_err_m < 27.0
    assert vel_err_m_s < 0.05


def test_howto_execute_simple_monte_carlo():
    """
    Goal: Generate a multi-variate dispersion of initial orbital states
    and propagate them in parallel to analyze trajectory divergence.

    While the Premium Nyx provides a lot more thorough Monte Carlo capabilities
    these dispersed states may be used for single segment propagation with dispersions.
    """

    almanac = MetaAlmanac("../data/02_config/ci_almanac.dhall").process()
    eme2k = almanac.frame_info(Frames.EME2000)

    orbit = Orbit.from_keplerian(
        sma_km=6800.0,
        ecc=1e-8,
        inc_deg=45.0,
        raan_deg=60.0,
        aop_deg=75.0,
        ta_deg=90.0,
        epoch=Epoch("2020-02-29 01:02:03 TDB"),
        frame=eme2k,
    )

    spacecraft = Spacecraft(
        orbit,
        Mass(dry_mass_kg=123.0, prop_mass_kg=12.3, extra_mass_kg=1.0),
        SRPData(area_m2=10.0, coeff_reflectivity=1.2),
        DragData(area_m2=10.0, coeff_drag=2.0),
    )

    # TEST Accessibility of each spacecraft field
    assert spacecraft.mass.total_mass_kg() == 123.0 + 12.3 + 1.0

    # Step 1: Define the Multivariate Normal Spacecraft, which allows dispersions using orbital elements.
    # NOTE This approach transforms the orbital element dispersions into the Cartesian state space
    # using the partial derivatives of the chosen orbital elements with respect to each Cartesian element.
    # All dispersions are applied simultaneously: this ensures the dispersions are properly multi-variate.
    disp = [
        StateDispersion.zero_mean(
            StateParameter.Element(OrbitalElement.SemiMajorAxis), 15.0
        ),
        # We disperse the eccentricity by orders of magnitude more than its nominal value (1e-6 vs 1e-8).
        # Because this is a rigorous physical transformation, negative eccentricity magnitudes
        # computed by the normal distribution will be strictly bounded to positive values,
        # mathematically necessitating a 180-degree flip in the Argument of Periapsis.
        StateDispersion.zero_mean(
            StateParameter.Element(OrbitalElement.Eccentricity), 1e-6
        ),
    ]

    # Step 2: Sample the Dispersions
    # Generate a requested number of dispersed spacecraft variants based on the defined distributions.
    mvn = MvnSpacecraft(spacecraft, disp)
    dispersed_spacecraft = mvn.sample(100, seed=123)
    dispersed_spacecraft_chk_seed = mvn.sample(100, seed=123)
    assert dispersed_spacecraft_chk_seed == dispersed_spacecraft, "seed not applied"

    # TEST Verify that the AoP is now a bimodal distribution
    dispersion_df = pl.DataFrame(
        {
            "SMA (km)": [s.orbit.sma_km() for s in dispersed_spacecraft],
            "Ecc": [s.orbit.ecc() for s in dispersed_spacecraft],
            "AoP (deg)": [s.orbit.aop_deg() for s in dispersed_spacecraft],
        }
    )

    median_aop_deg = dispersion_df["AoP (deg)"].quantile(0.5)
    # Check the 1th quantile above the median
    assert (
        dispersion_df.filter(pl.col("AoP (deg)") > median_aop_deg)[
            "AoP (deg)"
        ].quantile(0.1)
        > 344.0
    )
    # Check the 1th quantil below the median
    assert (
        165.0
        < dispersion_df.filter(pl.col("AoP (deg)") < median_aop_deg)[
            "AoP (deg)"
        ].quantile(0.9)
        < 166.0
    )

    # Step 3: Configure the propagator
    accel_models = AccelModels(
        point_masses=PointMasses(
            celestial_objects=[CelestialObjects.EARTH, CelestialObjects.MOON]
        ),
    )
    srp = SolarPressure([Frames.EARTH_J2000], almanac, estimate=False)
    dynamics = Dynamics(accel_models, force_models=ForceModels(srp))
    propagator = Propagator(dynamics, almanac, IntegratorMethod.RungeKutta89)

    # Step 4: Execute the single segment Monte Carlo
    # Calling `many_until_event` delegates the workload directly to the underlying Rust engine.
    # This intentionally bypasses the Python Global Interpreter Lock (GIL), saturating all
    # available CPU cores with parallel propagations.
    # We propagate each spacecraft until it hits its 25th occurrence of Local Solar Time equaling 6.0 hours.
    logging.info(
        f"Initiating parallel propagation of {len(dispersed_spacecraft)} dispersed states."
    )
    results = propagator.many_until_event(
        dispersed_spacecraft,
        Event(ScalarExpr.LocalSolarTime(), Condition.Equals(6.0), Duration("1 ms")),
        trajectory=False,
        trigger=25,
        max_duration=Duration("1 day"),
    )

    # Step 5: Analyze Results
    # Collect the distinct epochs and Longitude of the Ascending Node (LTAN) at the stopping condition.
    df_data = {"Epoch (UTC)": [], "LTAN (deg)": []}
    for rslt in results:
        df_data["Epoch (UTC)"] += [rslt.state.orbit.epoch.to_datetime()]
        df_data["LTAN (deg)"] += [rslt.state.orbit.ltan_deg()]

    df = pl.DataFrame(df_data)
    print(df.describe())
    # We can confirm that there were dispersions because all of the epochs for the 25th event were unique
    assert df["Epoch (UTC)"].unique().count() == len(df)
    assert (df["Epoch (UTC)"].max() - df["Epoch (UTC)"].min()).total_seconds() > 0.0
    # And the final LTAN isn't strictly equal
    assert 59.9 < df["LTAN (deg)"].quantile(0.1) < 60.0
    assert 60.0 < df["LTAN (deg)"].quantile(0.9) < 60.1


def test_howto_configure_solid_tides():
    """
    Goal: Configure the high-fidelity solid tides model for the Cislunar system.
    """

    almanac = MetaAlmanac("../data/02_config/ci_almanac.dhall").process()
    moon_j2k = almanac.frame_info(Frames.MOON_J2000)

    iau_earth = almanac.frame_info(Frames.IAU_EARTH_FRAME)
    iau_moon = almanac.frame_info(Frames.IAU_MOON_FRAME)

    # Define a highly elliptical Earth rbiter at its apoapse
    orbit = Orbit.from_keplerian(
        sma_km=10_000.0,
        ecc=0.7,
        inc_deg=35.0,
        raan_deg=45.0,
        aop_deg=65.0,
        ta_deg=180.0,
        epoch=Epoch("2025-08-25 11:55:44 UTC"),
        frame=almanac.frame_info(Frames.EARTH_J2000),
    )
    spacecraft = Spacecraft(orbit)

    # Step 1: Define the solid tides.
    # IMPORTANT: Solid tides are defined in the body fixed frame, so we must initialize the model with these frames.
    solid_tides = SolidTides.earth_moon_system(
        Frames.IAU_EARTH_FRAME, moon_frame=Frames.IAU_MOON_FRAME, almanac=almanac
    )
    print(solid_tides)

    # NOTE Nyx's Solid Tides modeling is fully configurable, allowing for its use in gas giant systems (e.g. Jupiter).
    # The following is a demonstration of configuring the solid tides model of the Earth/Moon system using
    # the exact tidal perturbers. Importantly, this does not initialize the frame data for the user, so the frames
    # MUST be apriori loaded with their gravitational parameters.
    iau_earth = almanac.frame_info(Frames.IAU_EARTH_FRAME)
    iau_moon = almanac.frame_info(Frames.IAU_MOON_FRAME)
    sun = almanac.frame_info(Frames.SUN_J2000)
    solid_tides_detailed = SolidTides(
        frame=iau_earth,
        k2=0.3019,
        k3=0.093,
        perturbers=[
            TidalPerturber(iau_moon, compute_degree_3=True),
            TidalPerturber(sun, False),
        ],
    )

    assert solid_tides_detailed == solid_tides

    # Define the acceleration mode and finalize the dynamics
    accel_models = AccelModels(
        point_masses=PointMasses(
            celestial_objects=[CelestialObjects.MOON, CelestialObjects.SUN]
        ),
        solid_tides=solid_tides,
    )
    dynamics = Dynamics(accel_models)

    # Step 2: Propagate
    # For high-precision gravity fields, Dormand-Prince 7-8 is a robust choice.
    propagator = Propagator(dynamics, almanac, IntegratorMethod.DormandPrince78)

    final_state = propagator.for_duration(spacecraft, Unit.Day * 1, trajectory=False)
    logging.info(
        f"Final position magnitude: {final_state.state.orbit.rmag_km():.3f} km"
    )

    # TEST Propagate without the solid tides, ceteris paribus.
    # Simply copy the previous accels and unset the solid tides
    third_body_accels = accel_models
    third_body_accels.solid_tides = None

    dynamics = Dynamics(third_body_accels)
    final_state_no_st = Propagator(
        dynamics, almanac, IntegratorMethod.DormandPrince78
    ).for_duration(spacecraft, Unit.Day * 1, trajectory=False)

    ric_diff = final_state.state.orbit.ric_difference(final_state_no_st.state.orbit)

    pos_err_m = ric_diff.rmag_km() * 1e3
    vel_err_m_s = ric_diff.vmag_km_s() * 1e3

    print(
        f"== RIC diff with/without cislunar solid tides ==\nPosition: {pos_err_m:.1f} m\nVelocity: {vel_err_m_s:.1e} m/s"
    )

    assert pos_err_m > 0.0
    assert vel_err_m_s > 0.0


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)

    test_howto_propagate_with_perturbations()
    test_howto_execute_simple_monte_carlo()
    test_howto_configure_spherical_harmonics()
    test_howto_configure_solid_tides()
