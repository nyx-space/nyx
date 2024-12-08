#![doc = include_str!("./README.md")]
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::{
    almanac::metaload::MetaFile,
    constants::{
        celestial_objects::{MOON, SUN},
        frames::{EARTH_J2000, IAU_EARTH_FRAME},
    },
};
use hifitime::{Epoch, TimeUnits, Unit};
use nyx::{
    cosmic::{eclipse::EclipseLocator, GuidanceMode, MetaAlmanac, Orbit, SrpConfig},
    dynamics::{
        guidance::{Ruggiero, Thruster},
        Harmonics, OrbitalDynamics, SolarPressure, SpacecraftDynamics,
    },
    io::{gravity::HarmonicsMem, ExportCfg},
    mc::{MonteCarlo, MvnSpacecraft, StateDispersion},
    md::{prelude::Objective, StateParameter},
    propagators::{ErrorControl, IntegratorOptions, Propagator},
    Spacecraft, State,
};
use std::{error::Error, sync::Arc};

fn main() -> Result<(), Box<dyn Error>> {
    pel::init();
    // Set up the dynamics like in the orbit raise.
    let almanac = Arc::new(MetaAlmanac::latest().map_err(Box::new)?);
    let epoch = Epoch::from_gregorian_utc_hms(2024, 2, 29, 12, 13, 14);

    // Define the GEO orbit, and we're just going to maintain it very tightly.
    let earth_j2000 = almanac.frame_from_uid(EARTH_J2000)?;
    let orbit = Orbit::try_keplerian(42164.0, 1e-5, 0., 163.0, 75.0, 0.0, epoch, earth_j2000)?;
    println!("{orbit:x}");

    let sc = Spacecraft::builder()
        .orbit(orbit)
        .dry_mass_kg(1000.0) // 1000 kg of dry mass
        .fuel_mass_kg(1000.0) // 1000 kg of fuel, totalling 2.0 tons
        .srp(SrpConfig::from_area(3.0 * 6.0)) // Assuming 1 kW/m^2 or 18 kW, giving a margin of 4.35 kW for on-propulsion consumption
        .thruster(Thruster {
            // "NEXT-STEP" row in Table 2
            isp_s: 4435.0,
            thrust_N: 0.472,
        })
        .mode(GuidanceMode::Thrust) // Start thrusting immediately.
        .build();

    // Set up the spacecraft dynamics like in the orbit raise example.

    let prop_time = 30.0 * Unit::Day;

    // Define the guidance law -- we're just using a Ruggiero controller as demonstrated in AAS-2004-5089.
    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42_164.0, 5.0), // 5 km
        Objective::within_tolerance(StateParameter::Eccentricity, 0.001, 5e-5),
        Objective::within_tolerance(StateParameter::Inclination, 0.05, 1e-2),
    ];

    let ruggiero_ctrl = Ruggiero::from_max_eclipse(objectives, sc, 0.2)?;
    println!("{ruggiero_ctrl}");

    let mut orbital_dyn = OrbitalDynamics::point_masses(vec![MOON, SUN]);

    let mut jgm3_meta = MetaFile {
        uri: "http://public-data.nyxspace.com/nyx/models/JGM3.cof.gz".to_string(),
        crc32: Some(0xF446F027), // Specifying the CRC32 avoids redownloading it if it's cached.
    };
    jgm3_meta.process(true)?;

    let harmonics = Harmonics::from_stor(
        almanac.frame_from_uid(IAU_EARTH_FRAME)?,
        HarmonicsMem::from_cof(&jgm3_meta.uri, 8, 8, true)?,
    );
    orbital_dyn.accel_models.push(harmonics);

    let srp_dyn = SolarPressure::default(EARTH_J2000, almanac.clone())?;
    let sc_dynamics = SpacecraftDynamics::from_model(orbital_dyn, srp_dyn)
        .with_guidance_law(ruggiero_ctrl.clone());

    println!("{sc_dynamics}");

    // Finally, let's use the Monte Carlo framework built into Nyx to propagate spacecraft.

    // Let's start by defining the dispersion.
    // The MultivariateNormal structure allows us to define the dispersions in any of the orbital parameters, but these are applied directly in the Cartesian state space.
    // Note that additional validation on the MVN is in progress -- https://github.com/nyx-space/nyx/issues/339.
    let mc_rv = MvnSpacecraft::new(
        sc,
        vec![StateDispersion::zero_mean(StateParameter::SMA, 3.0)],
    )?;

    let my_mc = MonteCarlo::new(
        sc, // Nominal state
        mc_rv,
        "03_geo_sk".to_string(), // Scenario name
        None, // No specific seed specified, so one will be drawn from the computer's entropy.
    );

    // Build the propagator setup.
    let setup = Propagator::rk89(
        sc_dynamics.clone(),
        IntegratorOptions::builder()
            .min_step(10.0_f64.seconds())
            .error_ctrl(ErrorControl::RSSCartesianStep)
            .build(),
    );

    let num_runs = 25;
    let rslts = my_mc.run_until_epoch(setup, almanac.clone(), sc.epoch() + prop_time, num_runs);

    assert_eq!(rslts.runs.len(), num_runs);

    // For all of the resulting trajectories, we'll want to compute the percentage of penumbra and umbra.

    rslts.to_parquet(
        "03_geo_sk.parquet",
        Some(vec![
            &EclipseLocator::cislunar(almanac.clone()).to_penumbra_event()
        ]),
        ExportCfg::default(),
        almanac,
    )?;

    Ok(())
}
