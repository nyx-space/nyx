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
    cosmic::{eclipse::EclipseLocator, GuidanceMode, Mass, MetaAlmanac, Orbit, SRPData},
    dynamics::{
        guidance::{GuidanceLaw, Ruggiero, Thruster},
        Harmonics, OrbitalDynamics, SolarPressure, SpacecraftDynamics,
    },
    io::{gravity::HarmonicsMem, ExportCfg},
    md::{prelude::Objective, StateParameter},
    propagators::{ErrorControl, IntegratorOptions, Propagator},
    Spacecraft,
};
use std::{error::Error, sync::Arc};

fn main() -> Result<(), Box<dyn Error>> {
    pel::init();

    // Dynamics models require planetary constants and ephemerides to be defined.
    // Let's start by grabbing those by using ANISE's latest MetaAlmanac.
    // This will automatically download the DE440s planetary ephemeris,
    // the daily-updated Earth Orientation Parameters, the high fidelity Moon orientation
    // parameters (for the Moon Mean Earth and Moon Principal Axes frames), and the PCK11
    // planetary constants kernels.
    // For details, refer to https://github.com/nyx-space/anise/blob/master/data/latest.dhall.
    // Note that we place the Almanac into an Arc so we can clone it cheaply and provide read-only
    // references to many functions.
    let almanac = Arc::new(MetaAlmanac::latest().map_err(Box::new)?);
    // Fetch the EME2000 frame from the Almabac
    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
    // Define the orbit epoch
    let epoch = Epoch::from_gregorian_utc_hms(2024, 2, 29, 12, 13, 14);

    // Build the spacecraft itself.
    // Using slide 6 of https://aerospace.org/sites/default/files/2018-11/Davis-Mayberry_HPSEP_11212018.pdf
    // for the "next gen" SEP characteristics.

    // GTO start
    let orbit = Orbit::keplerian(24505.9, 0.725, 7.05, 0.0, 0.0, 0.0, epoch, eme2k);

    let sc = Spacecraft::builder()
        .orbit(orbit)
        .mass(Mass::from_dry_and_prop_masses(1000.0, 1000.0)) // 1000 kg of dry mass and prop, totalling 2.0 tons
        .srp(SRPData::from_area(3.0 * 6.0)) // Assuming 1 kW/m^2 or 18 kW, giving a margin of 4.35 kW for on-propulsion consumption
        .thruster(Thruster {
            // "NEXT-STEP" row in Table 2
            isp_s: 4435.0,
            thrust_N: 0.472,
        })
        .mode(GuidanceMode::Thrust) // Start thrusting immediately.
        .build();

    let prop_time = 180.0 * Unit::Day;

    // Define the guidance law -- we're just using a Ruggiero controller as demonstrated in AAS-2004-5089.
    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42_165.0, 20.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.001, 5e-5),
        Objective::within_tolerance(StateParameter::Inclination, 0.05, 1e-2),
    ];

    // Ensure that we only thrust if we have more than 20% illumination.
    let ruggiero_ctrl = Ruggiero::from_max_eclipse(objectives, sc, 0.2).unwrap();
    println!("{ruggiero_ctrl}");

    // Define the high fidelity dynamics

    // Set up the spacecraft dynamics.

    // Specify that the orbital dynamics must account for the graviational pull of the Moon and the Sun.
    // The gravity of the Earth will also be accounted for since the spaceraft in an Earth orbit.
    let mut orbital_dyn = OrbitalDynamics::point_masses(vec![MOON, SUN]);

    // We want to include the spherical harmonics, so let's download the gravitational data from the Nyx Cloud.
    // We're using the JGM3 model here, which is the default in GMAT.
    let mut jgm3_meta = MetaFile {
        uri: "http://public-data.nyxspace.com/nyx/models/JGM3.cof.gz".to_string(),
        crc32: Some(0xF446F027), // Specifying the CRC32 avoids redownloading it if it's cached.
    };
    // And let's download it if we don't have it yet.
    jgm3_meta.process(true)?;

    // Build the spherical harmonics.
    // The harmonics must be computed in the body fixed frame.
    // We're using the long term prediction of the Earth centered Earth fixed frame, IAU Earth.
    let harmonics = Harmonics::from_stor(
        almanac.frame_info(IAU_EARTH_FRAME)?,
        HarmonicsMem::from_cof(&jgm3_meta.uri, 8, 8, true).unwrap(),
    );

    // Include the spherical harmonics into the orbital dynamics.
    orbital_dyn.accel_models.push(harmonics);

    // We define the solar radiation pressure, using the default solar flux and accounting only
    // for the eclipsing caused by the Earth.
    let srp_dyn = SolarPressure::default(EARTH_J2000, almanac.clone())?;

    // Finalize setting up the dynamics, specifying the force models (orbital_dyn) separately from the
    // acceleration models (SRP in this case). Use `from_models` to specify multiple accel models.
    let sc_dynamics = SpacecraftDynamics::from_model(orbital_dyn, srp_dyn)
        .with_guidance_law(ruggiero_ctrl.clone());

    println!("{orbit:x}");

    // We specify a minimum step in the propagator because the Ruggiero control would otherwise drive this step very low.
    let (final_state, traj) = Propagator::rk89(
        sc_dynamics.clone(),
        IntegratorOptions::builder()
            .min_step(10.0_f64.seconds())
            .error_ctrl(ErrorControl::RSSCartesianStep)
            .build(),
    )
    .with(sc, almanac.clone())
    .for_duration_with_traj(prop_time)?;

    let prop_usage = sc.mass.prop_mass_kg - final_state.mass.prop_mass_kg;
    println!("{:x}", final_state.orbit);
    println!("prop usage: {prop_usage:.3} kg");

    // Finally, export the results for analysis, including the penumbra percentage throughout the orbit raise.
    traj.to_parquet(
        "./03_geo_raise.parquet",
        Some(vec![
            &EclipseLocator::cislunar(almanac.clone()).to_penumbra_event()
        ]),
        ExportCfg::default(),
        almanac,
    )?;

    for status_line in ruggiero_ctrl.status(&final_state) {
        println!("{status_line}");
    }

    ruggiero_ctrl
        .achieved(&final_state)
        .expect("objective not achieved");

    Ok(())
}
