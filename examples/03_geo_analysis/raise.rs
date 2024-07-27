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
use hifitime::{Epoch, Unit};
use nyx::{
    cosmic::{GuidanceMode, MetaAlmanac, Orbit, SrpConfig},
    dynamics::{
        guidance::{GuidanceLaw, Ruggiero, Thruster},
        Harmonics, OrbitalDynamics, SolarPressure, SpacecraftDynamics,
    },
    io::{gravity::HarmonicsMem, ExportCfg},
    md::{prelude::Objective, StateParameter},
    propagators::Propagator,
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
    let eme2k = almanac.frame_from_uid(EARTH_J2000).unwrap();
    // Define the orbit epoch
    let epoch = Epoch::from_gregorian_utc_hms(2024, 2, 29, 12, 13, 14);

    // Build the spacecraft itself.
    // Using slide 6 of https://aerospace.org/sites/default/files/2018-11/Davis-Mayberry_HPSEP_11212018.pdf
    // for the "next gen" SEP characteristics.

    // GTO start
    let orbit = Orbit::keplerian(24505.9, 0.725, 7.05, 0.0, 0.0, 0.0, epoch, eme2k);

    let sc = Spacecraft::builder()
        .orbit(orbit)
        .dry_mass_kg(1000.0) // 1000 kg of dry mass
        .fuel_mass_kg(1500.0) // 1500 kg of fuel, totalling 2.5 tons
        .srp(SrpConfig::from_area(3.0 * 6.0)) // Assuming 1 kW/m^2 or 18 kW, giving a margin of 4.35 kW for on-propulsion consumption
        .thruster(Thruster {
            isp_s: 4435.0,
            thrust_N: 0.472,
        }) // "NEXT-STEP" row in Table 2
        .mode(GuidanceMode::Thrust) // Start thrusting immediately.
        .build();

    let prop_time = 120.0 * Unit::Day;

    // Define the guidance law -- we're just using a Ruggiero controller as demonstrated in AAS-2004-5089.
    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42_165.0, 20.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.001, 5e-5),
        // Objective::within_tolerance(StateParameter::Inclination, 0.05, 5e-3),
    ];

    // Define the efficiency thresholds for this controller
    let ηthresholds = [0.5, 0.75, 0.85];
    let ruggiero_ctrl = Ruggiero::with_ηthresholds(objectives, &ηthresholds, sc).unwrap();
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
    jgm3_meta.process()?;

    // Build the spherical harmonics.
    // The harmonics must be computed in the body fixed frame.
    // We're using the long term prediction of the Earth centered Earth fixed frame, IAU Earth.
    let harmonics = Harmonics::from_stor(
        almanac.frame_from_uid(IAU_EARTH_FRAME)?,
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

    println!("[qlaw_as_ruggiero_case_b] {:x}", orbit);

    let (final_state, traj) = Propagator::default(
        sc_dynamics.clone(),
        // PropOpts::<RSSCartesianStep>::with_tolerance(1e-10),
    )
    .with(sc, almanac.clone())
    .for_duration_with_traj(prop_time)?;

    let fuel_usage = sc.fuel_mass_kg - final_state.fuel_mass_kg;
    println!("[qlaw_as_ruggiero_case_b] {:x}", final_state.orbit);
    println!("[qlaw_as_ruggiero_case_b] fuel usage: {:.3} kg", fuel_usage);

    // Finally, export the results for analysis.
    traj.to_parquet_with_cfg("./03_geo_raise.parquet", ExportCfg::default(), almanac)?;

    for status_line in ruggiero_ctrl.status(&final_state) {
        println!("{status_line}");
    }

    assert!(
        ruggiero_ctrl.achieved(&final_state).unwrap(),
        "objective not achieved"
    );

    Ok(())
}
