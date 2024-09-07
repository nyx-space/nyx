#![doc = include_str!("./README.md")]
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::{
    almanac::metaload::MetaFile,
    constants::{
        celestial_objects::{MOON, SUN},
        frames::{EARTH_J2000, IAU_EARTH_FRAME, MOON_J2000},
    },
};
use hifitime::{Epoch, Unit};
use nyx::{
    cosmic::{eclipse::EclipseLocator, MetaAlmanac, Orbit, SrpConfig},
    dynamics::{Harmonics, OrbitalDynamics, SolarPressure, SpacecraftDynamics},
    io::{gravity::HarmonicsMem, ExportCfg},
    propagators::Propagator,
    Spacecraft, State,
};
use polars::{df, prelude::ParquetWriter};

use std::fs::File;
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
    // Define the orbit epoch
    let epoch = Epoch::from_gregorian_utc_hms(2024, 2, 29, 12, 13, 14);

    // Define the orbit.
    // First we need to fetch the Earth J2000 from information from the Almanac.
    // This allows the frame to include the gravitational parameters and the shape of the Earth,
    // defined as a tri-axial ellipoid. Note that this shape can be changed manually or in the Almanac
    // by loading a different set of planetary constants.
    let earth_j2000 = almanac.frame_from_uid(EARTH_J2000)?;

    // Placing this GEO bird just above Colorado.
    // In theory, the eccentricity is zero, but in practice, it's about 1e-5 to 1e-6 at best.
    let orbit = Orbit::try_keplerian(42164.0, 1e-5, 0., 163.0, 75.0, 0.0, epoch, earth_j2000)?;
    // Print in in Keplerian form.
    println!("{orbit:x}");

    let state_bf = almanac.transform_to(orbit, IAU_EARTH_FRAME, None)?;
    let (orig_lat_deg, orig_long_deg, orig_alt_km) = state_bf.latlongalt()?;

    // Nyx is used for high fidelity propagation, not Keplerian propagation as above.
    // Nyx only propagates Spacecraft at the moment, which allows it to account for acceleration
    // models such as solar radiation pressure.

    // Let's build a cubesat sized spacecraft, with an SRP area of 10 cm^2 and a mass of 9.6 kg.
    let sc = Spacecraft::builder()
        .orbit(orbit)
        .dry_mass_kg(9.60)
        .srp(SrpConfig {
            area_m2: 10e-4,
            cr: 1.1,
        })
        .build();
    println!("{sc:x}");

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
    let harmonics_21x21 = Harmonics::from_stor(
        almanac.frame_from_uid(IAU_EARTH_FRAME)?,
        HarmonicsMem::from_cof(&jgm3_meta.uri, 21, 21, true).unwrap(),
    );

    // Include the spherical harmonics into the orbital dynamics.
    orbital_dyn.accel_models.push(harmonics_21x21);

    // We define the solar radiation pressure, using the default solar flux and accounting only
    // for the eclipsing caused by the Earth and Moon.
    let srp_dyn = SolarPressure::new(vec![EARTH_J2000, MOON_J2000], almanac.clone())?;

    // Finalize setting up the dynamics, specifying the force models (orbital_dyn) separately from the
    // acceleration models (SRP in this case). Use `from_models` to specify multiple accel models.
    let dynamics = SpacecraftDynamics::from_model(orbital_dyn, srp_dyn);

    println!("{dynamics}");

    // Finally, let's propagate this orbit to the same epoch as above.
    // The first returned value is the spacecraft state at the final epoch.
    // The second value is the full trajectory where the step size is variable step used by the propagator.
    let (future_sc, trajectory) = Propagator::default(dynamics)
        .with(sc, almanac.clone())
        .until_epoch_with_traj(epoch + Unit::Century * 0.03)?;

    println!("=== High fidelity propagation ===");
    println!(
        "SMA changed by {:.3} km",
        orbit.sma_km()? - future_sc.orbit.sma_km()?
    );
    println!(
        "ECC changed by {:.6}",
        orbit.ecc()? - future_sc.orbit.ecc()?
    );
    println!(
        "INC changed by {:.3e} deg",
        orbit.inc_deg()? - future_sc.orbit.inc_deg()?
    );
    println!(
        "RAAN changed by {:.3} deg",
        orbit.raan_deg()? - future_sc.orbit.raan_deg()?
    );
    println!(
        "AOP changed by {:.3} deg",
        orbit.aop_deg()? - future_sc.orbit.aop_deg()?
    );
    println!(
        "TA changed by {:.3} deg",
        orbit.ta_deg()? - future_sc.orbit.ta_deg()?
    );

    // We also have access to the full trajectory throughout the propagation.
    println!("{trajectory}");

    println!("Spacecraft params after 3 years without active control:\n{future_sc:x}");

    // With the trajectory, let's build a few data products.

    // 1. Export the trajectory as a parquet file, which includes the Keplerian orbital elements.

    let analysis_step = Unit::Minute * 5;

    trajectory.to_parquet(
        "./03_geo_hf_prop.parquet",
        Some(vec![
            &EclipseLocator::cislunar(almanac.clone()).to_penumbra_event()
        ]),
        ExportCfg::builder().step(analysis_step).build(),
        almanac.clone(),
    )?;

    // 2. Compute the latitude, longitude, and altitude throughout the trajectory by rotating the spacecraft position into the Earth body fixed frame.

    // We iterate over the trajectory, grabbing a state every two minutes.
    let mut offset_s = vec![];
    let mut epoch_str = vec![];
    let mut longitude_deg = vec![];
    let mut latitude_deg = vec![];
    let mut altitude_km = vec![];

    for state in trajectory.every(analysis_step) {
        // Convert the GEO bird state into the body fixed frame, and keep track of its latitude, longitude, and altitude.
        // These define the GEO stationkeeping box.

        let this_epoch = state.epoch();

        offset_s.push((this_epoch - orbit.epoch).to_seconds());
        epoch_str.push(this_epoch.to_isoformat());

        let state_bf = almanac.transform_to(state.orbit, IAU_EARTH_FRAME, None)?;
        let (lat_deg, long_deg, alt_km) = state_bf.latlongalt()?;
        longitude_deg.push(long_deg);
        latitude_deg.push(lat_deg);
        altitude_km.push(alt_km);
    }

    println!(
        "Longitude changed by {:.3} deg -- Box is 0.1 deg E-W",
        orig_long_deg - longitude_deg.last().unwrap()
    );

    println!(
        "Latitude changed by {:.3} deg -- Box is 0.05 deg N-S",
        orig_lat_deg - latitude_deg.last().unwrap()
    );

    println!(
        "Altitude changed by {:.3} km -- Box is 30 km",
        orig_alt_km - altitude_km.last().unwrap()
    );

    // Build the station keeping data frame.
    let mut sk_df = df!(
        "Offset (s)" => offset_s.clone(),
        "Epoch (UTC)" => epoch_str.clone(),
        "Longitude E-W (deg)" => longitude_deg,
        "Latitude N-S (deg)" => latitude_deg,
        "Altitude (km)" => altitude_km,

    )?;

    // Create a file to write the Parquet to
    let file = File::create("./03_geo_lla.parquet").expect("Could not create file");

    // Create a ParquetWriter and write the DataFrame to the file
    ParquetWriter::new(file).finish(&mut sk_df)?;

    Ok(())
}
