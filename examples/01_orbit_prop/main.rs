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
use log::warn;
use nyx::{
    cosmic::{MetaAlmanac, Orbit, SrpConfig},
    dynamics::{Harmonics, OrbitalDynamics, SolarPressure, SpacecraftDynamics},
    io::{gravity::HarmonicsMem, ExportCfg},
    od::GroundStation,
    propagators::Propagator,
    Spacecraft, State,
};
use polars::{df, series::ChunkCompare};

use std::{error::Error, sync::Arc};

fn main() -> Result<(), Box<dyn Error>> {
    pel::init();
    // Dynamics models require planetary constants and ephemerides to be defined.
    // Let's start by grabbing those by using ANISE's latest MetaAlmanac.
    // This will automatically download the DE440s (or later) planetary ephemeris,
    // the daily-updated Earth Orientation Parameters, the high fidelity Moon orientation
    // parameters (for the Moon Mean Earth and Moon Principal Axes frames), and the PCK11
    // planetary constants kernels.
    // For details, refer to https://github.com/nyx-space/anise/blob/master/data/latest.dhall.
    // Note that we place the Almanac into an Arc so we can clone it cheaply across threads.
    let almanac = Arc::new(MetaAlmanac::latest().map_err(Box::new)?);
    // Define the orbit epoch
    let epoch = Epoch::from_gregorian_utc_hms(2024, 2, 29, 12, 13, 14);

    // Define the orbit in Keplerian because it's easy to visualize.
    // First we need to fetch the Earth J2000 from information from the Almanac.
    // This allows the frame to include the gravitational parameters and the shape of the Earth,
    // defined as a tri-axial ellipoid. Note that this shape can be changed manually or in the Almanac
    // by loading a different set of planetary constants.
    let earth_j2000 = almanac.frame_from_uid(EARTH_J2000)?;

    // Define the orbit from its altitude.
    let orbit =
        Orbit::try_keplerian_altitude(300.0, 0.015, 28.5, 65.2, 75.0, 0.0, epoch, earth_j2000)?;
    // Print in in Keplerian format
    println!("{orbit:x}");

    // There are two ways to propagate an orbit. We can make a quick approximation assuming only two-body
    // motion. This is a useful first order approximation but it isn't used in real-world applications.

    // This approach is a feature of ANISE.
    let future_orbit_tb = orbit.at_epoch(epoch + Unit::Day * 3)?;
    println!("{future_orbit_tb:x}");

    // Two body propagation relies solely on Kepler's laws, so none of the Keplerian orbital elements should
    // have changed apart from the true anomaly.

    println!(
        "SMA changed by {:.3e} km",
        orbit.sma_km()? - future_orbit_tb.sma_km()?
    );
    println!(
        "ECC changed by {:.3e}",
        orbit.ecc()? - future_orbit_tb.ecc()?
    );
    println!(
        "INC changed by {:.3e} deg",
        orbit.inc_deg()? - future_orbit_tb.inc_deg()?
    );
    println!(
        "RAAN changed by {:.3e} deg",
        orbit.raan_deg()? - future_orbit_tb.raan_deg()?
    );
    println!(
        "AOP changed by {:.3e} deg",
        orbit.aop_deg()? - future_orbit_tb.aop_deg()?
    );
    println!(
        "TA changed by {:.3} deg",
        orbit.ta_deg()? - future_orbit_tb.ta_deg()?
    );

    // Nyx is used for high fidelity propagation, not Keplerian propagation as above.
    // Nyx only propagates Spacecraft at the moment, which allows it to account for acceleration
    // models such as solar radiation pressure.

    // Let's build a cubesat sized spacecraft, with an SRP area of 10 cm^2 and a mass of 9.6 kg.
    let sc = Spacecraft::builder()
        .orbit(orbit)
        .dry_mass_kg(9.60)
        .srp(SrpConfig {
            area_m2: 100e-4,
            cr: 1.1,
        })
        .build();
    println!("{sc:x}");

    // Let's set up the spacecraft dynamics.

    // Specify that the orbital dynamics must account for the graviational pull of the Moon and the Sun.
    // The gravity of the Earth will also be accounted for since the spaceraft in an Earth orbit.
    let mut orbital_dyn = OrbitalDynamics::point_masses(vec![MOON, SUN]);

    // We want to include the spherical harmonics, so let's download the gravitational data from the Nyx Cloud.
    // We're using the JGM3 model here because it's pretty lightweight but it is the default in GMAT.
    let mut jgm3_meta = MetaFile {
        uri: "http://public-data.nyxspace.com/nyx/models/JGM3.cof.gz".to_string(),
        crc32: Some(0xF446F027), // Specifying the CRC32 avoids redownloading it if it's cached.
    };
    // And let's download it if we don't have it yet.
    jgm3_meta.process()?;

    // Now let's build the spherical harmonics.
    // The harmonics must be computed in the body fixed frame. We're using the long term prediction of the Earth body frame, IAU Earth.
    let harmonics_21x21 = Harmonics::from_stor(
        almanac.frame_from_uid(IAU_EARTH_FRAME)?,
        HarmonicsMem::from_cof(&jgm3_meta.uri, 21, 21, true).unwrap(),
    );

    // Include the spherical harmonics into the orbital dynamics.
    orbital_dyn.accel_models.push(harmonics_21x21);

    // We define the solar radiation pressure, using the default solar flux and accounting only
    // for the eclipsing caused by the Earth.
    let srp_dyn = SolarPressure::default(EARTH_J2000, almanac.clone())?;

    let dynamics = SpacecraftDynamics::from_model(orbital_dyn, srp_dyn);

    println!("{dynamics}");

    // Finally, let's propagate this orbit to the same epoch as above, and compute the differences.
    let (future_sc, trajectory) = Propagator::default(dynamics)
        .with(sc, almanac.clone())
        .until_epoch_with_traj(future_orbit_tb.epoch)?;

    println!("=== High fidelity propagation ===");
    println!(
        "SMA changed by {:.3e} km",
        orbit.sma_km()? - future_sc.orbit.sma_km()?
    );
    println!(
        "ECC changed by {:.3e}",
        orbit.ecc()? - future_sc.orbit.ecc()?
    );
    println!(
        "INC changed by {:.3e} deg",
        orbit.inc_deg()? - future_sc.orbit.inc_deg()?
    );
    println!(
        "RAAN changed by {:.3e} deg",
        orbit.raan_deg()? - future_sc.orbit.raan_deg()?
    );
    println!(
        "AOP changed by {:.3e} deg",
        orbit.aop_deg()? - future_sc.orbit.aop_deg()?
    );
    println!(
        "TA changed by {:.3} deg",
        orbit.ta_deg()? - future_sc.orbit.ta_deg()?
    );

    // We also have access to the full trajectory throughout the propagation.
    println!("{trajectory}");

    // With the trajectory, let's build a few data products (we want to analyze something after all).

    // 1. Export the trajectory as a CCSDS OEM version 2.0 file and as a parquet file, which includes the Keplerian orbital elements.

    trajectory.to_oem_file(
        "./01_cubesat_hf_prop.oem",
        ExportCfg::builder().step(Unit::Minute * 2).build(),
    )?;

    trajectory.to_parquet_with_cfg(
        "./01_cubesat_hf_prop.parquet",
        ExportCfg::builder().step(Unit::Minute * 2).build(),
        almanac.clone(),
    )?;

    // 2. Compare the difference in the radial-intrack-crosstrack frame between the high fidelity
    // and Keplerian propagation. The RIC frame is commonly used to compute the difference in position
    // and velocity of different spacecraft.
    // 3. Compute the azimuth, elevation, range, and range-rate data of that spacecraft as seen from Boulder, CO, USA.

    let boulder_station = GroundStation::from_point(
        "Boulder, CO, USA".to_string(),
        40.014984,
        -105.270546,
        1655.0,
        almanac.frame_from_uid(IAU_EARTH_FRAME)?,
    );

    // We iterate over the trajectory, grabbing a state every two minutes.
    let mut offset_s = vec![];
    let mut ric_x_km = vec![];
    let mut ric_y_km = vec![];
    let mut ric_z_km = vec![];
    let mut ric_vx_km_s = vec![];
    let mut ric_vy_km_s = vec![];
    let mut ric_vz_km_s = vec![];

    let mut azimuth_deg = vec![];
    let mut elevation_deg = vec![];
    let mut range_km = vec![];
    let mut range_rate_km_s = vec![];
    for state in trajectory.every(Unit::Minute * 2) {
        // Try to compute the Keplerian/two body state just in time.
        // This method sometimes fails to converge on an appropriate true anomaly
        // from the mean anomaly. If that happens, we just skip this state.
        // The high fidelity and Keplerian states diverge continuously, and we're curious
        // about the divergence.
        let this_epoch = state.epoch();
        match orbit.at_epoch(this_epoch) {
            Ok(tb_then) => {
                offset_s.push((this_epoch - orbit.epoch).to_seconds());
                // Compute the two body state just in time.
                let ric = state.orbit.ric_difference(&tb_then)?;
                ric_x_km.push(ric.radius_km.x);
                ric_y_km.push(ric.radius_km.y);
                ric_z_km.push(ric.radius_km.z);
                ric_vx_km_s.push(ric.velocity_km_s.x);
                ric_vy_km_s.push(ric.velocity_km_s.y);
                ric_vz_km_s.push(ric.velocity_km_s.z);

                // Compute the AER data
                let aer = almanac.azimuth_elevation_range_sez(
                    boulder_station.to_orbit(this_epoch, &almanac)?,
                    state.orbit,
                )?;
                azimuth_deg.push(aer.azimuth_deg);
                elevation_deg.push(aer.elevation_deg);
                range_km.push(aer.range_km);
                range_rate_km_s.push(aer.range_rate_km_s);
            }
            Err(e) => warn!("{} {e}", state.epoch()),
        };
    }

    let ric_df = df!(
        "Offset (s)" => offset_s.clone(),
        "RIC X (km)" => ric_x_km,
        "RIC Y (km)" => ric_y_km,
        "RIC Z (km)" => ric_z_km,
        "RIC VX (km/s)" => ric_vx_km_s,
        "RIC VY (km/s)" => ric_vy_km_s,
        "RIC VZ (km/s)" => ric_vz_km_s,
    )?;

    println!("RIC difference at start\n{}", ric_df.head(Some(10)));
    println!("RIC difference at end\n{}", ric_df.tail(Some(10)));

    let aer_df = df!("Offset (s)" => offset_s.clone(),
        "azimuth (deg)" => azimuth_deg,
        "elevation (deg)" => elevation_deg,
        "range (km)" => range_km,
        "range rate (km/s)" => range_rate_km_s,
    )?;

    // Finally, let's see when the spacecraft is visible, assuming 15 degrees minimum elevation.
    let mask = aer_df.column("elevation (deg)")?.gt(15.0)?;
    let cubesat_visible = aer_df.filter(&mask)?;

    println!("{cubesat_visible}");

    Ok(())
}
