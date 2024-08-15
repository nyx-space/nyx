#![doc = include_str!("./README.md")]
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::{
    almanac::metaload::MetaFile,
    constants::{
        celestial_objects::{JUPITER_BARYCENTER, MOON, SUN},
        frames::{EARTH_J2000, MOON_J2000},
    },
};
use hifitime::{Epoch, TimeUnits, Unit};
use nyx::{
    cosmic::{eclipse::EclipseLocator, Aberration, Frame, MetaAlmanac, SrpConfig},
    dynamics::{guidance::LocalFrame, OrbitalDynamics, SolarPressure, SpacecraftDynamics},
    io::{ConfigRepr, ExportCfg},
    mc::MonteCarlo,
    md::prelude::Traj,
    od::{
        msr::RangeDoppler,
        prelude::{TrackingArcSim, TrkConfig, KF},
        process::SpacecraftUncertainty,
        GroundStation, SpacecraftODProcess,
    },
    propagators::Propagator,
    Orbit, Spacecraft, State,
};

use std::{collections::BTreeMap, error::Error, path::PathBuf, str::FromStr, sync::Arc};

fn main() -> Result<(), Box<dyn Error>> {
    pel::init();
    // Dynamics models require planetary constants and ephemerides to be defined.
    // Let's start by grabbing those by using ANISE's latest MetaAlmanac.
    // For details, refer to https://github.com/nyx-space/anise/blob/master/data/latest.dhall.

    // Download the latest (of time of writing) LRO definitive ephemeris from the public Nyx Space cloud.
    // Note that the original file is in _big endian_ format, and my machine is little endian, so I've used the
    // `bingo` tool from https://naif.jpl.nasa.gov/naif/utilities_PC_Linux_64bit.html to convert the original file
    // to little endian and upload it to the cloud.
    // Refer to https://naif.jpl.nasa.gov/pub/naif/pds/data/lro-l-spice-6-v1.0/lrosp_1000/data/spk/?C=M;O=D for original file.
    let mut lro_def_ephem = MetaFile {
        uri: "http://public-data.nyxspace.com/nyx/examples/lrorg_2023349_2024075_v01_LE.bsp"
            .to_string(),
        crc32: Some(0xe76ce3b5),
    };
    lro_def_ephem.process()?;

    // Load this ephem in the general Almanac we're using for this analysis.
    let almanac = Arc::new(
        MetaAlmanac::latest()
            .map_err(Box::new)?
            .load_from_metafile(lro_def_ephem)?,
    );

    // Orbit determination requires a Trajectory structure, which can be saved as parquet file.
    // In our case, the trajectory comes from the BSP file, so we need to build a Trajectory from the almanac directly.
    // To query the Almanac, we need to build the LRO frame in the J2000 orientation in our case.
    // Inspecting the LRO BSP in the ANISE GUI shows us that NASA has assigned ID -85 to LRO.
    let lro_frame = Frame::from_ephem_j2000(-85);
    // To build the trajectory we need to provide a spacecraft template.
    let sc_template = Spacecraft::builder()
        .dry_mass_kg(1018.0) // Launch masses
        .fuel_mass_kg(900.0)
        .srp(SrpConfig {
            // SRP configuration is arbitrary, but we will be estimating it anyway.
            area_m2: 3.9 * 2.7,
            cr: 0.9,
        })
        .orbit(Orbit::zero(MOON_J2000)) // Setting a zero orbit here because it's just a template
        .build();
    // Now we can build the trajectory from the BSP file.
    // We'll arbitrarily set the tracking arc to 48 hours with a one minute time step.
    let trajectory = Traj::from_bsp(
        lro_frame,
        MOON_J2000,
        almanac.clone(),
        sc_template,
        5.seconds(),
        Some(Epoch::from_str("2024-01-01 00:00:00 UTC").unwrap()),
        Some(Epoch::from_str("2024-01-03 00:00:00 UTC").unwrap()),
        Aberration::LT,
        Some("LRO".to_string()),
    )?;

    println!("{trajectory}");

    // Load the Deep Space Network ground stations.
    // Nyx allows you to build these at runtime but it's pretty static so we can just load them from YAML.
    let ground_station_file: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "examples",
        "04_lro_od",
        "dsn-network.yaml",
    ]
    .iter()
    .collect();

    let devices = GroundStation::load_many(ground_station_file).unwrap();

    // Typical OD software requires that you specify your own tracking schedule or you'll have overlapping measurements.
    // Nyx can build a tracking schedule for you based on the first station with access.
    let trkconfg_yaml: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "examples",
        "04_lro_od",
        "tracking-cfg.yaml",
    ]
    .iter()
    .collect();

    let configs: BTreeMap<String, TrkConfig> = TrkConfig::load_named(trkconfg_yaml).unwrap();

    dbg!(&configs);

    // Build the tracking arc simulation to generate a "standard measurement".
    let mut trk = TrackingArcSim::<Spacecraft, RangeDoppler, _>::with_seed(
        devices, trajectory, configs, 12345,
    )
    .unwrap();

    trk.build_schedule(almanac.clone()).unwrap();
    let arc = trk.generate_measurements(almanac).unwrap();
    // Save the simulated tracking data
    arc.to_parquet_simple("./04_lro_simulated_tracking.parquet")?;

    println!("{arc}");

    Ok(())
}
