#![doc = include_str!("./README.md")]
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::{
    almanac::metaload::MetaFile,
    constants::{
        celestial_objects::{EARTH, SUN},
        frames::{EARTH_J2000, MOON_J2000, MOON_PA_FRAME},
    },
};
use hifitime::{Epoch, TimeUnits, Unit};
use nyx::{
    cosmic::{Aberration, Frame, MetaAlmanac, SrpConfig},
    dynamics::{
        guidance::LocalFrame, Harmonics, OrbitalDynamics, SolarPressure, SpacecraftDynamics,
    },
    io::{ConfigRepr, ExportCfg},
    md::prelude::{HarmonicsMem, Traj},
    od::{
        msr::RangeDoppler,
        prelude::{TrackingArcSim, TrkConfig, KF},
        process::{EkfTrigger, IterationConf, ODProcess, SpacecraftUncertainty},
        snc::SNC3,
        GroundStation,
    },
    propagators::Propagator,
    Orbit, Spacecraft,
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
        Some(Epoch::from_str("2024-01-01 00:00:00 UTC")?),
        Some(Epoch::from_str("2024-01-03 00:00:00 UTC")?),
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

    let devices = GroundStation::load_many(ground_station_file)?;

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

    let configs: BTreeMap<String, TrkConfig> = TrkConfig::load_named(trkconfg_yaml)?;

    dbg!(&configs);

    // Build the tracking arc simulation to generate a "standard measurement".
    let mut trk = TrackingArcSim::<Spacecraft, RangeDoppler, _>::with_seed(
        devices,
        trajectory.clone(),
        configs,
        12345,
    )?;

    trk.build_schedule(almanac.clone())?;
    let arc = trk.generate_measurements(almanac.clone())?;
    // Save the simulated tracking data
    arc.to_parquet_simple("./04_lro_simulated_tracking.parquet")?;

    // We'll note that in our case, we have continuous coverage of LRO when the vehicle is not behind the Moon.
    println!("{arc}");

    // Now that we have simulated measurement, we'll run the orbit determination.

    // First, we need to set up the propagator used in the OD.
    // Set up the spacecraft dynamics.

    // Specify that the orbital dynamics must account for the graviational pull of the Earth and the Sun.
    // The gravity of the Moon will also be accounted for since the spaceraft in a lunar orbit.
    let mut orbital_dyn = OrbitalDynamics::point_masses(vec![EARTH, SUN]);

    // We want to include the spherical harmonics, so let's download the gravitational data from the Nyx Cloud.
    // We're using the GRAIL JGGRX model.
    let mut jggrx_meta = MetaFile {
        uri: "http://public-data.nyxspace.com/nyx/models/Luna_jggrx_1500e_sha.tab.gz".to_string(),
        crc32: Some(0x6bcacda8), // Specifying the CRC32 avoids redownloading it if it's cached.
    };
    // And let's download it if we don't have it yet.
    jggrx_meta.process()?;

    // Build the spherical harmonics.
    // The harmonics must be computed in the body fixed frame.
    // We're using the long term prediction of the Earth centered Earth fixed frame, IAU Earth.
    let sph_harmonics = Harmonics::from_stor(
        almanac.frame_from_uid(MOON_PA_FRAME.with_orient(31008))?,
        HarmonicsMem::from_shadr(&jggrx_meta.uri, 200, 200, true)?,
    );

    // Include the spherical harmonics into the orbital dynamics.
    orbital_dyn.accel_models.push(sph_harmonics);

    // We define the solar radiation pressure, using the default solar flux and accounting only
    // for the eclipsing caused by the Earth and Moon.
    // Note that by default, enabling the SolarPressure model will also enable the estimation of the coefficient of reflectivity.
    let srp_dyn = SolarPressure::new(vec![EARTH_J2000, MOON_J2000], almanac.clone())?;

    // Finalize setting up the dynamics, specifying the force models (orbital_dyn) separately from the
    // acceleration models (SRP in this case). Use `from_models` to specify multiple accel models.
    let dynamics = SpacecraftDynamics::from_model(orbital_dyn, srp_dyn);

    println!("{dynamics}");

    // Now we can build the propagator.
    let setup = Propagator::default_dp78(dynamics);

    // For an OD arc, we need to start with an initial estimate and covariance.
    // The ephem published by NASA does not include the covariance. Instead, we'll make one up!

    let sc_seed = trajectory.first().with_stm();

    let sc = SpacecraftUncertainty::builder()
        .nominal(sc_seed)
        .frame(LocalFrame::RIC)
        .x_km(1.0)
        .y_km(1.0)
        .z_km(1.0)
        .vx_km_s(0.5e-1)
        .vy_km_s(0.5e-1)
        .vz_km_s(0.5e-1)
        .cr(0.2)
        .build();

    println!("{sc}");

    let initial_estimate = sc.to_estimate()?;

    let kf = KF::new(
        initial_estimate,
        SNC3::from_diagonal(2 * Unit::Minute, &[1e-14, 1e-14, 1e-14]),
    );

    let mut odp = ODProcess::ekf(
        setup.with(sc_seed, almanac.clone()),
        kf,
        EkfTrigger::new(30_000, Unit::Minute * 2),
        None,
        almanac,
    );

    odp.process_arc::<GroundStation>(&arc)?;
    // Let's run a smoother just to see that the filter won't run it if the RSS error is small.
    odp.iterate_arc::<GroundStation>(&arc, IterationConf::once())?;

    odp.to_parquet("./04_lro_od_results.parquet", ExportCfg::default())?;

    // In our case, we have the truth trajectory from NASA.
    // So we can compute the RIC state difference between the real LRO ephem and what we've just estimated.
    // Export the OD trajectory first.
    let od_trajectory = odp.to_traj()?;
    // Build the RIC difference.
    od_trajectory.ric_diff_to_parquet(
        &trajectory,
        "./04_lro_od_truth_error.parquet",
        ExportCfg::default(),
    )?;

    Ok(())
}
