// #![doc = include_str!("./README.md")]
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::{
    almanac::metaload::MetaFile,
    constants::{
        celestial_objects::{EARTH, JUPITER_BARYCENTER, MOON, SUN},
        frames::{EARTH_J2000, MOON_J2000, MOON_PA_FRAME},
    },
    prelude::Almanac,
};
use hifitime::{Epoch, TimeSeries, TimeUnits, Unit};
use nyx::{
    cosmic::{Aberration, Frame, Mass, MetaAlmanac, SRPData},
    dynamics::{
        guidance::LocalFrame, GravityField, OrbitalDynamics, SolarPressure, SpacecraftDynamics,
    },
    io::{ConfigRepr, ExportCfg},
    md::prelude::{GravityFieldData, Traj},
    od::{
        msr::MeasurementType,
        prelude::{KalmanVariant, TrackingArcSim, TrkConfig},
        process::{Estimate, NavSolution, ResidRejectCrit, SpacecraftUncertainty},
        snc::ProcessNoise3D,
        GroundStation, SpacecraftKalmanOD, SpacecraftKalmanScalarOD,
    },
    propagators::{IntegratorOptions, Propagator},
    Orbit, Spacecraft, State,
};

use std::{collections::BTreeMap, error::Error, path::PathBuf, str::FromStr, sync::Arc};

// TODO: Convert this to a Spacecraft Sequence

fn main() -> Result<(), Box<dyn Error>> {
    pel::init();

    // ====================== //
    // === ALMANAC SET UP === //
    // ====================== //

    // Dynamics models require planetary constants and ephemerides to be defined.
    // Let's start by grabbing those by using ANISE's MetaAlmanac.

    let data_folder: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "examples",
        "06_lunar_orbit_determination",
    ]
    .iter()
    .collect();

    let meta = data_folder.join("metaalmanac.dhall");

    // Load this ephem in the general Almanac we're using for this analysis.
    let almanac = MetaAlmanac::new(meta.to_string_lossy().as_ref())
        .map_err(Box::new)?
        .process(true)
        .map_err(Box::new)?;

    // Lock the almanac (an Arc is a read only structure).
    let almanac = Arc::new(almanac);

    // Build a nominal trajectory
    // TODO: Switch this to a sequence once the OD over a spacecraft sequence is implemented.

    let epoch = Epoch::from_gregorian_utc_at_noon(2024, 2, 29);
    let moon_j2000 = almanac.frame_info(MOON_J2000)?;

    // To build the trajectory we need to provide a spacecraft template.
    let orbiter = Spacecraft::builder()
        .mass(Mass::from_dry_and_prop_masses(1018.0, 900.0))
        .srp(SRPData {
            area_m2: 3.9 * 2.7,
            coeff_reflectivity: 0.96,
        })
        .orbit(Orbit::try_keplerian_altitude(
            150.0, 0.00212, 33.6, 45.0, 45.0, 0.0, epoch, moon_j2000,
        )?) // Setting a zero orbit here because it's just a template
        .build();

    // ========================== //
    // === BUILD NOMINAL TRAJ === //
    // ========================== //

    // Set up the spacecraft dynamics.

    // Specify that the orbital dynamics must account for the graviational pull of the Earth and the Sun.
    // The gravity of the Moon will also be accounted for since the spaceraft in a lunar orbit.
    let mut orbital_dyn = OrbitalDynamics::point_masses(vec![EARTH, SUN, JUPITER_BARYCENTER]);

    // We want to include the spherical harmonics, so let's download the gravitational data from the Nyx Cloud.
    // We're using the GRAIL JGGRX model.
    let mut jggrx_meta = MetaFile {
        uri: "http://public-data.nyxspace.com/nyx/models/Luna_jggrx_1500e_sha.tab.gz".to_string(),
        crc32: Some(0x6bcacda8), // Specifying the CRC32 avoids redownloading it if it's cached.
    };
    // And let's download it if we don't have it yet.
    jggrx_meta.process(true)?;

    // Build the spherical harmonics.
    // The harmonics must be computed in the body fixed frame.
    // We're using the long term prediction of the Moon principal axes frame.
    let moon_pa_frame = MOON_PA_FRAME.with_orient(31008);
    let sph_harmonics = GravityField::from_stor(
        almanac.frame_info(moon_pa_frame)?,
        GravityFieldData::from_shadr(&jggrx_meta.uri, 80, 80, true)?,
    );

    // Include the spherical harmonics into the orbital dynamics.
    orbital_dyn.accel_models.push(sph_harmonics);

    // We define the solar radiation pressure, using the default solar flux and accounting only
    // for the eclipsing caused by the Earth and Moon.
    // Note that by default, enabling the SolarPressure model will also enable the estimation of the coefficient of reflectivity.
    let srp_dyn = SolarPressure::new(vec![MOON_J2000], almanac.clone())?;

    // Finalize setting up the dynamics, specifying the force models (orbital_dyn) separately from the
    // acceleration models (SRP in this case). Use `from_models` to specify multiple accel models.
    let dynamics = SpacecraftDynamics::from_model(orbital_dyn, srp_dyn);

    println!("{dynamics}");

    let setup = Propagator::rk89(dynamics.clone(), IntegratorOptions::default());

    let truth_traj = setup
        .with(orbiter, almanac.clone())
        .for_duration_with_traj(Unit::Day * 2)?
        .1;

    // ==================== //
    // === OD SIMULATOR === //
    // ==================== //

    // Load the Deep Space Network ground stations.
    // Nyx allows you to build these at runtime but it's pretty static so we can just load them from YAML.
    let ground_station_file = data_folder.join("dsn-network.yaml");
    let devices = GroundStation::load_named(ground_station_file)?;

    let proc_devices = devices.clone();

    // Typical OD software requires that you specify your own tracking schedule or you'll have overlapping measurements.
    // Nyx can build a tracking schedule for you based on the first station with access.
    let configs: BTreeMap<String, TrkConfig> =
        TrkConfig::load_named(data_folder.join("tracking-cfg.yaml"))?;

    // Build the tracking arc simulation to generate a "standard measurement".
    let mut trk = TrackingArcSim::<Spacecraft, GroundStation>::with_seed(
        devices.clone(),
        truth_traj.clone(),
        configs,
        123, // Set a seed for reproducibility
    )?;

    trk.build_schedule(almanac.clone())?;
    let arc = trk.generate_measurements(almanac.clone())?;
    // Save the simulated tracking data
    arc.to_parquet_simple("./data/04_output/06_lunar_simulated_tracking.parquet")?;

    // We'll note that in our case, we have continuous coverage of LRO when the vehicle is not behind the Moon.
    println!("{arc}");

    // Now that we have simulated measurements, we'll run the orbit determination.

    // ===================== //
    // === OD ESTIMATION === //
    // ===================== //

    let sc = SpacecraftUncertainty::builder()
        .nominal(orbiter)
        .frame(LocalFrame::RIC)
        .x_km(0.5)
        .y_km(0.5)
        .z_km(0.5)
        .vx_km_s(5e-3)
        .vy_km_s(5e-3)
        .vz_km_s(5e-3)
        .build();

    // Build the filter initial estimate, which we will reuse in the filter.
    let initial_estimate = sc.to_estimate()?;

    println!("== FILTER STATE ==\n{orbiter:x}\n{initial_estimate}");

    // Build the SNC in the Moon J2000 frame, specified as a velocity noise over time.
    let process_noise = ProcessNoise3D::from_velocity_km_s(
        &[1e-12, 1e-12, 1e-12],
        1 * Unit::Hour,
        10 * Unit::Minute,
        None,
    );

    println!("{process_noise}");

    // We'll set up the OD process to reject measurements whose residuals are move than 3 sigmas away from what we expect.
    let odp = SpacecraftKalmanScalarOD::new(
        setup,
        KalmanVariant::ReferenceUpdate,
        Some(ResidRejectCrit::default()),
        proc_devices,
        almanac.clone(),
    )
    .with_process_noise(process_noise);

    let od_sol = odp.process_arc(initial_estimate, &arc)?;

    let final_est = od_sol.estimates.last().unwrap();

    println!("{final_est}");

    let ric_err = truth_traj
        .at(final_est.epoch())?
        .orbit
        .ric_difference(&final_est.orbital_state())?;
    println!("== RIC at end ==");
    println!("RIC Position (m): {:.3}", ric_err.radius_km * 1e3);
    println!("RIC Velocity (m/s): {:.3}", ric_err.velocity_km_s * 1e3);

    println!(
        "Num residuals rejected: #{}",
        od_sol.rejected_residuals().len()
    );
    println!(
        "Percentage within +/-3: {}",
        od_sol.residual_ratio_within_threshold(3.0).unwrap()
    );
    println!("Ratios normal? {}", od_sol.is_normal(None).unwrap());

    od_sol.to_parquet(
        "./data/04_output/06_lunar_od_results.parquet",
        ExportCfg::default(),
    )?;

    Ok(())
}
