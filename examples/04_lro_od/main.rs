#![doc = include_str!("./README.md")]
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::{
    almanac::metaload::MetaFile,
    constants::{
        celestial_objects::{EARTH, JUPITER_BARYCENTER, MOON, SUN},
        frames::{EARTH_J2000, MOON_J2000, MOON_PA_FRAME},
    },
};
use hifitime::{Epoch, TimeUnits, Unit};
use nyx::{
    cosmic::{Aberration, Frame, Mass, MetaAlmanac, SRPData},
    dynamics::{
        guidance::LocalFrame, Harmonics, OrbitalDynamics, SolarPressure, SpacecraftDynamics,
    },
    io::{ConfigRepr, ExportCfg},
    md::prelude::{HarmonicsMem, Traj},
    od::{
        msr::MeasurementType,
        prelude::{KalmanVariant, TrackingArcSim, TrkConfig},
        process::{Estimate, NavSolution, ResidRejectCrit, SpacecraftUncertainty},
        snc::ProcessNoise3D,
        GroundStation, SpacecraftKalmanOD,
    },
    propagators::Propagator,
    Orbit, Spacecraft, State,
};

use std::{collections::BTreeMap, error::Error, path::PathBuf, str::FromStr, sync::Arc};

fn main() -> Result<(), Box<dyn Error>> {
    pel::init();

    // ====================== //
    // === ALMANAC SET UP === //
    // ====================== //

    // Dynamics models require planetary constants and ephemerides to be defined.
    // Let's start by grabbing those by using ANISE's MetaAlmanac.

    let data_folder: PathBuf = [env!("CARGO_MANIFEST_DIR"), "examples", "04_lro_od"]
        .iter()
        .collect();

    let meta = data_folder.join("lro-dynamics.dhall");

    // Load this ephem in the general Almanac we're using for this analysis.
    let mut almanac = MetaAlmanac::new(meta.to_string_lossy().as_ref())
        .map_err(Box::new)?
        .process(true)
        .map_err(Box::new)?;

    let mut moon_pc = almanac.planetary_data.get_by_id(MOON)?;
    moon_pc.mu_km3_s2 = 4902.74987;
    almanac.planetary_data.set_by_id(MOON, moon_pc)?;

    let mut earth_pc = almanac.planetary_data.get_by_id(EARTH)?;
    earth_pc.mu_km3_s2 = 398600.436;
    almanac.planetary_data.set_by_id(EARTH, earth_pc)?;

    // Save this new kernel for reuse.
    // In an operational context, this would be part of the "Lock" process, and should not change throughout the mission.
    almanac
        .planetary_data
        .save_as(&data_folder.join("lro-specific.pca"), true)?;

    // Lock the almanac (an Arc is a read only structure).
    let almanac = Arc::new(almanac);

    // Orbit determination requires a Trajectory structure, which can be saved as parquet file.
    // In our case, the trajectory comes from the BSP file, so we need to build a Trajectory from the almanac directly.
    // To query the Almanac, we need to build the LRO frame in the J2000 orientation in our case.
    // Inspecting the LRO BSP in the ANISE GUI shows us that NASA has assigned ID -85 to LRO.
    let lro_frame = Frame::from_ephem_j2000(-85);

    // To build the trajectory we need to provide a spacecraft template.
    let sc_template = Spacecraft::builder()
        .mass(Mass::from_dry_and_prop_masses(1018.0, 900.0)) // Launch masses
        .srp(SRPData {
            // SRP configuration is arbitrary, but we will be estimating it anyway.
            area_m2: 3.9 * 2.7,
            coeff_reflectivity: 0.96,
        })
        .orbit(Orbit::zero(MOON_J2000)) // Setting a zero orbit here because it's just a template
        .build();
    // Now we can build the trajectory from the BSP file.
    // We'll arbitrarily set the tracking arc to 24 hours with a five second time step.
    let traj_as_flown = Traj::from_bsp(
        lro_frame,
        MOON_J2000,
        almanac.clone(),
        sc_template,
        5.seconds(),
        Some(Epoch::from_str("2024-01-01 00:00:00 UTC")?),
        Some(Epoch::from_str("2024-01-02 00:00:00 UTC")?),
        Aberration::LT,
        Some("LRO".to_string()),
    )?;

    println!("{traj_as_flown}");

    // ====================== //
    // === MODEL MATCHING === //
    // ====================== //

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
    let sph_harmonics = Harmonics::from_stor(
        almanac.frame_info(moon_pa_frame)?,
        HarmonicsMem::from_shadr(&jggrx_meta.uri, 80, 80, true)?,
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
    let setup = Propagator::default_dp78(dynamics.clone());

    // For reference, let's build the trajectory with Nyx's models from that LRO state.
    let (sim_final, traj_as_sim) = setup
        .with(*traj_as_flown.first(), almanac.clone())
        .until_epoch_with_traj(traj_as_flown.last().epoch())?;

    println!("SIM INIT:  {:x}", traj_as_flown.first());
    println!("SIM FINAL: {sim_final:x}");
    // Compute RIC difference between SIM and LRO ephem
    let sim_lro_delta = sim_final
        .orbit
        .ric_difference(&traj_as_flown.last().orbit)?;
    println!("{traj_as_sim}");
    println!(
        "SIM v LRO - RIC Position (m): {:.3}",
        sim_lro_delta.radius_km * 1e3
    );
    println!(
        "SIM v LRO - RIC Velocity (m/s): {:.3}",
        sim_lro_delta.velocity_km_s * 1e3
    );

    traj_as_sim.ric_diff_to_parquet(
        &traj_as_flown,
        "./04_lro_sim_truth_error.parquet",
        ExportCfg::default(),
    )?;

    // ==================== //
    // === OD SIMULATOR === //
    // ==================== //

    // After quite some time trying to exactly match the model, we still end up with an oscillatory difference on the order of 150 meters between the propagated state
    // and the truth LRO state.

    // Therefore, we will actually run an estimation from a dispersed LRO state.
    // The sc_seed is the true LRO state from the BSP.
    let sc_seed = *traj_as_flown.first();

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

    let devices = GroundStation::load_named(ground_station_file)?;

    let mut proc_devices = devices.clone();

    // Increase the noise in the devices to accept more measurements.
    for gs in proc_devices.values_mut() {
        if let Some(noise) = &mut gs
            .stochastic_noises
            .as_mut()
            .unwrap()
            .get_mut(&MeasurementType::Range)
        {
            *noise.white_noise.as_mut().unwrap() *= 3.0;
        }
    }

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

    // Build the tracking arc simulation to generate a "standard measurement".
    let mut trk = TrackingArcSim::<Spacecraft, GroundStation>::with_seed(
        devices.clone(),
        traj_as_flown.clone(),
        configs,
        123, // Set a seed for reproducibility
    )?;

    trk.build_schedule(almanac.clone())?;
    let arc = trk.generate_measurements(almanac.clone())?;
    // Save the simulated tracking data
    arc.to_parquet_simple("./04_lro_simulated_tracking.parquet")?;

    // We'll note that in our case, we have continuous coverage of LRO when the vehicle is not behind the Moon.
    println!("{arc}");

    // Now that we have simulated measurements, we'll run the orbit determination.

    // ===================== //
    // === OD ESTIMATION === //
    // ===================== //

    let sc = SpacecraftUncertainty::builder()
        .nominal(sc_seed)
        .frame(LocalFrame::RIC)
        .x_km(0.5)
        .y_km(0.5)
        .z_km(0.5)
        .vx_km_s(5e-3)
        .vy_km_s(5e-3)
        .vz_km_s(5e-3)
        .build();

    // Build the filter initial estimate, which we will reuse in the filter.
    let mut initial_estimate = sc.to_estimate()?;
    initial_estimate.covar *= 3.0;

    println!("== FILTER STATE ==\n{sc_seed:x}\n{initial_estimate}");

    // Build the SNC in the Moon J2000 frame, specified as a velocity noise over time.
    let process_noise = ProcessNoise3D::from_velocity_km_s(
        &[1e-10, 1e-10, 1e-10],
        1 * Unit::Hour,
        10 * Unit::Minute,
        None,
    );

    println!("{process_noise}");

    // We'll set up the OD process to reject measurements whose residuals are move than 3 sigmas away from what we expect.
    let odp = SpacecraftKalmanOD::new(
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

    let ric_err = traj_as_flown
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

    od_sol.to_parquet("./04_lro_od_results.parquet", ExportCfg::default())?;

    // In our case, we have the truth trajectory from NASA.
    // So we can compute the RIC state difference between the real LRO ephem and what we've just estimated.
    // Export the OD trajectory first.
    let od_trajectory = od_sol.to_traj()?;
    // Build the RIC difference.
    od_trajectory.ric_diff_to_parquet(
        &traj_as_flown,
        "./04_lro_od_truth_error.parquet",
        ExportCfg::default(),
    )?;

    Ok(())
}
