extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger as pel;

use anise::{
    almanac::Almanac,
    constants::{
        celestial_objects::{EARTH, SUN},
        frames::{EARTH_J2000, IAU_MOON_FRAME, MOON_J2000},
    },
};
use hifitime::{Epoch, TimeUnits};
use nyx::{
    cosmic::Aberration,
    dynamics::{OrbitalDynamics, SpacecraftDynamics},
    io::ExportCfg,
    md::prelude::{IntegratorOptions, Propagator},
    od::interlink::InterlinkTxSpacecraft,
    od::noise::link_specific,
    od::prelude::*,
    Orbit, Spacecraft, State,
};

use indexmap::{IndexMap, IndexSet};
use std::{collections::BTreeMap, error::Error, path::PathBuf, sync::Arc};

// Test the Cislunar Autonomous Positioning System (CAPS) similar to how it's flown on CAPSTONE.
/// Assume that the NRHO orbiter is the transmitter in a two-way communication with a low lunar orbiter.
/// This is a Spacecraft to Spacecraft Orbit Determination Process (S2SODP).
///
/// We test that with dispersions we can still converge on a better than the original dispersion.
/// NOTE: In this short tracking arc, we do not converge well because we don't have good enough visibility
/// of the crosstrack. This is reflected in the covariance.
fn main() -> Result<(), Box<dyn Error>> {
    pel::init();

    // ====================== //
    // === ALMANAC SET UP === //
    // ====================== //

    let manifest_dir =
        PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap_or(".".to_string()));

    let out = manifest_dir.join("data/04_output/");

    let almanac = Arc::new(
        Almanac::new(
            &manifest_dir
                .join("data/01_planetary/pck08.pca")
                .to_string_lossy(),
        )
        .unwrap()
        .load(
            &manifest_dir
                .join("data/01_planetary/de440s.bsp")
                .to_string_lossy(),
        )
        .unwrap(),
    );

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
    let moon_iau = almanac.frame_info(IAU_MOON_FRAME).unwrap();

    let epoch = Epoch::from_gregorian_tai(2021, 5, 29, 19, 51, 16, 852_000);
    let nrho = Orbit::cartesian(
        166_473.631_302_239_7,
        -274_715.487_253_382_7,
        -211_233.210_176_686_7,
        0.933_451_604_520_018_4,
        0.436_775_046_841_900_9,
        -0.082_211_021_250_348_95,
        epoch,
        eme2k,
    );

    let tx_nrho_sc = Spacecraft::from(nrho);

    let state_luna = almanac.transform_to(nrho, MOON_J2000, None).unwrap();
    println!("Start state (dynamics: Earth, Moon, Sun gravity):\n{state_luna}");

    let bodies = vec![EARTH, SUN];
    let dynamics = SpacecraftDynamics::new(OrbitalDynamics::point_masses(bodies));

    let setup = Propagator::rk89(
        dynamics,
        IntegratorOptions::builder().max_step(0.5.minutes()).build(),
    );

    /* == Propagate the NRHO vehicle == */
    let prop_time = 1.1 * state_luna.period().unwrap();

    let (nrho_final, mut tx_traj) = setup
        .with(tx_nrho_sc, almanac.clone())
        .for_duration_with_traj(prop_time)
        .unwrap();

    tx_traj.name = Some("NRHO Tx SC".to_string());

    println!("{tx_traj}");

    /* == Propagate an LLO vehicle == */
    let llo_orbit =
        Orbit::try_keplerian_altitude(110.0, 1e-4, 90.0, 0.0, 0.0, 0.0, epoch, moon_iau).unwrap();

    let llo_sc = Spacecraft::builder().orbit(llo_orbit).build();

    let (_, llo_traj) = setup
        .with(llo_sc, almanac.clone())
        .until_epoch_with_traj(nrho_final.epoch())
        .unwrap();

    // Export the subset of the first two hours.
    llo_traj
        .clone()
        .filter_by_offset(..2.hours())
        .to_parquet_simple(out.join("05_caps_llo_truth.pq"), almanac.clone())?;

    /* == Setup the interlink == */

    let mut measurement_types = IndexSet::new();
    measurement_types.insert(MeasurementType::Range);
    measurement_types.insert(MeasurementType::Doppler);

    let mut stochastics = IndexMap::new();

    let sa45_csac_allan_dev = 1e-11;

    stochastics.insert(
        MeasurementType::Range,
        StochasticNoise::from_hardware_range_km(
            sa45_csac_allan_dev,
            10.0.seconds(),
            link_specific::ChipRate::StandardT4B,
            link_specific::SN0::Average,
        ),
    );

    stochastics.insert(
        MeasurementType::Doppler,
        StochasticNoise::from_hardware_doppler_km_s(
            sa45_csac_allan_dev,
            10.0.seconds(),
            link_specific::CarrierFreq::SBand,
            link_specific::CN0::Average,
        ),
    );

    let interlink = InterlinkTxSpacecraft {
        traj: tx_traj,
        measurement_types,
        integration_time: None,
        timestamp_noise_s: None,
        ab_corr: Aberration::LT,
        stochastic_noises: Some(stochastics),
    };

    // Devices are the transmitter, which is our NRHO vehicle.
    let mut devices = BTreeMap::new();
    devices.insert("NRHO Tx SC".to_string(), interlink);

    let mut configs = BTreeMap::new();
    configs.insert(
        "NRHO Tx SC".to_string(),
        TrkConfig::builder()
            .strands(vec![Strand {
                start: epoch,
                end: nrho_final.epoch(),
            }])
            .build(),
    );

    let mut trk_sim =
        TrackingArcSim::with_seed(devices.clone(), llo_traj.clone(), configs, 0).unwrap();
    println!("{trk_sim}");

    let trk_data = trk_sim.generate_measurements(almanac.clone()).unwrap();
    println!("{trk_data}");

    trk_data
        .to_parquet_simple(out.clone().join("nrho_interlink_msr.pq"))
        .unwrap();

    // Run a truth OD where we estimate the LLO position
    let llo_uncertainty = SpacecraftUncertainty::builder()
        .nominal(llo_sc)
        .x_km(1.0)
        .y_km(1.0)
        .z_km(1.0)
        .vx_km_s(1e-3)
        .vy_km_s(1e-3)
        .vz_km_s(1e-3)
        .build();

    let mut proc_devices = devices.clone();

    // Define the initial estimate, randomized, seed for reproducibility
    let mut initial_estimate = llo_uncertainty.to_estimate_randomized(Some(0)).unwrap();
    // Inflate the covariance -- https://github.com/nyx-space/nyx/issues/339
    initial_estimate.covar *= 2.5;

    // Increase the noise in the devices to accept more measurements.

    for link in proc_devices.values_mut() {
        for noise in &mut link.stochastic_noises.as_mut().unwrap().values_mut() {
            *noise.white_noise.as_mut().unwrap() *= 3.0;
        }
    }

    let init_err = initial_estimate
        .orbital_state()
        .ric_difference(&llo_orbit)
        .unwrap();

    println!("initial estimate:\n{initial_estimate}");
    println!("RIC errors = {init_err}",);

    let odp = InterlinkKalmanOD::new(
        setup.clone(),
        KalmanVariant::ReferenceUpdate,
        Some(ResidRejectCrit::default()),
        proc_devices,
        almanac.clone(),
    );

    // Shrink the data to process.
    let arc = trk_data.filter_by_offset(..2.hours());

    let od_sol = odp.process_arc(initial_estimate, &arc).unwrap();

    println!("{od_sol}");

    od_sol
        .to_parquet(
            out.join("05_caps_interlink_od_sol.pq"),
            ExportCfg::default(),
        )
        .unwrap();

    let od_traj = od_sol.to_traj().unwrap();

    od_traj
        .ric_diff_to_parquet(
            &llo_traj,
            out.join("05_caps_interlink_llo_est_error.pq"),
            ExportCfg::default(),
        )
        .unwrap();

    let final_est = od_sol.estimates.last().unwrap();
    assert!(final_est.within_3sigma(), "should be within 3 sigma");

    println!("ESTIMATE\n{final_est:x}\n");
    let truth = llo_traj.at(final_est.epoch()).unwrap();
    println!("TRUTH\n{truth:x}");

    let final_err = truth
        .orbit
        .ric_difference(&final_est.orbital_state())
        .unwrap();
    println!("ERROR {final_err}");

    // Build the residuals versus reference plot.
    let rvr_sol = odp
        .process_arc(initial_estimate, &arc.resid_vs_ref_check())
        .unwrap();

    rvr_sol
        .to_parquet(
            out.join("05_caps_interlink_resid_v_ref.pq"),
            ExportCfg::default(),
        )
        .unwrap();

    let final_rvr = rvr_sol.estimates.last().unwrap();

    println!("RMAG error {:.3} m", final_err.rmag_km() * 1e3);
    println!(
        "Pure prop error {:.3} m",
        final_rvr
            .orbital_state()
            .ric_difference(&final_est.orbital_state())
            .unwrap()
            .rmag_km()
            * 1e3
    );

    Ok(())
}
