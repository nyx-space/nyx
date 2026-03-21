use nyx_space::dynamics::orbital::OrbitalDynamics;
use nyx_space::dynamics::SpacecraftDynamics;
use nyx_space::md::prelude::Propagator;
use nyx_space::od::msr::MeasurementType;
use nyx_space::od::noise::{StochasticNoise, WhiteNoise};
use nyx_space::od::process::KalmanODProcess;
use nyx_space::od::simulator::TrackingArcSim;
use nyx_space::od::xyz::XyzDevice;
use nyx_space::Spacecraft;
use anise::prelude::Almanac;
use nalgebra::{Const, OMatrix};
use std::sync::Arc;
use nyx_space::od::kalman::KalmanVariant;
use std::collections::BTreeMap;
use nyx_space::od::prelude::KfEstimate;
use hifitime::Unit;
use nyx_space::od::prelude::TrkConfig;
use anise::constants::frames::EARTH_J2000;
use nyx_space::cosmic::Orbit;
use nyx_space::time::Epoch;
use nyx_space::State;

#[test]
fn test_xyz_filtering() {
    let almanac = crate::test_almanac_arcd();

    let eme2k = almanac.frame_info(EARTH_J2000).unwrap();
    let dt = Epoch::from_gregorian_utc_hms(2020, 1, 1, 4, 0, 0);
    let sim_sc = Spacecraft::from(Orbit::keplerian(
        22000.0, 0.01, 30.0, 80.0, 40.0, 170.0, dt, eme2k,
    ));

    let dyns = SpacecraftDynamics::new(OrbitalDynamics::point_masses(vec![anise::constants::celestial_objects::EARTH]));
    let sim_prop = Propagator::default(dyns.clone());
    let (_, sim_traj) = sim_prop.with(sim_sc, almanac.clone()).for_duration_with_traj(6.0 * Unit::Hour).unwrap();

    let noise_x = StochasticNoise {
        white_noise: Some(WhiteNoise::new(1e-3, 1.0 * Unit::Second).unwrap()),
        bias: None,
    };
    let noise_y = StochasticNoise {
        white_noise: Some(WhiteNoise::new(1e-3, 1.0 * Unit::Second).unwrap()),
        bias: None,
    };
    let noise_z = StochasticNoise {
        white_noise: Some(WhiteNoise::new(1e-3, 1.0 * Unit::Second).unwrap()),
        bias: None,
    };

    let device = XyzDevice::new("GPS".to_string())
        .with_noise(MeasurementType::X, noise_x)
        .with_noise(MeasurementType::Y, noise_y)
        .with_noise(MeasurementType::Z, noise_z);

    let mut devices = BTreeMap::new();
    devices.insert(device.name.clone(), device.clone());

    let configs = BTreeMap::from([(
        device.name.clone(),
        TrkConfig::from_sample_rate(1.0 * Unit::Minute),
    )]);

    let mut arc_sim = TrackingArcSim::with_seed(
        devices.clone(),
        sim_traj.clone(),
        configs,
        12345
    ).unwrap();

    let msr_arc = arc_sim.generate_measurements(almanac.clone()).unwrap();

    let mut est_sc = sim_sc;
    est_sc.orbit.radius_km.x += 1.0;
    est_sc.orbit.radius_km.y -= 1.0;
    est_sc.orbit.radius_km.z += 1.0;
    est_sc.orbit.velocity_km_s.x += 1e-3;
    est_sc.orbit.velocity_km_s.y -= 1e-3;
    est_sc.orbit.velocity_km_s.z += 1e-3;

    let od = KalmanODProcess::<SpacecraftDynamics, Const<3>, Const<3>, XyzDevice>::new(
        sim_prop,
        KalmanVariant::DeviationTracking,
        None,
        devices,
        almanac.clone(),
    );

    let mut initial_cov = OMatrix::<f64, Const<9>, Const<9>>::zeros();
    initial_cov[(0, 0)] = 1.0;
    initial_cov[(1, 1)] = 1.0;
    initial_cov[(2, 2)] = 1.0;
    initial_cov[(3, 3)] = 1e-6;
    initial_cov[(4, 4)] = 1e-6;
    initial_cov[(5, 5)] = 1e-6;

    let initial_estimate = KfEstimate::from_covar(est_sc, initial_cov);

    let od_res = od.process_arc(initial_estimate, &msr_arc).unwrap();

    let final_state = od_res.estimates.last().unwrap().nominal_state;
    let true_state = sim_traj.at(final_state.epoch()).unwrap();

    let diff_x = (final_state.orbit.radius_km.x - true_state.orbit.radius_km.x).abs();
    let diff_y = (final_state.orbit.radius_km.y - true_state.orbit.radius_km.y).abs();
    let diff_z = (final_state.orbit.radius_km.z - true_state.orbit.radius_km.z).abs();

    println!("Diff X: {}, Y: {}, Z: {}", diff_x, diff_y, diff_z);

    assert!(diff_x < 0.1);
    assert!(diff_y < 0.1);
    assert!(diff_z < 0.1);
}
