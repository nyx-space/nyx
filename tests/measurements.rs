extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn nil_measurement() {
    use hifitime::julian::*;
    use nyx::celestia::{State, ECEF};
    use nyx::od::Measurement;
    use nyx::od::ranging::GroundStation;
    use std::f64::EPSILON;
    // Let's create a station and make it estimate the range and range rate of something which is strictly in the same spot.

    let lat = -7.906_635_7;
    let long = 345.5975;
    let height = 56.0e-3;
    let dt = ModifiedJulian::j2000();

    let station = GroundStation::from_noise_values("nil".to_string(), 0.0, lat, long, height, 0.0, 0.0);

    let at_station = State::<ECEF>::from_geodesic(lat, long, height, dt);

    let meas = station.measure(at_station, dt.into_instant());

    let h_tilde = *meas.sensitivity();
    assert!(h_tilde[(0, 0)].is_nan(), "expected NaN");
    assert!(h_tilde[(1, 0)].is_nan(), "expected NaN");
    assert!(h_tilde[(2, 0)].is_nan(), "expected NaN");
    assert!(h_tilde[(0, 1)].is_nan(), "expected NaN");
    assert!(h_tilde[(1, 1)].is_nan(), "expected NaN");
    assert!(h_tilde[(2, 1)].is_nan(), "expected NaN");
    assert!(h_tilde[(3, 1)].is_nan(), "expected NaN");
    assert!(h_tilde[(4, 1)].is_nan(), "expected NaN");
    assert!(h_tilde[(5, 1)].is_nan(), "expected NaN");

    assert!(meas.observation()[(0, 0)] - 0.0 < EPSILON, "observation is not range=0");
    assert!(meas.observation()[(1, 0)].is_nan(), "observation is not range=0");
}
