extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

#[test]
fn eci() {
    use hifitime::instant::{Era, Instant};
    use na::Matrix3;
    use nyx::celestia::{CoordinateFrame, ECI};

    assert_eq!(
        ECI::to_inertial(Instant::from_precise_seconds(1.0, Era::Present)),
        Matrix3::identity()
    );
}

#[test]
fn theta_gmst() {
    // Vallado's example
    use hifitime::datetime::Datetime;
    use hifitime::TimeSystem;
    use nyx::celestia::ECEF;
    use std::f64::EPSILON;
    let dt = Datetime::new(1992, 8, 20, 12, 14, 0, 0).expect("wut?");
    assert!(
        (ECEF::theta_gmst(dt.into_instant()) - 152.5787878104796).abs() < EPSILON,
        "wrong θ GMST computed"
    );
}
