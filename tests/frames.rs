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
