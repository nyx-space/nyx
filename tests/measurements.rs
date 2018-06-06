extern crate hifitime;
extern crate nyx_space as nyx;

use hifitime::julian::ModifiedJulian;
#[test]
fn nil_measurement() {
    use nyx::celestia::{State, ECEF};
    use nyx::od::Measurement;
    use nyx::od::ranging::{GroundStation, StdMeasurement};
    // Let's create a station and make it estimate the range and range rate of something which is strictly in the same spot.

    let lat = -7.906_635_7;
    let long = 345.5975;
    let height = 56.0e-3;

    let station = GroundStation::from_noise_values("nil", 0.0, lat, long, height, 0.0, 0.0);

    let at_station = State::<ECEF>::from_geodesic(lat, long, height, dt);

    let dt = ModifiedJulian::j2000();

    let meas = station.measure(at_station, dt);

    println!("{:?}", meas);
}
