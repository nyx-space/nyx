extern crate hifitime;
extern crate nyx_space as nyx;

#[test]
fn nil_measurement() {
    use hifitime::julian::*;
    use nyx::celestia::{State, ECEF};
    use nyx::od::Measurement;
    use nyx::od::ranging::{GroundStation, StdMeasurement};
    // Let's create a station and make it estimate the range and range rate of something which is strictly in the same spot.

    let lat = -7.906_635_7;
    let long = 345.5975;
    let height = 56.0e-3;
    let dt = ModifiedJulian::j2000();

    let station = GroundStation::from_noise_values("nil".to_string(), 0.0, lat, long, height, 0.0, 0.0);

    let at_station = State::<ECEF>::from_geodesic(lat, long, height, dt);

    let meas = station.measure(at_station, dt.into_instant());

    println!("{:?}", meas);
}
