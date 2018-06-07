extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::hifitime::instant::*;
use self::hifitime::julian::*;
use self::hifitime::SECONDS_PER_DAY;
use self::na::Vector6;
use self::nyx::celestia::{State, EARTH, ECI};
use self::nyx::dynamics::celestial::{TwoBody, TwoBodyWithStm};
use self::nyx::dynamics::Dynamics;
use self::nyx::od::kalman::{Estimate, FilterError, KF};
use self::nyx::od::ranging::{GroundStation, StdMeasurement};
use self::nyx::od::{Linearization, Measurement};
use self::nyx::propagators::{error_ctrl, Options, Propagator, RK4Fixed};
use std::sync::mpsc;
use std::sync::mpsc::{Receiver, Sender};

#[test]
fn fixed_step_perfect_stations() {
    use std::thread;
    // Tests that we can generate measurements on one side and get a proper estimate on the other.

    // Define the ground stations.
    let dss65_madrid = GroundStation::from_noise_values("Madrid", 0.0, 40.427222, 4.250556, 0.834939, 0.0, 0.0);
    let dss34_canberra = GroundStation::from_noise_values("Canberra", 0.0, -35.398333, 148.981944, 0.691750, 0.0, 0.0);
    let dss13_goldstone = GroundStation::from_noise_values("Goldstone", 0.0, 35.247164, 243.205, 1.07114904, 0.0, 0.0);

    // Define the propagator information.
    let prop_time = SECONDS_PER_DAY;
    let step_size = 10.0;
    let opts = Options::with_fixed_step(step_size);

    // Define the information for the channels.
    let (truth_tx, truth_rx): (Sender<(f64, Vector6<f64>)>, Receiver<(f64, Vector6<f64>)>) = mpsc::channel();
    let (meas_tx, meas_rx): (Sender<StdMeasurement>, Receiver<StdMeasurement>) = mpsc::channel();

    // Define state information.
    let dt = ModifiedJulian { days: 21545.0 };
    let initial_state = State::from_cartesian_eci(-2436.45, -2436.45, 6891.037, 5.088611, -5.088611, 0.0, dt);

    // Generate the truth data on one thread.
    thread::spawn(move || {
        let mut prop = Propagator::new::<RK4Fixed>(&opts);
        let mut dyn = TwoBody::from_state_vec::<EARTH>(initial_state.to_cartesian_vec());
        dyn.tx_chan = Some(&truth_tx);
        let (final_t, final_state_vec) = prop.until_time_elapsed(prop_time, &mut dyn, error_ctrl::rss_step_pos_vel);
        dyn.set_state(final_t, &final_state_vec);
    });

    // Receive the states on the main thread, and populate the measurement channel.
    loop {
        if let Ok((t, state_vec)) = truth_rx.recv() {
            // Convert the state to ECI.
            let this_dt =
                ModifiedJulian::from_instant(dt.into_instant() + Instant::from_precise_seconds(t, Era::Present).duration());
            let rx_state = State::from_cartesian_vec::<EARTH, ModifiedJulian>(&state_vec, this_dt, ECI {});
            let meas_dss13 = dss13_goldstone.measure(rx_state, this_dt.into_instant());
            if meas_dss13.visible() {}
        }
    }
}
