extern crate csv;
extern crate hifitime;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::nyx::celestia::State;
use self::nyx::od::ui::{Estimate, Filter, KfEstimate, NyxError, KF};

mod measurements;
mod multi_body;
mod robust;
mod spacecraft;
mod srif;
mod two_body;

use self::na::{Matrix2, Matrix2x6, Matrix6, Vector2};
use std::f64::EPSILON;

macro_rules! f64_nil {
    ($x:expr, $msg:expr) => {
        assert!($x.abs() < EPSILON, $msg)
    };
}

#[test]
fn empty_estimate() {
    let empty = KfEstimate::zeros(State::zeros());
    f64_nil!(
        empty.state_deviation.norm(),
        "expected state norm to be nil"
    );
    f64_nil!(empty.covar.norm(), "expected covar norm to be nil");
    f64_nil!(empty.stm.norm(), "expected STM norm to be nil");
    assert_eq!(empty.predicted, true, "expected predicted to be true");
}

#[test]
fn csv_serialize_empty_estimate() {
    use std::io;
    let empty = KfEstimate::zeros(State::zeros());
    let mut wtr = csv::Writer::from_writer(io::stdout());
    wtr.serialize(empty).expect("could not write to stdout");
}

#[test]
fn filter_errors() {
    let initial_estimate = KfEstimate::zeros(State::zeros());
    let measurement_noise = Matrix2::zeros();
    let real_obs = Vector2::zeros();
    let computed_obs = Vector2::zeros();
    let sensitivity = Matrix2x6::zeros();
    let stm = Matrix6::zeros();

    let mut ckf = KF::no_snc(initial_estimate, measurement_noise);
    match ckf.time_update(State::zeros()) {
        Ok(_) => panic!("expected the time update to fail"),
        Err(e) => {
            assert_eq!(e, NyxError::StateTransitionMatrixNotUpdated);
        }
    }
    match ckf.measurement_update(State::zeros(), real_obs, computed_obs) {
        Ok(_) => panic!("expected the measurement update to fail"),
        Err(e) => {
            assert_eq!(e, NyxError::StateTransitionMatrixNotUpdated);
        }
    }
    ckf.update_stm(stm);
    match ckf.measurement_update(State::zeros(), real_obs, computed_obs) {
        Ok(_) => panic!("expected the measurement update to fail"),
        Err(e) => {
            assert_eq!(e, NyxError::SensitivityNotUpdated);
        }
    }
    ckf.update_h_tilde(sensitivity);
    match ckf.measurement_update(State::zeros(), real_obs, computed_obs) {
        Ok(_) => panic!("expected the measurement update to fail"),
        Err(e) => {
            assert_eq!(e, NyxError::SingularKalmanGain);
        }
    }
}
