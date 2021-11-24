extern crate csv;
extern crate nalgebra as na;
extern crate nyx_space as nyx;

use self::nyx::cosmic::Orbit;
use self::nyx::od::ui::{Estimate, Filter, KfEstimate, NyxError, KF};
use self::nyx::State;

mod measurements;
mod multi_body;
mod robust;
mod spacecraft;
mod two_body;
mod xhat_dev;

use self::nyx::linalg::{Matrix2, Matrix2x6, Vector2};
use std::f64::EPSILON;

macro_rules! f64_nil {
    ($x:expr, $msg:expr) => {
        assert!($x.abs() < EPSILON, $msg)
    };
}

#[test]
fn empty_estimate() {
    let empty = KfEstimate::zeros(Orbit::zeros());
    f64_nil!(
        empty.state_deviation.norm(),
        "expected state norm to be nil"
    );
    f64_nil!(empty.covar.norm(), "expected covar norm to be nil");
    f64_nil!(
        empty.stm.diagonal().norm() - 6.0_f64.sqrt(),
        "expected STM norm to be sqrt(dim(STM))"
    );
    assert!(empty.predicted, "expected predicted to be true");
}

#[test]
fn csv_serialize_empty_estimate() {
    use std::io;
    let empty = KfEstimate::zeros(Orbit::zeros());
    let mut wtr = csv::Writer::from_writer(io::stdout());
    wtr.serialize(empty).expect("could not write to stdout");
}

#[test]
fn filter_errors() {
    let initial_estimate = KfEstimate::zeros(Orbit::zeros());
    let measurement_noise = Matrix2::zeros();
    let real_obs = &Vector2::zeros();
    let computed_obs = &Vector2::zeros();
    let sensitivity = Matrix2x6::zeros();

    let mut ckf = KF::no_snc(initial_estimate, measurement_noise);
    match ckf.measurement_update(Orbit::zeros(), real_obs, computed_obs) {
        Ok(_) => panic!("expected the measurement update to fail"),
        Err(e) => {
            assert_eq!(e, NyxError::SensitivityNotUpdated);
        }
    }
    ckf.update_h_tilde(sensitivity);
    match ckf.measurement_update(Orbit::zeros(), real_obs, computed_obs) {
        Ok(_) => panic!("expected the measurement update to fail"),
        Err(e) => {
            assert_eq!(e, NyxError::SingularKalmanGain);
        }
    }
}
