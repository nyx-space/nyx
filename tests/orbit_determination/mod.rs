extern crate nalgebra as na;
extern crate nyx_space as nyx;

use nyx::od::ODError;
use nyx::Spacecraft;

use self::nyx::od::prelude::{Estimate, Filter, KalmanVariant, KfEstimate, KF};
use self::nyx::State;

mod measurements;
mod multi_body;
mod resid_reject;
mod robust;
mod robust_az_el;
mod simulator;
mod spacecraft;
mod trackingarc;
mod two_body;
mod xhat_dev;

use self::nyx::linalg::{Matrix2, SMatrix, Vector2};

macro_rules! f64_nil {
    ($x:expr, $msg:expr) => {
        assert!($x.abs() < f64::EPSILON, $msg)
    };
}

#[test]
fn empty_estimate() {
    let empty = KfEstimate::zeros(Spacecraft::zeros().with_stm());
    f64_nil!(
        empty.state_deviation.norm(),
        "expected state norm to be nil"
    );
    f64_nil!(empty.covar.norm(), "expected covar norm to be nil");
    f64_nil!(
        empty.stm.diagonal().norm() - 9.0_f64.sqrt(),
        "expected STM norm to be sqrt(dim(STM))"
    );
    assert!(empty.predicted, "expected predicted to be true");
}

#[test]
fn filter_errors() {
    let initial_estimate = KfEstimate::zeros(Spacecraft::zeros());
    let measurement_noise = Matrix2::zeros();
    let real_obs = Vector2::zeros();
    let computed_obs = Vector2::zeros();
    let sensitivity = SMatrix::<f64, 2, 9>::zeros();

    let mut ckf = KF::new(initial_estimate, KalmanVariant::DeviationTracking);

    match ckf.measurement_update(
        Spacecraft::zeros().with_stm(),
        real_obs,
        computed_obs,
        measurement_noise,
        sensitivity,
        None,
    ) {
        Ok(_) => panic!("expected the measurement update to fail"),
        Err(e) => {
            assert_eq!(e, ODError::SingularNoiseRk);
        }
    }
}
