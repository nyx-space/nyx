extern crate nyx_space as nyx;
extern crate rand;
extern crate rand_distr;

use rand::thread_rng;
use rand_distr::{Distribution, Normal};

use nyx::dynamics::Dynamics;
use nyx::linalg::{Const, OMatrix, OVector};
use nyx::time::Epoch;
use nyx::{NyxError, State};

use std::fmt;

#[derive(Copy, Clone, PartialEq)]
struct BiasDriftState {
    epoch: Epoch,
    bias: f64,
    drift: f64,
}

impl fmt::Display for BiasDriftState {
    // Prints as Cartesian in floating point with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "bias = {}\tdrift = {}",
            format!("{:.*}", decimals, self.bias),
            format!("{:.*}", decimals, self.drift),
        )
    }
}

impl fmt::LowerExp for BiasDriftState {
    // Prints as Cartesian in scientific notation with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "bias = {}\tdrift = {}",
            format!("{:.*e}", decimals, self.bias),
            format!("{:.*e}", decimals, self.drift),
        )
    }
}

impl State for BiasDriftState {
    type Size = Const<2>;
    type VecLength = Const<2>;

    fn zeros() -> Self {
        Self {
            epoch: Epoch::from_jde_et(0.0),
            bias: 0.0,
            drift: 0.0,
        }
    }

    fn as_vector(&self) -> Result<OVector<f64, Self::VecLength>, NyxError> {
        Ok(OVector::<f64, Const<2>>::new(self.bias, self.drift))
    }

    fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, NyxError> {
        // NOTE: The STM is dependent on the model here, not on the state!
        unimplemented!()
    }

    fn reset_stm(&mut self) {}

    fn set(
        &mut self,
        epoch: Epoch,
        vector: &OVector<f64, Self::VecLength>,
    ) -> Result<(), NyxError> {
        self.set_epoch(epoch);
        self.bias = vector[0];
        self.drift = vector[1];
        println!("{}", self);
        Ok(())
    }

    fn epoch(&self) -> Epoch {
        self.epoch
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.epoch = epoch;
    }

    fn add(self, other: OVector<f64, Self::Size>) -> Self {
        let mut me = self;
        me.bias += other[0];
        me.drift += other[1];

        me
    }
}

#[derive(Clone)]
struct BiasDriftCoupledGM {
    tau: f64,
    omega: f64,
    zeta: f64,
}

impl Dynamics for BiasDriftCoupledGM {
    type StateType = BiasDriftState;
    type HyperdualSize = Const<0>;

    fn eom(
        &self,
        _delta_t: f64,
        state_vec: &OVector<f64, <Self::StateType as State>::VecLength>,
        _state_ctx: &Self::StateType,
    ) -> Result<OVector<f64, <Self::StateType as State>::VecLength>, NyxError> {
        let a_mat = OMatrix::<f64, Const<2>, Const<2>>::new(
            -1.0 / self.tau,
            1.0,
            -self.omega.powi(2),
            -2.0 * self.zeta * self.omega,
        );

        // let a_mat = OMatrix::<f64, Const<2>, Const<2>>::new(0.0, 1.0, 0.0, 0.0);

        // Ignore the noises for now -- would need a mut thread_rng for reproducibility
        let b_mat = OMatrix::<f64, Const<2>, Const<2>>::identity();

        let d = Normal::new(0.0, 0.9).unwrap();

        let w_vec =
            OVector::<f64, Const<2>>::new(d.sample(&mut thread_rng()), d.sample(&mut thread_rng()));

        Ok(a_mat * state_vec + b_mat * w_vec)
    }

    #[allow(clippy::type_complexity)]
    fn dual_eom(
        &self,
        _delta_t: f64,
        _osculating_state: &Self::StateType,
    ) -> Result<
        (
            OVector<f64, <Self::StateType as State>::Size>,
            OMatrix<f64, <Self::StateType as State>::Size, <Self::StateType as State>::Size>,
        ),
        NyxError,
    > {
        unimplemented!()
    }
}

#[test]
fn gm_test() {
    use nyx::propagators::Propagator;
    use nyx::time::TimeUnit;
    use std::sync::Arc;
    let init = BiasDriftState::zeros();
    let model = BiasDriftCoupledGM {
        tau: 900.0,
        omega: 2.0,
        zeta: 0.1,
    };
    let fstate = Propagator::default(Arc::new(model))
        .with(init)
        .for_duration(60 * TimeUnit::Second)
        .unwrap();
    println!("Final state = {}", fstate);
}
