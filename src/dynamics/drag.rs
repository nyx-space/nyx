use super::Dynamics;
use super::na::{U3, U6, Vector3, Vector6, VectorN};

/// `BasicDrag` implements the basic drag model as defined in Vallado, 4th ed., page 551, with an important caveat.
///
/// **WARNING:** This basic model assumes that the velocity of the spacecraft is identical to the velocity of the upper atmosphere,
/// This is a **bad** assumption. **Do not** use this model for high fidelity simulations.
#[derive(Copy, Clone)]
pub struct BasicDrag {
    pub rho: f64,  // atmospheric density in kg/m^3
    pub cd: f64,   // In Earth's atmosphere, this can be set to 2.2. Spheres are between 2.0 and 2.1.
    pub area: f64, // in m^2
    pub mass: f64, // in kg
}

impl Dynamics for BasicDrag {
    type StateSize = U6;
    /// NOTE: No state is associated with BasicDrag, always return zero time
    fn time(&self) -> f64 {
        0.0
    }

    /// NOTE: No state is associated with BasicDrag, always return zero
    fn state(&self) -> VectorN<f64, Self::StateSize> {
        Vector6::zeros()
    }

    /// NOTE: Nothing happens in this `set_state` since there is no state of the drag.
    fn set_state(&mut self, _new_t: f64, _new_state: &VectorN<f64, Self::StateSize>) {}

    /// This provides a **DELTA** of the state, which must be added to the result of the TwoBody propagator being used.
    /// However, the provided `state` must be the position and velocity.
    fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let velocity = state.fixed_rows::<U3>(3).into_owned();
        let drag = -0.5 * (self.cd * self.area / self.mass) * self.rho * velocity.norm() * velocity;
        Vector6::from_iterator(drag.iter().chain(Vector3::zeros().iter()).cloned())
    }
}
