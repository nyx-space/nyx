/// Provides different methods for controlling the error computation of the integrator.
pub mod error_ctrl;

pub mod events;

// Re-Export
mod rk;
pub use self::rk::*;
mod dormand;
pub use self::dormand::*;
mod fehlberg;
pub use self::fehlberg::*;
mod verner;
pub use self::verner::*;
mod propagator;
pub use self::propagator::*;

/// The `RK` trait defines a Runge Kutta integrator.
pub trait RK
where
    Self: Sized,
{
    /// Returns the order of this integrator (as u8 because there probably isn't an order greater than 255).
    /// The order is used for the adaptive step size only to compute the error between estimates.
    fn order() -> u8;

    /// Returns the stages of this integrator (as usize because it's used as indexing)
    fn stages() -> usize;

    /// Returns a pointer to a list of f64 corresponding to the A coefficients of the Butcher table for that RK.
    /// This module only supports *implicit* integrators, and as such, `Self.a_coeffs().len()` must be of
    /// size (order+1)*(order)/2.
    /// *Warning:* this RK trait supposes that the implementation is consistent, i.e. c_i = \sum_j a_{ij}.
    fn a_coeffs() -> &'static [f64];
    /// Returns a pointer to a list of f64 corresponding to the b_i and b^*_i coefficients of the
    /// Butcher table for that RK. `Self.a_coeffs().len()` must be of size (order+1)*2.
    fn b_coeffs() -> &'static [f64];
}

/// Stores the details of the previous integration step of a given propagator. Access as `my_prop.clone().latest_details()`.
#[derive(Clone, Debug)]
pub struct IntegrationDetails {
    /// step size used
    pub step: f64,
    /// error in the previous integration step
    pub error: f64,
    /// number of attempts needed by an adaptive step size to be within the tolerance
    pub attempts: u8,
}
