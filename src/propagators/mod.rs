/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

extern crate rayon;
use std::fmt;

/// Provides different methods for controlling the error computation of the integrator.
pub mod error_ctrl;
pub use self::error_ctrl::*;

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

use crate::time::Duration;

/// The `RK` trait defines a Runge Kutta integrator.
#[allow(clippy::upper_case_acronyms)]
pub trait RK
where
    Self: Sized,
{
    /// Returns the order of this integrator (as u8 because there probably isn't an order greater than 255).
    /// The order is used for the adaptive step size only to compute the error between estimates.
    const ORDER: u8;

    /// Returns the stages of this integrator (as usize because it's used as indexing)
    const STAGES: usize;

    /// Returns a pointer to a list of f64 corresponding to the A coefficients of the Butcher table for that RK.
    /// This module only supports *implicit* integrators, and as such, `Self.a_coeffs().len()` must be of
    /// size (order+1)*(order)/2.
    /// *Warning:* this RK trait supposes that the implementation is consistent, i.e. c_i = \sum_j a_{ij}.
    const A_COEFFS: &'static [f64];
    /// Returns a pointer to a list of f64 corresponding to the b_i and b^*_i coefficients of the
    /// Butcher table for that RK. `Self.a_coeffs().len()` must be of size (order+1)*2.
    const B_COEFFS: &'static [f64];
}

/// Stores the details of the previous integration step of a given propagator. Access as `my_prop.clone().latest_details()`.
#[derive(Copy, Clone, Debug)]
pub struct IntegrationDetails {
    /// step size used
    pub step: Duration,
    /// error in the previous integration step
    pub error: f64,
    /// number of attempts needed by an adaptive step size to be within the tolerance
    pub attempts: u8,
}

impl fmt::Display for IntegrationDetails {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "IntegrationDetails {{step: {}, error: {:.3e}, attempts: {}}}",
            self.step, self.error, self.attempts
        )
    }
}
