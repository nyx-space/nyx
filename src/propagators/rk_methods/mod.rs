/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

mod rk;
pub use self::rk::*;
mod dormand;
pub use self::dormand::*;
mod fehlberg;
pub use self::fehlberg::*;
mod verner;
pub use self::verner::*;

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
