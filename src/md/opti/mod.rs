/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

pub mod convert_impulsive;
pub mod multipleshooting;
pub use multipleshooting::{ctrlnodes, multishoot};
/// Uses a Levenberg Marquardt minimizer to solve the damped least squares problem.
pub mod minimize_lm;
pub mod optimizer;
/// Uses a [Newton Raphson](https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization) method where the Jacobian is computed via finite differencing.
pub mod raphson_finite_diff;
/// Uses a [Newton Raphson](https://en.wikipedia.org/wiki/Newton%27s_method_in_optimization) method where the Jacobian is computed via hyperdual numbers.
pub mod raphson_hyperdual;
pub mod solution;
pub mod target_variable;

#[derive(Copy, Clone, Debug, PartialEq)]
pub enum DiffMethod {
    /// Slower, but more commonly used
    FiniteDiff,
    /// Significantly faster, but requires the automatic differentiation to be coded
    AutoDiff,
}
