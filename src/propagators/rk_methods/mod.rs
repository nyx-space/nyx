/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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
use std::str::FromStr;

use serde::Deserialize;
use serde::Serialize;

use crate::io::ConfigError;

use self::rk::*;
mod dormand;
use self::dormand::*;
mod verner;
use self::verner::*;

use super::PropagationError;

/// The `RK` trait defines a Runge Kutta integrator.
#[allow(clippy::upper_case_acronyms)]
trait RK
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

/// Enum of supported integration methods, all of which are part of the Runge Kutta family of ordinary differential equation (ODE) solvers.
/// Nomenclature: X-Y means that this is an X order solver with a Y order error correction step.
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize, Default)]
pub enum IntegratorMethod {
    /// Runge Kutta 8-9 is the recommended integrator for most application.
    #[default]
    RungeKutta89,
    /// `Dormand78` is a [Dormand-Prince integrator](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method). Coefficients taken from GMAT `src/base/propagator/PrinceDormand78.cpp`.
    DormandPrince78,
    /// `Dormand45` is a [Dormand-Prince integrator](https://en.wikipedia.org/wiki/Dormand%E2%80%93Prince_method).
    DormandPrince45,
    /// Runge Kutta 4 is a fixed step solver.
    RungeKutta4,
    /// Runge Kutta 4-5 [Cash Karp integrator](https://en.wikipedia.org/wiki/Cash%E2%80%93Karp_method).
    CashKarp45,
    /// Verner56 is an RK Verner integrator of order 5-6. Coefficients taken from [here (PDF)](http://people.math.sfu.ca/~jverner/classify.1992.ps).
    Verner56,
}

impl IntegratorMethod {
    /// Returns the order of this integrator (as u8 because there probably isn't an order greater than 255).
    /// The order is used for the adaptive step size only to compute the error between estimates.
    pub const fn order(self) -> u8 {
        match self {
            Self::RungeKutta89 => RK89::ORDER,
            Self::DormandPrince78 => Dormand78::ORDER,
            Self::DormandPrince45 => Dormand45::ORDER,
            Self::RungeKutta4 => RK4Fixed::ORDER,
            Self::CashKarp45 => CashKarp45::ORDER,
            Self::Verner56 => Verner56::ORDER,
        }
    }

    /// Returns the stages of this integrator, i.e. how many times the derivatives will be called
    pub const fn stages(self) -> usize {
        match self {
            Self::RungeKutta89 => RK89::STAGES,
            Self::DormandPrince78 => Dormand78::STAGES,
            Self::DormandPrince45 => Dormand45::STAGES,
            Self::RungeKutta4 => RK4Fixed::STAGES,
            Self::CashKarp45 => CashKarp45::STAGES,
            Self::Verner56 => Verner56::STAGES,
        }
    }

    /// Returns a pointer to a list of f64 corresponding to the A coefficients of the Butcher table for that RK.
    /// This module only supports *implicit* integrators, and as such, `Self.a_coeffs().len()` must be of
    /// size (order+1)*(order)/2.
    /// *Warning:* this RK trait supposes that the implementation is consistent, i.e. c_i = \sum_j a_{ij}.
    pub const fn a_coeffs(self) -> &'static [f64] {
        match self {
            Self::RungeKutta89 => RK89::A_COEFFS,
            Self::DormandPrince78 => Dormand78::A_COEFFS,
            Self::DormandPrince45 => Dormand45::A_COEFFS,
            Self::RungeKutta4 => RK4Fixed::A_COEFFS,
            Self::CashKarp45 => CashKarp45::A_COEFFS,
            Self::Verner56 => Verner56::A_COEFFS,
        }
    }
    /// Returns a pointer to a list of f64 corresponding to the b_i and b^*_i coefficients of the
    /// Butcher table for that RK. `Self.a_coeffs().len()` must be of size (order+1)*2.
    pub const fn b_coeffs(self) -> &'static [f64] {
        match self {
            Self::RungeKutta89 => RK89::B_COEFFS,
            Self::DormandPrince78 => Dormand78::B_COEFFS,
            Self::DormandPrince45 => Dormand45::B_COEFFS,
            Self::RungeKutta4 => RK4Fixed::B_COEFFS,
            Self::CashKarp45 => CashKarp45::B_COEFFS,
            Self::Verner56 => Verner56::B_COEFFS,
        }
    }
}

impl FromStr for IntegratorMethod {
    type Err = PropagationError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "rungekutta89" => Ok(Self::RungeKutta89),
            "dormandprince78" => Ok(Self::DormandPrince78),
            "dormandprince45" => Ok(Self::DormandPrince45),
            "rungekutta4" => Ok(Self::RungeKutta4),
            "cashkarp45" => Ok(Self::CashKarp45),
            "verner56" => Ok(Self::Verner56),
            _ => {
                let valid = [
                    "RungeKutta89",
                    "DormandPrince78",
                    "DormandPrince45",
                    "RungeKutta4",
                    "CashKarp45",
                    "Verner56",
                ];
                let valid_msg = valid.join(",");
                Err(PropagationError::PropConfigError {
                    source: ConfigError::InvalidConfig {
                        msg: format!("unknow integration method `{s}`, must be one of {valid_msg}"),
                    },
                })
            }
        }
    }
}

#[cfg(test)]
mod ut_propagator {
    use std::str::FromStr;

    use super::IntegratorMethod;

    #[test]
    fn from_str_ok() {
        let valid = [
            "RungeKutta89",
            "DormandPrince78",
            "DormandPrince45",
            "RungeKutta4",
            "CashKarp45",
            "Verner56",
        ];
        for method in valid {
            assert!(IntegratorMethod::from_str(method.to_uppercase().as_str()).is_ok());
        }
        assert!(IntegratorMethod::from_str("blah").is_err());
    }
}
