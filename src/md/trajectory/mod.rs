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

pub(crate) mod spline;
mod traj;
mod traj_it;

pub(crate) use traj::interpolate;
pub use traj::Traj;

use super::StateParameter;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::time::{Duration, Epoch};
use crate::{NyxError, Orbit, Spacecraft, State};

use std::error::Error;
use std::fmt;

#[derive(Clone, PartialEq, Debug)]
pub enum TrajError {
    EventNotFound {
        start: Epoch,
        end: Epoch,
        event: String,
    },
    NoInterpolationData(Epoch),
    CreationError(String),
    OutOfSpline {
        req_epoch: Epoch,
        req_dur: Duration,
        spline_dur: Duration,
    },
}

impl fmt::Display for TrajError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::EventNotFound { start, end, event } => {
                write!(f, "Event {} not found between {} and {}", event, start, end)
            }
            Self::CreationError(reason) => write!(f, "Failed to create trajectory: {}", reason),
            Self::NoInterpolationData(e) => write!(f, "No interpolation data at {}", e),
            Self::OutOfSpline {
                req_epoch,
                req_dur,
                spline_dur,
            } => {
                write!(f, "Probable bug: Requested epoch {}, corresponding to an offset of {} in a spline of duration {}", req_epoch, req_dur, spline_dur)
            }
        }
    }
}

impl Error for TrajError {}

pub trait InterpState: State
where
    Self: Sized,
    DefaultAllocator: Allocator<f64, Self::Size>
        + Allocator<f64, Self::Size, Self::Size>
        + Allocator<f64, Self::VecLength>,
{
    /// Return the parameters in order
    fn params() -> Vec<StateParameter>;

    /// Return the requested parameter and its time derivative
    fn value_and_deriv(&self, param: &StateParameter) -> Result<(f64, f64), NyxError> {
        Ok((self.value(param)?, self.deriv(param)?))
    }

    /// Return the time derivative requested parameter
    fn deriv(&self, param: &StateParameter) -> Result<f64, NyxError> {
        Ok(self.value_and_deriv(param)?.1)
    }

    /// Sets the requested parameter
    fn set_value_and_deriv(
        &mut self,
        param: &StateParameter,
        value: f64,
        value_dt: f64,
    ) -> Result<(), NyxError>;
}

impl InterpState for Orbit {
    fn params() -> Vec<StateParameter> {
        vec![
            StateParameter::X,
            StateParameter::Y,
            StateParameter::Z,
            StateParameter::VX,
            StateParameter::VY,
            StateParameter::VZ,
        ]
    }
    fn value_and_deriv(&self, param: &StateParameter) -> Result<(f64, f64), NyxError> {
        match *param {
            StateParameter::X => Ok((self.x, self.vx)),
            StateParameter::Y => Ok((self.y, self.vy)),
            StateParameter::Z => Ok((self.z, self.vz)),
            StateParameter::VX => Ok((self.vx, 0.0)),
            StateParameter::VY => Ok((self.vy, 0.0)),
            StateParameter::VZ => Ok((self.vz, 0.0)),
            _ => Err(NyxError::StateParameterUnavailable),
        }
    }

    fn set_value_and_deriv(
        &mut self,
        param: &StateParameter,
        value: f64,
        _: f64,
    ) -> Result<(), NyxError> {
        match *param {
            StateParameter::X => {
                self.x = value;
            }
            StateParameter::Y => {
                self.y = value;
            }
            StateParameter::Z => {
                self.z = value;
            }
            StateParameter::VX => {
                self.vx = value;
            }
            StateParameter::VY => {
                self.vy = value;
            }
            StateParameter::VZ => {
                self.vz = value;
            }

            _ => return Err(NyxError::StateParameterUnavailable),
        }
        Ok(())
    }
}

impl InterpState for Spacecraft {
    fn params() -> Vec<StateParameter> {
        vec![
            StateParameter::X,
            StateParameter::Y,
            StateParameter::Z,
            StateParameter::VX,
            StateParameter::VY,
            StateParameter::VZ,
            StateParameter::FuelMass,
        ]
    }
    fn value_and_deriv(&self, param: &StateParameter) -> Result<(f64, f64), NyxError> {
        match *param {
            StateParameter::X => Ok((self.orbit.x, self.orbit.vx)),
            StateParameter::Y => Ok((self.orbit.y, self.orbit.vy)),
            StateParameter::Z => Ok((self.orbit.z, self.orbit.vz)),
            StateParameter::VX => Ok((self.orbit.vx, 0.0)),
            StateParameter::VY => Ok((self.orbit.vy, 0.0)),
            StateParameter::VZ => Ok((self.orbit.vz, 0.0)),
            StateParameter::FuelMass => Ok((self.fuel_mass_kg, 0.0)),
            _ => Err(NyxError::StateParameterUnavailable),
        }
    }

    fn set_value_and_deriv(
        &mut self,
        param: &StateParameter,
        value: f64,
        _: f64,
    ) -> Result<(), NyxError> {
        match *param {
            StateParameter::X => {
                self.orbit.x = value;
            }
            StateParameter::Y => {
                self.orbit.y = value;
            }
            StateParameter::Z => {
                self.orbit.z = value;
            }
            StateParameter::VX => {
                self.orbit.vx = value;
            }
            StateParameter::VY => {
                self.orbit.vy = value;
            }
            StateParameter::VZ => {
                self.orbit.vz = value;
            }
            StateParameter::Cr => self.cr = value,
            StateParameter::Cd => self.cd = value,
            StateParameter::FuelMass => self.fuel_mass_kg = value,
            _ => return Err(NyxError::StateParameterUnavailable),
        }
        Ok(())
    }
}
