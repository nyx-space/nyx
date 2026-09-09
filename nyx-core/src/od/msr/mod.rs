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

use std::str::FromStr;

use serde::{Deserialize, Serialize};

pub mod measurement;
pub mod sensitivity;
mod trackingdata;
pub mod two_way;
mod types;

pub use crate::od::ground_station::DopplerConfig;
pub use measurement::Measurement;
pub use trackingdata::TrackingDataArc;
pub use types::MeasurementType;

#[cfg(feature = "python")]
mod python;

#[cfg(feature = "python")]
use pyo3::prelude::*;

use crate::io::InputOutputError;

#[cfg_attr(feature = "python", pyclass(from_py_object))]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Default, Serialize, Deserialize, der::Enumerated)]
#[repr(u8)]
pub enum IntegrationRef {
    Start = 0,
    #[default]
    Middle = 1,
    End = 2,
}

impl FromStr for IntegrationRef {
    type Err = InputOutputError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.trim().to_lowercase().as_str() {
            "middle" => Ok(Self::Middle),
            "start" => Ok(Self::Start),
            "end" => Ok(Self::End),
            _ => Err(InputOutputError::UnsupportedData {
                which: format!("`{s}` is not a valid integration reference"),
            }),
        }
    }
}
