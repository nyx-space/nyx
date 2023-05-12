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

use crate::cosmic::Cosm;
use crate::io::{ConfigError, ConfigRepr, Configurable};
#[cfg(feature = "python")]
use pyo3::prelude::*;
use serde_derive::{Deserialize, Serialize};
use std::sync::Arc;

/// Reject measurements with a residual ratio greater than the provided value.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "python", pyclass)]
pub struct FltResid {
    /// Minimum number of accepted measurements before applying the rejection criteria.
    pub min_accepted: usize,
    /// Number of sigmas for a measurement to be considered an outlier.
    pub num_sigmas: f64,
}

#[cfg(feature = "python")]
impl FltResid {
    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }
}

impl Default for FltResid {
    fn default() -> Self {
        Self {
            min_accepted: 10,
            num_sigmas: 3.0,
        }
    }
}

impl ConfigRepr for FltResid {}

impl Configurable for FltResid {
    type IntermediateRepr = Self;

    fn from_config(cfg: Self, _cosm: Arc<Cosm>) -> Result<Self, ConfigError>
    where
        Self: Sized,
    {
        Ok(cfg)
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, ConfigError> {
        Ok(*self)
    }
}
