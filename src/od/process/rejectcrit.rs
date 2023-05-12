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
use serde_derive::{Deserialize, Serialize};
use std::sync::Arc;

/// Defines how to reject measurements based on the value of their residual ratios.
/// The count is the minimum number of measurements to be ingested prior to applying the rejection criteria.
#[derive(Copy, Clone, Debug, Serialize, Deserialize)]
pub enum RejectCriteria {
    /// Reject measurements with a residual ratio greater than the provided value. This is a good default option when set to 3.0.
    ResidualRatio { count: usize, value: f64 },
    /// Reject measurements if their residual ratio is greater than the current z-score of the residual ratio multiplied by the provided value.
    ZScoreMultiplier { count: usize, value: f64 },
    /// Accept all measurements
    None,
}

impl ConfigRepr for RejectCriteria {}

impl Configurable for RejectCriteria {
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
