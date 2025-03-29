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

use crate::io::ConfigRepr;
use serde_derive::{Deserialize, Serialize};

/// Reject measurements if the prefit is greater than the provided sigmas deviation from the measurement noise.
///
/// # Important
/// Some software, like ODTK, processes each measurement as a scalar. Nyx can process the measurements together.
/// As such, if the prefit on range is bad, then the Doppler measurement with the same time stamp will also be rejected.
/// This can lead to better convergence of the filter, and more appropriate results.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct ResidRejectCrit {
    /// Number of sigmas for a measurement to be considered an outlier.
    pub num_sigmas: f64,
}

impl Default for ResidRejectCrit {
    /// By default, a measurement is rejected if its prefit residual is greater the 3-sigma value of the measurement noise at that time step.
    /// This corresponds to [1 chance in in 370](https://en.wikipedia.org/wiki/68%E2%80%9395%E2%80%9399.7_rule).
    fn default() -> Self {
        Self { num_sigmas: 3.0 }
    }
}

impl ConfigRepr for ResidRejectCrit {}
