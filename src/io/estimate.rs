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

use anise::astro::orbit::Orbit;
use serde::{Deserialize, Serialize};

use super::{matrices::Matrix6Serde, ConfigRepr};

/// Enables serializing and deserializing of an orbit estimate.
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub struct OrbitEstimateSerde {
    /// Expected nominal orbit
    pub nominal: Orbit,
    /// Covariance _must_ be specified as in Cartesian format
    pub covar: Matrix6Serde,
}

impl ConfigRepr for OrbitEstimateSerde {}
