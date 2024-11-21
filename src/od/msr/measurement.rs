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

use super::MeasurementType;
use hifitime::Epoch;
use std::collections::HashMap;

/// A type-agnostic simultaneous measurement storage structure. Allows storing any number of simultaneous measurement of a given taker.
#[derive(Clone, Debug)]
pub struct Measurement {
    /// Tracker alias which made this measurement
    pub tracker: String,
    /// Epoch of the measurement
    pub epoch: Epoch,
    /// All measurements made simultaneously
    pub data: HashMap<MeasurementType, f64>,
}
