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
use indexmap::{IndexMap, IndexSet};
use nalgebra::{allocator::Allocator, DefaultAllocator, DimName, OVector};
use std::fmt;

/// A type-agnostic simultaneous measurement storage structure. Allows storing any number of simultaneous measurement of a given taker.
#[derive(Clone, Debug, PartialEq)]
pub struct Measurement {
    /// Tracker alias which made this measurement
    pub tracker: String,
    /// Epoch of the measurement
    pub epoch: Epoch,
    /// All measurements made simultaneously
    pub data: IndexMap<MeasurementType, f64>,
}

impl Measurement {
    pub fn new(tracker: String, epoch: Epoch) -> Self {
        Self {
            tracker,
            epoch,
            data: IndexMap::new(),
        }
    }

    pub fn push(&mut self, msr_type: MeasurementType, msr_value: f64) {
        self.data.insert(msr_type, msr_value);
    }

    pub fn with(mut self, msr_type: MeasurementType, msr_value: f64) -> Self {
        self.push(msr_type, msr_value);
        self
    }

    /// Builds an observation vector for this measurement provided a set of measurement types.
    /// If the requested measurement type is not available, then that specific row is set to zero.
    /// The caller must set the appropriate sensitivity matrix rows to zero.
    pub fn observation<S: DimName>(&self, types: &IndexSet<MeasurementType>) -> OVector<f64, S>
    where
        DefaultAllocator: Allocator<S>,
    {
        let mut obs = OVector::zeros();
        for (i, t) in types.iter().enumerate() {
            if let Some(msr_value) = self.data.get(t) {
                obs[i] = *msr_value;
            }
        }
        obs
    }

    /// Returns a vector specifying which measurement types are available.
    pub fn availability(&self, types: &IndexSet<MeasurementType>) -> Vec<bool> {
        let mut rtn = Vec::with_capacity(types.len());
        for (i, t) in types.iter().enumerate() {
            if self.data.contains_key(t) {
                rtn[i] = true;
            }
        }
        rtn
    }
}

impl fmt::Display for Measurement {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let msrs = self
            .data
            .iter()
            .map(|(msr_type, msr_value)| format!("{msr_type:?} = {msr_value} {}", msr_type.unit()))
            .collect::<Vec<String>>()
            .join(", ");

        write!(f, "{} measured {} on {}", self.tracker, msrs, self.epoch)
    }
}
