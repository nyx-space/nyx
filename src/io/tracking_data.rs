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

use crate::{
    linalg::{allocator::Allocator, DefaultAllocator, DimName, OVector},
    State,
};
use hifitime::Epoch;
use serde::{Deserialize, Serialize};

use crate::od::SimMeasurement;

#[derive(Debug, Serialize, Deserialize)]
pub struct DynamicMeasurement {
    measurement_type: String,
    epoch: Epoch,
    data: Vec<f64>,
}

impl DynamicMeasurement {
    /// Initializes a new dynamic measurement from a concrete measurement.
    /// TODO: Is this ever needed? The measurements will be serialized directly, so I'm not sure I'll ever need this function.
    pub fn new<M: SimMeasurement>(measurement: M) -> DynamicMeasurement
    where
        Self: Sized,
        DefaultAllocator: Allocator<f64, M::MeasurementSize>
            + Allocator<f64, <M::State as State>::Size>
            + Allocator<f64, <M::State as State>::Size, <M::State as State>::Size>
            + Allocator<f64, <M::State as State>::VecLength>
            + Allocator<f64, M::MeasurementSize, <M::State as State>::Size>,
    {
        let mut data = Vec::new();
        for obs in measurement.observation().iter() {
            data.push(*obs);
        }
        // for data in measurement.
        DynamicMeasurement {
            measurement_type: std::any::type_name::<M>().to_string(),
            data,
            epoch: measurement.epoch(),
        }
    }

    /// Attempts to convert this dynamic measurement into the concrete turbofish measurement
    pub fn try_into<M: SimMeasurement>(self) -> Result<M, Box<dyn std::error::Error>>
    where
        Self: Sized,
        DefaultAllocator: Allocator<f64, M::MeasurementSize>
            + Allocator<f64, <M::State as State>::Size>
            + Allocator<f64, <M::State as State>::Size, <M::State as State>::Size>
            + Allocator<f64, <M::State as State>::VecLength>
            + Allocator<f64, M::MeasurementSize, <M::State as State>::Size>,
    {
        let expected_type = std::any::type_name::<M>().to_string();
        if self.measurement_type != expected_type {
            return Err(format!(
                "Expected measurement of type {}, but got {}",
                expected_type, self.measurement_type
            )
            .into());
        }
        // Try to build the vector of the observation
        if self.data.len() != M::MeasurementSize::USIZE {
            panic!();
        }
        let obs = OVector::<f64, M::MeasurementSize>::from_iterator(self.data);
        let measurement = M::from_observation(self.epoch, obs);
        Ok(measurement)
    }
}
