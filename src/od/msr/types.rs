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

use anise::astro::AzElRange;
use arrow::datatypes::{DataType, Field};
use serde_derive::{Deserialize, Serialize};
use std::collections::HashMap;

use crate::od::ODError;

#[derive(Copy, Clone, Debug, Hash, Serialize, Deserialize, PartialEq, Eq)]
pub enum MeasurementType {
    #[serde(rename = "range_km")]
    Range,
    #[serde(rename = "doppler_km_s")]
    Doppler,
    #[serde(rename = "azimuth_deg")]
    Azimuth,
    #[serde(rename = "elevation_deg")]
    Elevation,
}

impl MeasurementType {
    /// Returns the expected unit of this measurement type
    pub fn unit(self) -> &'static str {
        match self {
            Self::Range => "km",
            Self::Doppler => "km/s",
            Self::Azimuth | Self::Elevation => "deg",
        }
    }

    /// Returns the fields for this kind of measurement. The metadata includes a `unit` field with the unit.
    /// Column is nullable in case there is no such measurement at a given epoch.
    pub fn to_field(&self) -> Field {
        let mut meta = HashMap::new();
        meta.insert("unit".to_string(), self.unit().to_string());

        Field::new(
            format!("{self:?} ({})", self.unit()),
            DataType::Float64,
            true,
        )
        .with_metadata(meta)
    }

    /// Computes the one way measurement from an AER object and the noise of this measurement type, returned in the units of this measurement type.
    pub fn compute_one_way(self, aer: AzElRange, noise: f64) -> Result<f64, ODError> {
        match self {
            Self::Range => Ok(aer.range_km + noise),
            Self::Doppler => Ok(aer.range_rate_km_s + noise),
            Self::Azimuth => Ok(aer.azimuth_deg + noise),
            Self::Elevation => Ok(aer.elevation_deg + noise),
        }
    }

    /// Computes the two way measurement from two AER values and the noise of this measurement type, returned in the units of this measurement type.
    /// Two way is modeled by averaging the measurement in between both times, and adding the noise divided by sqrt(2).
    pub fn compute_two_way(
        self,
        aer_t0: AzElRange,
        aer_t1: AzElRange,
        noise: f64,
    ) -> Result<f64, ODError> {
        match self {
            Self::Range => {
                let range_km = (aer_t1.range_km + aer_t0.range_km) * 0.5;
                Ok(range_km + noise / 2.0_f64.sqrt())
            }
            Self::Doppler => {
                let doppler_km_s = (aer_t1.range_rate_km_s + aer_t0.range_rate_km_s) * 0.5;
                Ok(doppler_km_s + noise / 2.0_f64.sqrt())
            }
            Self::Azimuth => {
                let az_deg = (aer_t1.azimuth_deg + aer_t0.azimuth_deg) * 0.5;
                Ok(az_deg + noise / 2.0_f64.sqrt())
            }
            Self::Elevation => {
                let el_deg = (aer_t1.elevation_deg + aer_t0.elevation_deg) * 0.5;
                Ok(el_deg + noise / 2.0_f64.sqrt())
            }
        }
    }
}
