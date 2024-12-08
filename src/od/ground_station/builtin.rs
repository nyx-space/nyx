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

use super::*;

impl GroundStation {
    pub fn dss65_madrid(
        elevation_mask: f64,
        range_noise_km: StochasticNoise,
        doppler_noise_km_s: StochasticNoise,
        iau_earth: Frame,
    ) -> Self {
        let mut measurement_types = IndexSet::new();
        measurement_types.insert(MeasurementType::Range);
        measurement_types.insert(MeasurementType::Doppler);

        let mut stochastics = IndexMap::new();
        stochastics.insert(MeasurementType::Range, range_noise_km);
        stochastics.insert(MeasurementType::Doppler, doppler_noise_km_s);

        Self {
            name: "Madrid".to_string(),
            elevation_mask_deg: elevation_mask,
            latitude_deg: 40.427_222,
            longitude_deg: 4.250_556,
            height_km: 0.834_939,
            frame: iau_earth,
            measurement_types,
            integration_time: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            stochastic_noises: Some(stochastics),
        }
    }

    pub fn dss34_canberra(
        elevation_mask: f64,
        range_noise_km: StochasticNoise,
        doppler_noise_km_s: StochasticNoise,
        iau_earth: Frame,
    ) -> Self {
        let mut measurement_types = IndexSet::new();
        measurement_types.insert(MeasurementType::Range);
        measurement_types.insert(MeasurementType::Doppler);

        let mut stochastics = IndexMap::new();
        stochastics.insert(MeasurementType::Range, range_noise_km);
        stochastics.insert(MeasurementType::Doppler, doppler_noise_km_s);

        Self {
            name: "Canberra".to_string(),
            elevation_mask_deg: elevation_mask,
            latitude_deg: -35.398_333,
            longitude_deg: 148.981_944,
            height_km: 0.691_750,
            frame: iau_earth,
            measurement_types,
            integration_time: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            stochastic_noises: Some(stochastics),
        }
    }

    pub fn dss13_goldstone(
        elevation_mask: f64,
        range_noise_km: StochasticNoise,
        doppler_noise_km_s: StochasticNoise,
        iau_earth: Frame,
    ) -> Self {
        let mut measurement_types = IndexSet::new();
        measurement_types.insert(MeasurementType::Range);
        measurement_types.insert(MeasurementType::Doppler);

        let mut stochastics = IndexMap::new();
        stochastics.insert(MeasurementType::Range, range_noise_km);
        stochastics.insert(MeasurementType::Doppler, doppler_noise_km_s);

        Self {
            name: "Goldstone".to_string(),
            elevation_mask_deg: elevation_mask,
            latitude_deg: 35.247_164,
            longitude_deg: 243.205,
            height_km: 1.071_149_04,
            frame: iau_earth,
            measurement_types,
            integration_time: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            stochastic_noises: Some(stochastics),
        }
    }
}
