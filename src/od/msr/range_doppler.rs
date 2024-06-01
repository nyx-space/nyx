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

use crate::cosmic::Orbit;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, OMatrix, OVector, Vector2, U2};
use crate::od::msr::RangeMsr;
use crate::od::{EstimateFrom, Measurement};
use crate::{Spacecraft, TimeTagged};
use anise::astro::AzElRange;
use arrow::datatypes::{DataType, Field};
use hifitime::{Epoch, Unit};
use std::collections::HashMap;

/// A simultaneous range and Doppler measurement in units of km and km/s, available both in one way and two way measurement.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RangeDoppler {
    /// Epoch of the observation
    pub epoch: Epoch,
    /// Observation vector in km and km/s
    pub obs: Vector2<f64>,
}

impl RangeDoppler {
    /// Initialize a new one-way range and Doppler measurement from the provided states and the effective noises.
    ///
    /// # Panics
    /// + If the epochs of the two states differ.
    /// + If the frames of the two states differ.
    pub fn one_way(
        aer: AzElRange,
        timestamp_noise_s: f64,
        range_noise_km: f64,
        doppler_noise_km_s: f64,
    ) -> Self {
        Self {
            epoch: aer.epoch + timestamp_noise_s * Unit::Second,
            obs: Vector2::new(
                aer.range_km + range_noise_km,
                aer.range_rate_km_s + doppler_noise_km_s,
            ),
        }
    }

    /// Initialize a new two-way range and Doppler measurement from the provided states as times t_1 and t_2 and the effective noises.
    ///
    /// The measurement is time-tagged at realization, i.e. at the end of the integration time (plus timestamp noise).
    ///
    /// # Noise
    /// The measurements are not considered to be independent distributed variables. As such, the noises are reduced by a factor of sqrt(2).
    ///
    /// # Panics
    /// + If the epochs of the two states differ.
    /// + If the frames of the two states differ.
    /// + If both epochs are identical.
    pub fn two_way(
        aer_t0: AzElRange,
        aer_t1: AzElRange,
        timestamp_noise_s: f64,
        range_noise_km: f64,
        doppler_noise_km_s: f64,
    ) -> Self {
        assert_eq!(
            aer_t0.epoch, aer_t1.epoch,
            "AER data have different t_1: {} != {}",
            aer_t0.epoch, aer_t1.epoch
        );

        let range_km = (aer_t1.range_km - aer_t0.range_km) * 0.5;
        let doppler_km_s = (aer_t1.range_rate_km_s - aer_t0.range_rate_km_s) * 0.5;

        // Time tagged at the realization of this measurement, i.e. at the end of the integration time.
        let epoch = aer_t1.epoch + timestamp_noise_s * Unit::Second;

        let obs = Vector2::new(
            range_km + range_noise_km / 2.0_f64.sqrt(),
            doppler_km_s + doppler_noise_km_s / 2.0_f64.sqrt(),
        );

        debug!(
            "two way msr @ {epoch}:\naer_t0 = {}\naer_t1 = {}{obs}",
            aer_t0, aer_t1
        );

        Self { epoch, obs }
    }
}

impl TimeTagged for RangeDoppler {
    fn epoch(&self) -> Epoch {
        self.epoch
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.epoch = epoch
    }
}

impl Measurement for RangeDoppler {
    type MeasurementSize = U2;

    /// Returns this measurement as a vector of Range and Range Rate
    ///
    /// **Units:** km, km/s
    fn observation(&self) -> Vector2<f64> {
        self.obs
    }

    fn fields() -> Vec<Field> {
        let mut meta = HashMap::new();
        meta.insert("unit".to_string(), "km/s".to_string());

        vec![
            RangeMsr::fields()[0].clone(),
            Field::new("Doppler (km/s)", DataType::Float64, false).with_metadata(meta),
        ]
    }

    fn from_observation(epoch: Epoch, obs: OVector<f64, Self::MeasurementSize>) -> Self {
        Self { epoch, obs }
    }
}

impl EstimateFrom<Spacecraft, RangeDoppler> for Spacecraft {
    fn extract(from: Spacecraft) -> Self {
        from
    }

    fn sensitivity(
        msr: &RangeDoppler,
        receiver: Self,
        transmitter: Orbit,
    ) -> OMatrix<f64, <RangeDoppler as Measurement>::MeasurementSize, Self::Size>
    where
        DefaultAllocator:
            Allocator<f64, <RangeDoppler as Measurement>::MeasurementSize, Self::Size>,
    {
        let delta_r = receiver.orbit.radius_km - transmitter.radius_km;
        let delta_v = receiver.orbit.velocity_km_s - transmitter.velocity_km_s;
        let ρ = msr.observation()[0];
        let ρ_dot = msr.observation()[1];
        let m11 = delta_r.x / ρ;
        let m12 = delta_r.y / ρ;
        let m13 = delta_r.z / ρ;
        let m21 = delta_v.x / ρ - ρ_dot * delta_r.x / ρ.powi(2);
        let m22 = delta_v.y / ρ - ρ_dot * delta_r.y / ρ.powi(2);
        let m23 = delta_v.z / ρ - ρ_dot * delta_r.z / ρ.powi(2);

        let items = &[
            m11, m12, m13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, m21, m22, m23, m11, m12, m13, 0.0, 0.0,
            0.0,
        ];

        OMatrix::<f64, <RangeDoppler as Measurement>::MeasurementSize, Self::Size>::from_row_slice(
            items,
        )
    }
}
