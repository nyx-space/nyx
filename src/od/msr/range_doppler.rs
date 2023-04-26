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
use arrow::datatypes::{DataType, Field};
use hifitime::{Epoch, Unit};
use nalgebra::Matrix2x6;
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
        tx: Orbit,
        rx: Orbit,
        timestamp_noise_s: f64,
        range_noise_km: f64,
        doppler_noise_km_s: f64,
    ) -> Self {
        assert_eq!(tx.frame, rx.frame, "tx & rx in different frames");
        assert_eq!(tx.epoch, rx.epoch, "tx & rx states have different times");

        let range_vec_km = tx.radius() - rx.radius();
        let doppler_km_s = range_vec_km.dot(&(tx.velocity() - rx.velocity())) / range_vec_km.norm();

        Self {
            epoch: tx.epoch + timestamp_noise_s * Unit::Second,
            obs: Vector2::new(
                range_vec_km.norm() + range_noise_km,
                doppler_km_s + doppler_noise_km_s,
            ),
        }
    }

    /// Initialize a new two-way range and Doppler measurement from the provided states as times t_1 and t_2 and the effective noises.
    ///
    /// # Panics
    /// + If the epochs of the two states differ.
    /// + If the frames of the two states differ.
    /// + If both epochs are identical.
    pub fn two_way(
        tx: (Orbit, Orbit),
        rx: (Orbit, Orbit),
        timestamp_noise_s: f64,
        range_noise_km: f64,
        doppler_noise_km_s: f64,
    ) -> Self {
        assert_eq!(tx.0.frame, tx.1.frame, "both tx in different frames");
        assert_ne!(tx.0.epoch, tx.1.epoch, "tx states have identical times");
        assert_ne!(rx.0.epoch, rx.1.epoch, "rx states have identical times");
        assert_eq!(
            tx.0.epoch, rx.0.epoch,
            "tx and rx have different t_1: {} != {}",
            tx.0.epoch, rx.0.epoch
        );
        assert_eq!(
            tx.1.epoch, rx.1.epoch,
            "tx and rx have different t_1: {} != {}",
            tx.1.epoch, rx.1.epoch
        );

        // Compute the time difference.
        let delta_t = tx.1.epoch - tx.0.epoch;

        // Compute the one way range and range rate measurements where the noise is only apply to one of the measurements.
        let one_way_t_0 = Self::one_way(tx.0, rx.0, 0.0, 0.0, 0.0);
        let one_way_t_1 = Self::one_way(
            tx.1,
            rx.1,
            timestamp_noise_s,
            range_noise_km,
            doppler_noise_km_s,
        );

        // Compute the observation.
        let obs = (one_way_t_1.obs - one_way_t_0.obs) / delta_t.to_seconds();

        Self {
            epoch: one_way_t_1.epoch - delta_t * 0.5,
            obs,
        }
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

impl EstimateFrom<Spacecraft, RangeDoppler> for Orbit {
    fn extract(from: Spacecraft) -> Self {
        from.orbit
    }

    fn sensitivity(
        msr: &RangeDoppler,
        receiver: Self,
        transmitter: Self,
    ) -> OMatrix<f64, <RangeDoppler as Measurement>::MeasurementSize, Self::Size>
    where
        DefaultAllocator:
            Allocator<f64, <RangeDoppler as Measurement>::MeasurementSize, Self::Size>,
    {
        <Orbit as EstimateFrom<Orbit, RangeDoppler>>::sensitivity(msr, receiver, transmitter)
    }
}

impl EstimateFrom<Orbit, RangeDoppler> for Orbit {
    fn extract(from: Orbit) -> Self {
        from
    }

    fn sensitivity(
        msr: &RangeDoppler,
        receiver: Self,
        transmitter: Self,
    ) -> OMatrix<f64, <RangeDoppler as Measurement>::MeasurementSize, Self::Size>
    where
        DefaultAllocator:
            Allocator<f64, <RangeDoppler as Measurement>::MeasurementSize, Self::Size>,
    {
        let delta_r = receiver.radius() - transmitter.radius();
        let delta_v = receiver.velocity() - transmitter.velocity();
        let ρ = msr.observation()[0];
        let ρ_dot = msr.observation()[1];
        let m11 = delta_r.x / ρ;
        let m12 = delta_r.y / ρ;
        let m13 = delta_r.z / ρ;
        let m21 = delta_v.x / ρ - ρ_dot * delta_r.x / ρ.powi(2);
        let m22 = delta_v.y / ρ - ρ_dot * delta_r.y / ρ.powi(2);
        let m23 = delta_v.z / ρ - ρ_dot * delta_r.z / ρ.powi(2);

        Matrix2x6::new(m11, m12, m13, 0.0, 0.0, 0.0, m21, m22, m23, m11, m12, m13)
    }
}

impl EstimateFrom<Spacecraft, RangeDoppler> for Spacecraft {
    fn extract(from: Spacecraft) -> Self {
        from
    }

    fn sensitivity(
        _msr: &RangeDoppler,
        _receiver: Self,
        _transmitter: Orbit,
    ) -> OMatrix<f64, <RangeDoppler as Measurement>::MeasurementSize, Self::Size>
    where
        DefaultAllocator:
            Allocator<f64, <RangeDoppler as Measurement>::MeasurementSize, Self::Size>,
    {
        todo!("cannot yet estimate a full spacecraft state")
    }
}
