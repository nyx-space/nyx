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

use std::collections::HashMap;

use hifitime::{Epoch, Unit};

use crate::cosmic::Orbit;
use crate::linalg::{DimName, Matrix2x6, OVector, Vector2, U2, U6, U7};
use crate::od::msr::{RangeMsr, RangeRate};
use crate::od::Measurement;
use crate::TimeTagged;
use arrow::datatypes::{DataType, Field};

/// A simultaneous range and Doppler measurement in units of km and km/s.
/// This is a two-way measurement only.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RangeDoppler {
    /// Epoch of the observation
    pub epoch: Epoch,
    /// Observation vector in km and km/s
    pub obs: Vector2<f64>,
}

// TODO: Init using the state vectors of both the transmitter and receiver at two different epochs. If epochs are identical, treat it as an instantaneous measurement.
// Also allow for an initialization that only has one epoch for instantaneous measurements.
// The sensitivity should be computed on request only and using the nominal receiver state.

impl RangeDoppler {
    /// Initialize a new one-way range and Doppler measurement from the provided states and the effective noises.
    ///
    /// # Panics
    /// + If the epochs of the two states differ.
    /// + If the frames of the two states differ.
    pub fn one_way(
        tx: Orbit,
        rx: Orbit,
        clock_noise_s: f64,
        range_noise_km: f64,
        doppler_noise_km_s: f64,
    ) -> Self {
        assert_eq!(tx.frame, rx.frame, "tx & rx in different frames");
        assert_eq!(tx.epoch, rx.epoch, "tx & rx states have different times");

        let range_vec_km = tx.radius() - rx.radius();
        let doppler_km_s = range_vec_km.dot(&(tx.velocity() - rx.velocity())) / range_vec_km.norm();

        Self {
            epoch: tx.epoch + clock_noise_s * Unit::Second,
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
        clock_noise_s: f64,
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
        let delta_t_s = (tx.1.epoch - rx.1.epoch).to_seconds();

        // Compute the one way range and range rate measurements without any noise.
        let one_way_t_1 = Self::one_way(tx.0, rx.0, 0.0, 0.0, 0.0);
        let one_way_t_2 = Self::one_way(tx.1, rx.1, 0.0, 0.0, 0.0);

        let obs = (one_way_t_2.obs - one_way_t_1.obs) / delta_t_s;

        // Add the noise.
        Self {
            epoch: tx.0.epoch + clock_noise_s * Unit::Second,
            obs: Vector2::new(obs[0] + range_noise_km, obs[1] + doppler_noise_km_s),
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

// impl SimMeasurement for RangeDoppler {
//     type State = Orbit;

//     fn sensitivity(&self, _nominal: Orbit) -> Matrix2x6<f64> {
//         todo!()
//     }
// }
