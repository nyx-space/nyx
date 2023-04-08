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
use crate::linalg::{DimName, Matrix2x6, OVector, Vector2, U2, U6, U7};
use crate::od::Measurement;
use crate::time::Epoch;
use crate::TimeTagged;
use arrow::datatypes::Field;
use hyperdual::linalg::norm;
use hyperdual::{hyperspace_from_vector, OHyperdual};
use serde::ser::SerializeSeq;
use serde::{Serialize, Serializer};

use super::{RangeMsr, RangeRate};

/// Stores a standard measurement of range (km) and range rate (km/s)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StdMeasurement {
    pub epoch: Epoch,
    pub obs: Vector2<f64>,
}

impl StdMeasurement {
    pub fn range_km(&self) -> f64 {
        self.obs[(0, 0)]
    }
    pub fn doppler_km_s(&self) -> f64 {
        self.obs[(1, 0)]
    }

    fn compute_sensitivity(
        state: &OVector<OHyperdual<f64, U7>, U6>,
        range_noise: f64,
        range_rate_noise: f64,
    ) -> (Vector2<f64>, Matrix2x6<f64>) {
        // Extract data from hyperspace
        let range_vec = state.fixed_rows::<3>(0).into_owned();
        let velocity_vec = state.fixed_rows::<3>(3).into_owned();

        // Code up math as usual
        let delta_v_vec = velocity_vec / norm(&range_vec);
        let range = norm(&range_vec) + OHyperdual::from(range_noise);
        let range_rate = range_vec.dot(&delta_v_vec) + OHyperdual::from(range_rate_noise);

        // Extract result into Vector2 and Matrix2x6
        let mut fx = Vector2::zeros();
        let mut pmat = Matrix2x6::zeros();
        for i in 0..U2::dim() {
            fx[i] = if i == 0 {
                range.real()
            } else {
                range_rate.real()
            };
            for j in 1..U7::dim() {
                pmat[(i, j - 1)] = if i == 0 { range[j] } else { range_rate[j] };
            }
        }

        (fx, pmat)
    }

    /// Generate noiseless measurement
    pub fn noiseless(dt: Epoch, tx: Orbit, rx: Orbit) -> StdMeasurement {
        Self::raw(dt, tx, rx, 0.0, 0.0)
    }

    /// Generate a new measurement with the provided noise distribution.
    pub fn new(
        epoch: Epoch,
        tx: Orbit,
        rx: Orbit,
        range_noise: f64,
        doppler_noise: f64,
    ) -> StdMeasurement {
        Self::raw(epoch, tx, rx, range_noise, doppler_noise)
    }

    /// Generate a new measurement with the provided noise values.
    pub fn raw(
        epoch: Epoch,
        tx: Orbit,
        rx: Orbit,
        range_noise: f64,
        doppler_noise: f64,
    ) -> StdMeasurement {
        assert_eq!(tx.frame, rx.frame, "tx & rx in different frames");
        assert_eq!(tx.epoch, rx.epoch, "tx & rx states have different times");

        let hyperstate = hyperspace_from_vector(&(rx - tx).to_cartesian_vec());
        let (obs, _h_tilde) = Self::compute_sensitivity(&hyperstate, range_noise, doppler_noise);

        StdMeasurement { epoch, obs }
    }

    pub fn real(dt: Epoch, range_km: f64, doppler_km_s: f64) -> Self {
        Self {
            epoch: dt,
            obs: Vector2::new(range_km, doppler_km_s),
        }
    }
}

impl Measurement for StdMeasurement {
    type MeasurementSize = U2;

    /// Returns this measurement as a vector of Range and Range Rate
    ///
    /// **Units:** km, km/s
    fn observation(&self) -> Vector2<f64> {
        self.obs
    }

    fn fields() -> Vec<Field> {
        vec![
            RangeMsr::fields()[0].clone(),
            RangeRate::fields()[0].clone(),
        ]
    }

    fn from_observation(epoch: Epoch, obs: OVector<f64, Self::MeasurementSize>) -> Self {
        Self { epoch, obs }
    }
}

impl TimeTagged for StdMeasurement {
    fn epoch(&self) -> Epoch {
        self.epoch
    }

    fn set_epoch(&mut self, dt: Epoch) {
        self.epoch = dt
    }
}

impl Serialize for StdMeasurement {
    /// Serializes the observation vector at the given time.
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(3))?;
        seq.serialize_element(&self.epoch.to_mjd_tai_days())?;
        let obs = self.observation();
        seq.serialize_element(&obs[(0, 0)])?;
        seq.serialize_element(&obs[(1, 0)])?;
        seq.end()
    }
}
