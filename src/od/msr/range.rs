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
use crate::linalg::{DimName, Matrix1x6, OVector, Vector1, U1, U6, U7};
use crate::od::{Measurement, SimMeasurement};
use crate::time::Epoch;
use crate::TimeTagged;
use arrow::datatypes::{DataType, Field};
use hyperdual::linalg::norm;
use hyperdual::{hyperspace_from_vector, OHyperdual};
use serde::ser::SerializeSeq;
use serde::{Serialize, Serializer};
use std::collections::HashMap;

/// Stores a standard measurement of range (km)
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RangeMsr {
    pub dt: Epoch,
    pub obs: Vector1<f64>,
    visible: bool,
    h_tilde: Matrix1x6<f64>,
}

impl RangeMsr {
    pub fn range(&self) -> f64 {
        self.obs[(0, 0)]
    }

    fn compute_sensitivity(
        state: &OVector<OHyperdual<f64, U7>, U6>,
    ) -> (Vector1<f64>, Matrix1x6<f64>) {
        // Extract data from hyperspace
        let range_vec = state.fixed_rows::<3>(0).into_owned();

        // Code up math as usual
        let range = norm(&range_vec);

        // Extract result into Vector2 and Matrix2x6
        let fx = Vector1::new(range.real());
        let mut pmat = Matrix1x6::zeros();

        for j in 1..U7::dim() {
            pmat[(j - 1)] = range[j];
        }

        (fx, pmat)
    }

    pub fn new(tx: Orbit, rx: Orbit, visible: bool) -> RangeMsr {
        assert_eq!(tx.frame, rx.frame, "tx and rx in different frames");
        assert_eq!(tx.epoch, rx.epoch, "tx and rx states have different times");

        let dt = tx.epoch;
        let hyperstate = hyperspace_from_vector(&(rx - tx).to_cartesian_vec());
        let (obs, h_tilde) = Self::compute_sensitivity(&hyperstate);

        RangeMsr {
            dt,
            obs,
            visible,
            h_tilde,
        }
    }
}

impl Measurement for RangeMsr {
    type MeasurementSize = U1;

    /// Returns this measurement as a vector of Range and Range Rate
    ///
    /// **Units:** km
    fn observation(&self) -> Vector1<f64> {
        self.obs
    }

    fn fields() -> Vec<Field> {
        let mut meta = HashMap::new();
        meta.insert("unit".to_string(), "km".to_string());
        vec![Field::new("Range (km)", DataType::Float64, false).with_metadata(meta)]
    }

    fn from_observation(epoch: Epoch, obs: OVector<f64, Self::MeasurementSize>) -> Self {
        Self {
            dt: epoch,
            obs,
            visible: true,
            h_tilde: Matrix1x6::zeros(),
        }
    }
}

impl SimMeasurement for RangeMsr {
    type State = Orbit;

    fn sensitivity(&self, _nominal: Orbit) -> Matrix1x6<f64> {
        self.h_tilde
    }
}

impl TimeTagged for RangeMsr {
    fn epoch(&self) -> Epoch {
        self.dt
    }

    fn set_epoch(&mut self, dt: Epoch) {
        self.dt = dt
    }
}

impl Serialize for RangeMsr {
    /// Serializes the observation vector at the given time.
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        let mut seq = serializer.serialize_seq(Some(3))?;
        seq.serialize_element(&self.dt.to_mjd_tai_days())?;
        let obs = self.observation();
        seq.serialize_element(&obs[(0, 0)])?;
        seq.end()
    }
}
