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

use pyo3::pymethods;

use crate::od::msr::MeasurementType;

use super::Measurement;
use hifitime::Epoch;

#[pymethods]
impl Measurement {
    #[new]
    fn py_new(tracker: String, epoch: Epoch) -> Self {
        Self::new(tracker, epoch)
    }

    /// Returns the floating point value of this observation if this measurement contains the provided measurement type
    #[pyo3(name = "observation")]
    fn py_observation(&self, msr_type: MeasurementType) -> Option<f64> {
        self.data.get(&msr_type).copied()
    }

    fn __str__(&self) -> String {
        format!("{self}")
    }

    fn __repr__(&self) -> String {
        format!("{self} @ {self:p}")
    }
}
