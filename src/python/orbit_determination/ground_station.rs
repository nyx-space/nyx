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

use crate::io::ConfigRepr;
pub use crate::od::simulator::TrkConfig;
pub use crate::{io::ConfigError, od::prelude::GroundStation};

use pyo3::prelude::*;
use pyo3::types::PyType;

#[pymethods]
impl GroundStation {
    #[cfg(feature = "python")]
    #[classmethod]
    fn load(_cls: &PyType, path: &str) -> Result<Self, ConfigError> {
        <Self as ConfigRepr>::load(path)
    }

    #[classmethod]
    fn load_many(_cls: &PyType, path: &str) -> Result<Vec<Self>, ConfigError> {
        <Self as ConfigRepr>::load_many(path)
    }

    #[classmethod]
    fn load_named(_cls: &PyType, path: &str) -> Result<HashMap<String, Self>, ConfigError> {
        <Self as ConfigRepr>::load_named(path)
    }

    // TODO: Once the switch to RangeDoppler is done, add this method. This requires making the other structure available to Python too
    // and that isn't worth doing for StdMeasurement.
    // fn measure(&mut self, epoch: Epoch, state: &State) -> Result<Measurement, NyxError> {
    //     self.measure(epoch, state)
    // }

    // Manual getter/setters -- waiting on https://github.com/PyO3/pyo3/pull/2786

    #[getter]
    fn get_name(&self) -> PyResult<String> {
        Ok(self.name.clone())
    }

    #[setter]
    fn set_name(&mut self, name: String) -> PyResult<()> {
        self.name = name;
        Ok(())
    }

    #[getter]
    fn get_elevation_mask_deg(&self) -> PyResult<f64> {
        Ok(self.elevation_mask_deg)
    }

    #[setter]
    fn set_elevation_mask_deg(&mut self, mask_deg: f64) -> PyResult<()> {
        self.elevation_mask_deg = mask_deg;
        Ok(())
    }

    #[getter]
    fn get_latitude_deg(&self) -> PyResult<f64> {
        Ok(self.latitude_deg)
    }

    #[setter]
    fn set_latitude_deg(&mut self, lat_deg: f64) -> PyResult<()> {
        self.latitude_deg = lat_deg;
        Ok(())
    }

    #[getter]
    fn get_longitude_deg(&self) -> PyResult<f64> {
        Ok(self.longitude_deg)
    }

    #[setter]
    fn set_longitude_deg(&mut self, long_deg: f64) -> PyResult<()> {
        self.longitude_deg = long_deg;
        Ok(())
    }

    #[getter]
    fn get_height_km(&self) -> PyResult<f64> {
        Ok(self.height_km)
    }

    #[setter]
    fn set_height_km(&mut self, height_km: f64) -> PyResult<()> {
        self.height_km = height_km;
        Ok(())
    }

    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self}")
    }
}
