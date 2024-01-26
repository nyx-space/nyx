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

use std::collections::BTreeMap;

use crate::cosmic::{Cosm, Orbit};
use crate::io::{ConfigRepr, ParseSnafu};
use crate::od::simulator::TrackingDeviceSim;
pub use crate::od::simulator::TrkConfig;
use crate::od::ODError;
use crate::python::PythonError;
use crate::time::Duration;
use crate::NyxError;
pub use crate::{io::ConfigError, od::noise::GaussMarkov, od::prelude::GroundStation};

use crate::python::cosmic::Cosm as CosmPy;
use crate::python::pyo3utils::pyany_to_value;

use pyo3::class::basic::CompareOp;
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList, PyType};
use pythonize::pythonize;
use snafu::ResultExt;

#[pymethods]
impl GroundStation {
    #[new]
    fn new(
        name: String,
        elevation_mask_deg: f64,
        latitude_deg: f64,
        longitude_deg: f64,
        height_km: f64,
        frame: String,
        light_time_correction: bool,
        integration_time: Option<Duration>,
        timestamp_noise_s: Option<GaussMarkov>,
        range_noise_km: Option<GaussMarkov>,
        doppler_noise_km_s: Option<GaussMarkov>,
    ) -> Result<Self, NyxError> {
        let cosm = Cosm::de438();
        let frame_obj = cosm.try_frame(&frame)?;
        Ok(Self {
            name,
            elevation_mask_deg,
            latitude_deg,
            longitude_deg,
            height_km,
            frame: frame_obj,
            integration_time,
            light_time_correction,
            timestamp_noise_s,
            range_noise_km,
            doppler_noise_km_s,
        })
    }

    fn __getnewargs__(
        &self,
    ) -> Result<
        (
            String,
            f64,
            f64,
            f64,
            f64,
            String,
            bool,
            Option<Duration>,
            Option<GaussMarkov>,
            Option<GaussMarkov>,
            Option<GaussMarkov>,
        ),
        NyxError,
    > {
        Ok((
            self.name.clone(),
            self.elevation_mask_deg,
            self.latitude_deg,
            self.longitude_deg,
            self.height_km,
            self.frame.to_string(),
            self.light_time_correction,
            self.integration_time,
            self.timestamp_noise_s,
            self.range_noise_km,
            self.doppler_noise_km_s,
        ))
    }

    #[classmethod]
    fn load(_cls: &PyType, path: &str) -> Result<Self, ConfigError> {
        <Self as ConfigRepr>::load(path)
    }

    #[classmethod]
    fn load_many(_cls: &PyType, path: &str) -> Result<Vec<Self>, ConfigError> {
        <Self as ConfigRepr>::load_many(path)
    }

    #[classmethod]
    fn load_named(_cls: &PyType, path: &str) -> Result<BTreeMap<String, Self>, ConfigError> {
        <Self as ConfigRepr>::load_named(path)
    }

    #[classmethod]
    /// Loads the SpacecraftDynamics from its YAML representation
    fn loads(_cls: &PyType, data: &PyAny) -> Result<BTreeMap<String, Self>, ConfigError> {
        if let Ok(as_list) = data.downcast::<PyList>() {
            let mut as_map = BTreeMap::new();
            for item in as_list.iter() {
                // Check that the item is a dictionary
                let next: Self =
                    serde_yaml::from_value(pyany_to_value(item)?).with_context(|_| ParseSnafu)?;
                as_map.insert(next.name.clone(), next);
            }
            Ok(as_map)
        } else if let Ok(as_dict) = data.downcast::<PyDict>() {
            let mut as_map = BTreeMap::new();
            for item_as_list in as_dict.items() {
                let k_any = item_as_list.get_item(0).or_else(|_| {
                    Err(ConfigError::InvalidConfig {
                        msg: "could not get key from provided dictionary item".to_string(),
                    })
                })?;
                let v_any = item_as_list.get_item(1).or_else(|_| {
                    Err(ConfigError::InvalidConfig {
                        msg: "could not get key from provided dictionary item".to_string(),
                    })
                })?;
                // Check that it's a string, or abort here
                if let Ok(as_str) = k_any.extract::<String>() {
                    // Try to convert the underlying data
                    match pyany_to_value(v_any) {
                        Ok(value) => {
                            let next: Self =
                                serde_yaml::from_value(value).with_context(|_| ParseSnafu)?;
                            as_map.insert(as_str, next);
                        }
                        Err(_) => {
                            // Maybe this was to be parsed in full as a single item
                            let me: Self = serde_yaml::from_value(pyany_to_value(as_dict)?)
                                .with_context(|_| ParseSnafu)?;
                            as_map.clear();
                            as_map.insert(me.name.clone(), me);
                            return Ok(as_map);
                        }
                    }
                } else {
                    return Err(ConfigError::InvalidConfig {
                        msg:
                            "keys for `loads` must be strings, otherwise, use GroundStation(**data)"
                                .to_string(),
                    });
                }
            }
            Ok(as_map)
        } else {
            Err(ConfigError::InvalidConfig {
                msg: "config must be dict, list, or str".to_string(),
            })
        }
    }

    fn dumps(&self, py: Python) -> Result<PyObject, NyxError> {
        pythonize(py, &self).map_err(|e| NyxError::CustomError { msg: e.to_string() })
    }

    /// Perform a one-way measurement of the given orbit at the epoch stored in that orbit instance.
    /// Returns the range in kilometers and the Doppler measurement in kilometers per second.
    fn measure(&mut self, orbit: Orbit) -> Result<Option<(f64, f64)>, ODError> {
        match self.measure_instantaneous(orbit, None, Cosm::de438())? {
            Some(msr) => Ok(Some((msr.obs[0], msr.obs[1]))),
            None => Ok(None),
        }
    }

    /// Computes the azimuth and elevation of the provided object seen from this ground station, both in degrees.
    fn compute_azimuth_elevation(&self, receiver: Orbit, cosm: &CosmPy) -> (f64, f64) {
        let (az_deg, el_deg, _, _) = self.azimuth_elevation_of(receiver, &cosm.inner);

        (az_deg, el_deg)
    }

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

    fn __richcmp__(&self, other: &Self, op: CompareOp) -> Result<bool, PythonError> {
        match op {
            CompareOp::Eq => Ok(self == other),
            CompareOp::Ne => Ok(self != other),
            _ => Err(PythonError::OperationError { op }),
        }
    }
}
