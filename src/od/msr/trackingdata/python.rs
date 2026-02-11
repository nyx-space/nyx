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

use super::{MeasurementType, TrackingDataArc};
use crate::io::{ExportCfg, InputOutputError};
use pyo3::prelude::*;
use pyo3::types::PyType;
use std::collections::HashMap;

#[pymethods]
impl TrackingDataArc {
    /// Initializes a new Almanac from a file path to CCSDS OEM file, after converting to to SPICE SPK/BSP
    ///
    /// :type path: str
    /// :type aliases: dict
    /// :rtype: nyx_space.od.TrackingDataArc
    #[classmethod]
    #[pyo3(name = "from_ccsds_tdm")]
    fn py_from_ccsds_tdm_file(
        _cls: Bound<'_, PyType>,
        path: &str,
        aliases: Option<HashMap<String, String>>,
    ) -> Result<Self, InputOutputError> {
        TrackingDataArc::from_tdm(path, aliases)
    }

    #[pyo3(name = "write_ccsds_tdm")]
    fn py_write_ccsds_tdm(
        &self,
        spacecraft_name: String,
        aliases: Option<HashMap<String, String>>,
        path: &str,
    ) -> Result<String, InputOutputError> {
        Ok(self
            .clone()
            .to_tdm_file(spacecraft_name, aliases, path, ExportCfg::default())?
            .to_str()
            .unwrap_or("woah_bug_building_path")
            .to_string())
    }

    #[pyo3(name = "unique_aliases")]
    fn py_unique_aliases(&self) -> Vec<String> {
        self.unique_aliases().iter().cloned().collect()
    }
    #[pyo3(name = "unique_types")]
    fn py_unique_types(&self) -> Vec<MeasurementType> {
        self.unique_types().iter().cloned().collect()
    }

    fn __str__(&self) -> String {
        format!("{self}")
    }
}
