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
use anise::ephemerides::ephemeris::Ephemeris;
use hifitime::Epoch;
use log::info;
use nyx_space::{
    Spacecraft,
    io::InputOutputError,
    md::{
        Trajectory,
        trajectory::{ExportCfg, TrajError},
    },
};

use pyo3::prelude::*;

#[pyclass(from_py_object, name = "Trajectory")]
#[derive(Clone, Debug)]
pub struct PyTrajectory {
    pub inner: Trajectory,
}

#[pymethods]
impl PyTrajectory {
    #[new]
    fn py_new(path: String, template: Option<Spacecraft>) -> Result<Self, InputOutputError> {
        if path.to_lowercase().ends_with(".e") {
            info!("loading {path} as STK .e");
            let ephem = Ephemeris::from_stk_e_file(path).map_err(|e| {
                InputOutputError::UnsupportedData {
                    which: e.to_string(),
                }
            })?;
            // Rebuild a trajectory by applying the template
            let template = template.unwrap_or_default();
            let mut inner = Trajectory::default();
            for record in &ephem {
                inner.states.push(template.with_orbit(record.orbit));
            }
            inner.name = Some(ephem.object_id);

            Ok(Self { inner })
        } else if path.to_lowercase().ends_with(".oem") {
            info!("loading {path} as CCSDS OEM");
            let inner = Trajectory::from_oem_file(path, template).map_err(|e| {
                InputOutputError::UnsupportedData {
                    which: e.to_string(),
                }
            })?;
            Ok(Self { inner })
        } else {
            info!("loading {path} as Nyx Parquet");
            let inner = Trajectory::from_parquet(path)?;
            Ok(Self { inner })
        }
    }

    /// Export the difference in RIC from of this trajectory compare to the "other" trajectory in parquet format.
    ///
    /// # Notes
    /// + The RIC frame accounts for the transport theorem by performing a finite differencing of the RIC frame.
    fn ric_diff_to_parquet(
        &self,
        other: &Self,
        path: &str,
        cfg: ExportCfg,
    ) -> Result<String, TrajError> {
        self.inner
            .ric_diff_to_parquet(&other.inner, path, cfg)
            .map(|path| path.to_string_lossy().to_string())
    }

    /// Evaluate the trajectory at this specific epoch.
    fn at(&self, epoch: Epoch) -> Result<Spacecraft, TrajError> {
        self.inner.at(epoch)
    }

    fn to_parquet(&self, path: &str, cfg: ExportCfg) -> Result<String, TrajError> {
        self.inner
            .to_parquet(path, cfg)
            .map(|path| path.to_string_lossy().to_string())
            .map_err(|e| TrajError::TrajGeneric { err: e.to_string() })
    }
    /// Export this spacecraft trajectory estimate to an ANISE Ephemeris
    fn to_ephemeris(&self, object_id: String, cfg: ExportCfg) -> Ephemeris {
        self.inner.to_ephemeris(object_id, cfg)
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
}
