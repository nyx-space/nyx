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

use super::{MvnSpacecraft, StateDispersion};
use crate::md::StateParameter;
use crate::Spacecraft;
use anise::analysis::prelude::OrbitalElement;
use nalgebra::{SMatrix, SVector};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyType;

fn parse_state_parameter(name: &str) -> PyResult<StateParameter> {
    match name.to_lowercase().as_str() {
        "x" => Ok(StateParameter::Element(OrbitalElement::X)),
        "y" => Ok(StateParameter::Element(OrbitalElement::Y)),
        "z" => Ok(StateParameter::Element(OrbitalElement::Z)),
        "vx" => Ok(StateParameter::Element(OrbitalElement::VX)),
        "vy" => Ok(StateParameter::Element(OrbitalElement::VY)),
        "vz" => Ok(StateParameter::Element(OrbitalElement::VZ)),
        "sma" | "semimajoraxis" => Ok(StateParameter::Element(OrbitalElement::SemiMajorAxis)),
        "ecc" | "eccentricity" => Ok(StateParameter::Element(OrbitalElement::Eccentricity)),
        "inc" | "inclination" => Ok(StateParameter::Element(OrbitalElement::Inclination)),
        "raan" => Ok(StateParameter::Element(OrbitalElement::RAAN)),
        "aop" | "argumentofperiapsis" => Ok(StateParameter::Element(OrbitalElement::AoP)),
        "ta" | "trueanomaly" => Ok(StateParameter::Element(OrbitalElement::TrueAnomaly)),
        "rmag" => Ok(StateParameter::Element(OrbitalElement::Rmag)),
        "vmag" => Ok(StateParameter::Element(OrbitalElement::Vmag)),
        "cr" => Ok(StateParameter::Cr),
        "cd" => Ok(StateParameter::Cd),
        "drymass" => Ok(StateParameter::DryMass),
        "propmass" => Ok(StateParameter::PropMass),
        "totalmass" => Ok(StateParameter::TotalMass),
        "isp" => Ok(StateParameter::Isp),
        "thrust" => Ok(StateParameter::Thrust),
        "bdotr" => Ok(StateParameter::BdotR),
        "bdott" => Ok(StateParameter::BdotT),
        "bltof" => Ok(StateParameter::BLTOF),
        "epoch" => Ok(StateParameter::Epoch),
        "guidancemode" => Ok(StateParameter::GuidanceMode),
        _ => Err(PyValueError::new_err(format!("Unknown state parameter: {}", name))),
    }
}

#[pymethods]
impl MvnSpacecraft {
    #[new]
    fn py_new(template: Spacecraft, dispersions: Vec<StateDispersion>) -> PyResult<Self> {
        MvnSpacecraft::new(template, dispersions).map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "from_spacecraft_cov")]
    fn py_from_spacecraft_cov(
        _cls: &Bound<'_, PyType>,
        template: Spacecraft,
        cov: Vec<Vec<f64>>,
        mean: Vec<f64>,
    ) -> PyResult<Self> {
        if cov.len() != 9 || cov.iter().any(|row| row.len() != 9) {
            return Err(PyValueError::new_err(
                "Covariance matrix must be 9x9 (rows and columns)",
            ));
        }
        if mean.len() != 9 {
            return Err(PyValueError::new_err("Mean vector must be length 9"));
        }

        let mut cov_mat = SMatrix::<f64, 9, 9>::zeros();
        for (i, row) in cov.iter().enumerate() {
            for (j, val) in row.iter().enumerate() {
                cov_mat[(i, j)] = *val;
            }
        }

        let mut mean_vec = SVector::<f64, 9>::zeros();
        for (i, val) in mean.iter().enumerate() {
            mean_vec[i] = *val;
        }

        MvnSpacecraft::from_spacecraft_cov(template, cov_mat, mean_vec)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "zero_mean")]
    fn py_zero_mean(
        _cls: &Bound<'_, PyType>,
        template: Spacecraft,
        dispersions: Vec<StateDispersion>,
    ) -> PyResult<Self> {
        MvnSpacecraft::zero_mean(template, dispersions)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }
}

#[pymethods]
impl StateDispersion {
    #[new]
    fn py_new(param: String, std_dev: Option<f64>, mean: Option<f64>) -> PyResult<Self> {
        let param_enum = parse_state_parameter(&param)?;
        let builder = StateDispersion::builder().param(param_enum);
        Ok(match (std_dev, mean) {
            (Some(s), Some(m)) => builder.std_dev(s).mean(m).build(),
            (Some(s), None) => builder.std_dev(s).build(),
            (None, Some(m)) => builder.mean(m).build(),
            (None, None) => builder.build(),
        })
    }

    #[classmethod]
    #[pyo3(name = "zero_mean")]
    fn py_zero_mean(_cls: &Bound<'_, PyType>, param: String, std_dev: f64) -> PyResult<Self> {
        let param_enum = parse_state_parameter(&param)?;
        Ok(StateDispersion::zero_mean(param_enum, std_dev))
    }
}
