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

use crate::mc::GaussianGenerator;
use crate::md::StateParameter;
use crate::Orbit;
use crate::{NyxError, Spacecraft};
use pyo3::{prelude::*, py_run};
use rand::SeedableRng;
use rand_distr::Distribution;
use rand_pcg::Pcg64Mcg;

/// Monte Carlo
pub(crate) fn register_mc(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "_nyx_space.monte_carlo")?;

    sm.add_class::<StateParameter>()?;
    sm.add_function(wrap_pyfunction!(generate_orbits, sm)?)?;
    sm.add_function(wrap_pyfunction!(generate_spacecraft, sm)?)?;

    py_run!(
        py,
        sm,
        "import sys; sys.modules['nyx_space.monte_carlo'] = sm"
    );
    parent_module.add_submodule(sm)?;
    Ok(())
}

/// Generates orbits from the provided template orbit, the parameters to disperse, and whether these are absolute standard deviations or a percentage of the parameter's value.
#[pyfunction]
#[pyo3(text_signature = "(orbit, parameters, count, kind='abs', seed=None, /)")]
fn generate_orbits(
    orbit: Orbit,
    parameters: Vec<(StateParameter, f64)>,
    count: usize,
    kind: String,
    seed: Option<u64>,
) -> Result<Vec<Orbit>, NyxError> {
    let generator = match kind.as_str() {
        "abs" => GaussianGenerator::from_3std_devs(orbit, &parameters)?,
        "prct" => GaussianGenerator::from_3std_dev_prcts(orbit, &parameters)?,
        _ => {
            return Err(NyxError::CustomError {
                msg: format!(
                    "Unknown kind of distribution: {} (should be 'abs' or 'prct')",
                    kind
                ),
            })
        }
    };

    let rng = match seed {
        Some(seed) => Pcg64Mcg::new(seed.into()),
        None => Pcg64Mcg::from_entropy(),
    };

    // Generate multithreaded
    Ok(generator
        .sample_iter(rng)
        .take(count)
        .map(|disp_state| disp_state.state)
        .collect::<Vec<Orbit>>())
}

/// Generates spacecraft from the provided template spacecraft, the parameters to disperse, and whether these are absolute standard deviations or a percentage of the parameter's value.
#[pyfunction]
#[pyo3(text_signature = "(spacecraft, parameters, count, kind='abs', seed=None, /)")]
fn generate_spacecraft(
    spacecraft: Spacecraft,
    parameters: Vec<(StateParameter, f64)>,
    count: usize,
    kind: String,
    seed: Option<u64>,
) -> Result<Vec<Spacecraft>, NyxError> {
    let generator = match kind.as_str() {
        "abs" => GaussianGenerator::from_3std_devs(spacecraft, &parameters)?,
        "prct" => GaussianGenerator::from_3std_dev_prcts(spacecraft, &parameters)?,
        _ => {
            return Err(NyxError::CustomError {
                msg: format!(
                    "Unknown kind of distribution: {} (should be 'abs' or 'prct')",
                    kind
                ),
            })
        }
    };

    let rng = match seed {
        Some(seed) => Pcg64Mcg::new(seed.into()),
        None => Pcg64Mcg::from_entropy(),
    };

    // Generate multithreaded
    Ok(generator
        .sample_iter(rng)
        .take(count)
        .map(|disp_state| disp_state.state)
        .collect::<Vec<Spacecraft>>())
}
