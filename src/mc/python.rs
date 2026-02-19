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
use nalgebra::{SMatrix, SVector};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;
use pyo3::types::PyType;
use rand::SeedableRng;
use rand_distr::Distribution;
use rand_pcg::Pcg64Mcg;

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

        let cov_mat = SMatrix::<f64, 9, 9>::from_fn(|r, c| cov[r][c]);
        let mean_vec = SVector::<f64, 9>::from_vec(mean);

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

    /// Samples the multivariate distribution to generate a list of spacecraft states.
    ///
    /// The Pseudo-Random Number Generator (PRNG) used is the Permuted Congruential Generator (PCG).
    /// PCG is an excellent choice for Monte Carlo simulations because:
    /// 1. **Statistical Quality**: It passes difficult statistical tests (like TestU01), ensuring the generated numbers are random enough for high-fidelity simulations.
    /// 2. **Performance**: It is very fast and efficient, which is crucial when generating a large number of samples.
    /// 3. **Small State**: It has a small state size and is easy to seed, making it ideal for reproducible simulations.
    /// 4. **Reproducibility**: By providing a seed, the exact same sequence of spacecraft states can be generated, allowing for debugging and validation of Monte Carlo runs.
    #[pyo3(signature = (count, seed=None))]
    fn sample(&self, count: usize, seed: Option<u64>) -> Vec<Spacecraft> {
        let mut rng = match seed {
            Some(s) => Pcg64Mcg::seed_from_u64(s),
            None => Pcg64Mcg::from_rng(&mut rand::rng()),
        };

        let mut samples = Vec::with_capacity(count);
        for _ in 0..count {
            let dispersed_state = Distribution::sample(self, &mut rng);
            samples.push(dispersed_state.state);
        }
        samples
    }
}

#[pymethods]
impl StateDispersion {
    #[new]
    fn py_new(param: StateParameter, std_dev: Option<f64>, mean: Option<f64>) -> PyResult<Self> {
        let builder = StateDispersion::builder().param(param);
        Ok(match (std_dev, mean) {
            (Some(s), Some(m)) => builder.std_dev(s).mean(m).build(),
            (Some(s), None) => builder.std_dev(s).build(),
            (None, Some(m)) => builder.mean(m).build(),
            (None, None) => builder.build(),
        })
    }

    #[classmethod]
    #[pyo3(name = "zero_mean")]
    fn py_zero_mean(
        _cls: &Bound<'_, PyType>,
        param: StateParameter,
        std_dev: f64,
    ) -> PyResult<Self> {
        Ok(StateDispersion::zero_mean(param, std_dev))
    }
}
