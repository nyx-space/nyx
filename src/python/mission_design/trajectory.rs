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

use hifitime::Epoch;
use pyo3::prelude::*;

use crate::{md::ui::Traj as TrajRs, NyxError, Spacecraft};

/// A structure that stores a spacecraft structure generated from a propagation
#[pyclass]
pub(crate) struct Traj {
    pub(crate) inner: TrajRs<Spacecraft>,
}

#[pymethods]
impl Traj {
    fn at(&self, epoch: Epoch) -> Result<Spacecraft, NyxError> {
        self.inner.at(epoch)
    }

    fn to_parquet(&self, path: String) -> Result<String, NyxError> {
        self.inner
            .to_parquet(path, None)
            .map_err(|e| NyxError::CustomError(e.to_string()))
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
}
