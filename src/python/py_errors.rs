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

use pyo3::class::basic::CompareOp;
use pyo3::{exceptions::PyException, prelude::*};

#[derive(Snafu, Debug)]
pub(crate) enum PythonError {
    #[snafu(display("operation {op:?} not available on this type"))]
    OperationError { op: CompareOp },
}

impl From<PythonError> for PyErr {
    fn from(err: PythonError) -> PyErr {
        PyException::new_err(err.to_string())
    }
}
