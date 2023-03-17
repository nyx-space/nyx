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

use hifitime::{Epoch, Unit};
use pyo3::prelude::*;

use crate::{
    md::{ui::Traj as TrajRs, Event},
    NyxError, Spacecraft, State,
};

/// A structure that stores a spacecraft structure generated from a propagation
#[pyclass]
pub(crate) struct Traj {
    pub(crate) inner: TrajRs<Spacecraft>,
}

#[pymethods]
impl Traj {
    /// Returns the state at the provided epoch, or raises an exception if the epoch is outside of the bounds of the trajectory
    fn at(&self, epoch: Epoch) -> Result<Spacecraft, NyxError> {
        self.inner.at(epoch)
    }

    /// Return the first state of the trajectory
    fn first(&self) -> Spacecraft {
        *self.inner.first()
    }

    /// Return the last state of the trajectory
    fn last(&self) -> Spacecraft {
        *self.inner.last()
    }

    /// Finds a specific event in a trajectory.
    ///
    /// If a start or end epoch is provided (or both are provided), this function will return a list of a single event.
    /// If none are provided, this function will search the whole trajectory for the event and return all of the states where such event happens.
    fn find(
        &self,
        event: Event,
        start: Option<Epoch>,
        end: Option<Epoch>,
    ) -> Result<Vec<Spacecraft>, NyxError> {
        if start.is_some() || end.is_some() {
            let start = if let Some(start) = start {
                start
            } else {
                self.inner.first().epoch()
            };

            let end = if let Some(end) = end {
                end
            } else {
                self.inner.last().epoch()
            };

            Ok(vec![self.inner.find_bracketed(start, end, &event)?])
        } else {
            self.inner.find_all(&event)
        }
    }

    /// Find the minimum and maximum of the provided event through the trajectory with a specified time unit precision.
    pub fn find_minmax(
        &self,
        event: Event,
        precision: Unit,
    ) -> Result<(Spacecraft, Spacecraft), NyxError> {
        self.inner.find_minmax(&event, precision)
    }

    /// Saves this trajectory to a parquet file.
    fn to_parquet(&self, path: String) -> Result<String, NyxError> {
        self.inner
            .to_parquet(path, None)
            .map_err(|e| NyxError::CustomError(e.to_string()))
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
}
