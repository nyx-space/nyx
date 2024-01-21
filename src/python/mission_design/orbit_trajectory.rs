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

use hifitime::{Duration, Epoch, Unit};
use pyo3::prelude::*;

use super::SpacecraftTraj;
use crate::md::trajectory::ExportCfg;
use crate::python::cosmic::Frame;
use crate::{
    md::{
        prelude::{Cosm, Traj as TrajRs},
        Event, EventEvaluator,
    },
    NyxError, Orbit, Spacecraft, State,
};

use std::collections::HashMap;

/// An orbit trajectory. Cannot be pickled in Python, so you must export it to Parquet first and use the TrajectoryLoader.
#[pyclass]
pub(crate) struct OrbitTraj {
    pub(crate) inner: TrajRs<Orbit>,
}

#[pymethods]
impl OrbitTraj {
    /// Convert this orbit trajectory into a spacecraft trajectory by copying the provided template and setting its orbit state to that of each state of the trajectory
    pub fn upcast(&self, template: Spacecraft) -> SpacecraftTraj {
        SpacecraftTraj {
            inner: self.inner.upcast(template),
        }
    }

    /// Returns the state at the provided epoch, or raises an exception if the epoch is outside of the bounds of the trajectory
    fn at(&self, epoch: Epoch) -> Result<Orbit, NyxError> {
        self.inner.at(epoch)
    }

    /// Return the first state of the trajectory
    fn first(&self) -> Orbit {
        *self.inner.first()
    }

    /// Return the last state of the trajectory
    fn last(&self) -> Orbit {
        *self.inner.last()
    }

    /// Copies this object and resamples it with the provided step size
    fn resample(&self, step: Duration) -> Result<Self, NyxError> {
        let inner = self.inner.resample(step)?;

        Ok(Self { inner })
    }

    /// Copies this object and rebuilds it with the provided epochs
    fn rebuild(&self, epochs: Vec<Epoch>) -> Result<Self, NyxError> {
        let inner = self.inner.rebuild(&epochs)?;

        Ok(Self { inner })
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
    ) -> Result<Vec<Orbit>, NyxError> {
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

            Ok(vec![self.inner.find_bracketed(start, end, &event)?.state])
        } else {
            Ok(self
                .inner
                .find(&event)?
                .iter()
                .map(|details| details.state)
                .collect::<Vec<Orbit>>())
        }
    }

    /// Find the minimum and maximum of the provided event through the trajectory with a specified time unit precision.
    pub fn find_minmax(&self, event: Event, precision: Unit) -> Result<(Orbit, Orbit), NyxError> {
        self.inner.find_minmax(&event, precision)
    }

    /// Saves this trajectory to a parquet file, optionally adding the event columns to append and metadata.
    /// Set the groundtrack parameter to a body fixed frame to export this trajectory with latitude, longitude, and height columns in that body fixed frame.
    #[pyo3(text_signature = "(path, events=None, metadata=None, groundtrack=None")]
    fn to_parquet(
        &self,
        path: String,
        events: Option<Vec<Event>>,
        metadata: Option<HashMap<String, String>>,
        groundtrack: Option<Frame>,
    ) -> Result<String, NyxError> {
        let maybe = match events {
            None => match groundtrack {
                None => {
                    let cfg = match metadata {
                        None => ExportCfg::default(),
                        Some(_) => ExportCfg {
                            metadata,
                            ..Default::default()
                        },
                    };
                    self.inner.to_parquet(path, None, cfg)
                }
                Some(py_frame) => self.inner.to_groundtrack_parquet(
                    path,
                    py_frame.inner,
                    None,
                    metadata,
                    Cosm::de438(),
                ),
            },
            Some(events) => {
                let events = events
                    .iter()
                    .map(|e| e as &dyn EventEvaluator<Orbit>)
                    .collect::<Vec<&dyn EventEvaluator<Orbit>>>();

                match groundtrack {
                    None => {
                        let cfg = match metadata {
                            None => ExportCfg::default(),
                            Some(_) => ExportCfg {
                                metadata,
                                ..Default::default()
                            },
                        };
                        self.inner.to_parquet(path, Some(events), cfg)
                    }
                    Some(py_frame) => self.inner.to_groundtrack_parquet(
                        path,
                        py_frame.inner,
                        Some(events),
                        metadata,
                        Cosm::de438(),
                    ),
                }
            }
        };

        match maybe {
            Ok(path) => Ok(format!("{}", path.to_str().unwrap())),
            Err(e) => Err(NyxError::CustomError { msg: e.to_string() }),
        }
    }

    /// Allows converting the source trajectory into the (almost) equivalent trajectory in another frame.
    /// This simply converts each state into the other frame and may lead to aliasing due to the Nyquist–Shannon sampling theorem.
    fn to_frame(&self, new_frame: String) -> Result<Self, NyxError> {
        let cosm = Cosm::de438();

        let frame = cosm.try_frame(&new_frame)?;

        let conv_traj = self.inner.to_frame(frame, cosm)?;

        Ok(Self { inner: conv_traj })
    }

    /// Compute the RIC difference between this trajectory and another, writing the output to a parquet file.
    fn ric_diff_to_parquet(
        &self,
        other: &Self,
        path: String,
        cfg: Option<ExportCfg>,
    ) -> Result<String, NyxError> {
        match self.inner.ric_diff_to_parquet(
            &other.inner,
            path,
            cfg.unwrap_or_else(|| ExportCfg::default()),
        ) {
            Ok(path) => Ok(format!("{}", path.to_str().unwrap())),
            Err(e) => Err(NyxError::CustomError { msg: e.to_string() }),
        }
    }

    fn __add__(&self, rhs: &Self) -> Result<Self, NyxError> {
        let inner = (self.inner.clone() + rhs.inner.clone())?;

        Ok(Self { inner })
    }

    fn __str__(&self) -> String {
        format!("{}", self.inner)
    }
}
