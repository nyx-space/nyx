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
use anise::analysis::event::Event;
use anise::frames::Frame;
use anise::{ephemerides::ephemeris::Ephemeris, prelude::Almanac};
use hifitime::{Duration, Epoch};
use log::info;
use nyx_space::dynamics::SpacecraftDynamics;
use nyx_space::propagators::{IntegratorMethod, IntegratorOptions, Propagator as CorePropagator};
use nyx_space::{
    dynamics::sequence::Dynamics,
    io::InputOutputError,
    md::{
        trajectory::{ExportCfg, TrajError},
        Trajectory,
    },
    propagators::PropagationError,
    Spacecraft,
};
use rayon::prelude::*;
use std::sync::Arc;

use pyo3::prelude::*;

#[pyclass(unsendable, get_all, dict)]
pub struct PropagationResult {
    pub state: Spacecraft,
    pub trajectory: Option<PyTrajectory>,
}

impl PropagationResult {
    fn __str__(&self) -> String {
        format!(
            "{} (includes trajectory: {})",
            self.state,
            self.trajectory.is_some()
        )
    }

    fn __repr__(&self) -> String {
        format!("{self:p}")
    }
}

#[pyclass(from_py_object)]
#[derive(Clone)]
pub struct Propagator {
    #[pyo3(get, set)]
    pub dynamics: Dynamics,
    #[pyo3(get, set)]
    pub method: IntegratorMethod,
    #[pyo3(get, set)]
    pub options: IntegratorOptions,
    pub almanac: Arc<Almanac>,
}

impl Propagator {
    pub(crate) fn build(&self) -> Result<CorePropagator<SpacecraftDynamics>, PropagationError> {
        Ok(CorePropagator::new(
            self.dynamics
                .build(self.almanac.clone())
                .map_err(|msg| PropagationError::PropGenericError { msg })?,
            self.method,
            self.options,
        ))
    }
}

#[pymethods]
impl Propagator {
    #[pyo3(signature=(dynamics, almanac, method=IntegratorMethod::RungeKutta89, options=IntegratorOptions::default()))]
    #[new]
    fn py_new(
        dynamics: Dynamics,
        almanac: Almanac,
        method: IntegratorMethod,
        options: IntegratorOptions,
    ) -> Self {
        Self {
            dynamics,
            almanac: Arc::new(almanac),
            method,
            options,
        }
    }

    /// Propagates the initialization state until the desired epoch, optionally not building the trajectory
    #[pyo3(signature = (spacecraft, epoch, trajectory=true))]
    fn until_epoch(
        &self,
        spacecraft: Spacecraft,
        epoch: Epoch,
        trajectory: bool,
    ) -> Result<PropagationResult, PropagationError> {
        let setup = self.build()?;

        let mut prop = setup.with(spacecraft, self.almanac.clone());

        let (state, trajectory) = if trajectory {
            let (state, traj) = prop.until_epoch_with_traj(epoch)?;
            (state, Some(PyTrajectory { inner: traj }))
        } else {
            let state = prop.until_epoch(epoch)?;
            (state, None)
        };

        Ok(PropagationResult { state, trajectory })
    }

    /// Propagates the initialization state for the desired duration, optionally not building the trajectory
    #[pyo3(signature = (spacecraft, duration, trajectory=true))]
    fn for_duration(
        &self,
        spacecraft: Spacecraft,
        duration: Duration,
        trajectory: bool,
    ) -> Result<PropagationResult, PropagationError> {
        let setup = self.build()?;

        let mut prop = setup.with(spacecraft, self.almanac.clone());

        let (state, trajectory) = if trajectory {
            let (state, traj) = prop.for_duration_with_traj(duration)?;
            (state, Some(PyTrajectory { inner: traj }))
        } else {
            let state = prop.for_duration(duration)?;
            (state, None)
        };

        Ok(PropagationResult { state, trajectory })
    }

    /// Propagates the initialization state until the specified event has occurred `trigger` times, or until `max_duration` is reached.
    ///
    /// This method monitors the provided `event` during propagation. Once the event condition is met
    /// `trigger` number of times (e.g., set `trigger` to 1 for the first occurrence), the propagation stops
    /// at the end of that integration step.
    ///
    /// A root-finding algorithm (Brent's method) is then used to locate the exact time of the event
    /// within the final integration step. The returned state corresponds to this precise event time,
    /// interpolated from the trajectory.
    ///
    /// # Arguments
    ///
    /// * `max_duration` - The maximum duration to propagate if the event is not triggered the requested number of times.
    /// * `event` - The event definition (scalar expression and condition) to monitor.
    /// * `trigger` - The 1-based index of the event occurrence to stop at (e.g. 1 for the first crossing, 2 for the second).
    ///
    /// # Returns
    ///
    /// A tuple containing:
    /// 1. The interpolated state exactly at the moment the $n$-th event occurred.
    /// 2. The full trajectory recorded up to the end of the propagation step where the event occurred, unless explicitly ignored (but it is still built)
    ///
    /// # Errors
    ///
    /// * `PropagationError::NthEventError`: Returned if `max_duration` is reached before the event was triggered `trigger` times.
    /// * `PropagationError::TrajectoryEvent`: Returned if the interpolation of the event state fails.
    /// * `PropagationError::Analysis`: Returned if the event evaluation fails during the search.
    #[pyo3(signature = (spacecraft, event, max_duration, trigger=1, event_frame=None, trajectory=true))]
    fn until_event(
        &self,
        spacecraft: Spacecraft,
        event: &Event,
        max_duration: Duration,
        trigger: usize,
        event_frame: Option<Frame>,
        trajectory: bool,
    ) -> Result<PropagationResult, PropagationError> {
        let setup = self.build()?;

        let mut prop = setup.with(spacecraft, self.almanac.clone());

        let (state, traj) = prop.until_nth_event(max_duration, event, event_frame, trigger)?;

        Ok(PropagationResult {
            state,
            trajectory: if trajectory {
                Some(PyTrajectory { inner: traj })
            } else {
                None
            },
        })
    }

    /// Propagates the initialization state until the desired epoch, optionally not building the trajectory
    #[pyo3(signature = (spacecraft, epoch, trajectory=true))]
    fn many_until_epoch(
        &self,
        py: Python,
        spacecraft: Vec<Spacecraft>,
        epoch: Epoch,
        trajectory: bool,
    ) -> Result<Vec<PropagationResult>, PropagationError> {
        let setup = Arc::new(CorePropagator::new(
            self.dynamics
                .build(self.almanac.clone())
                .map_err(|msg| PropagationError::PropGenericError { msg })?,
            self.method,
            self.options,
        ));

        Ok(py.detach(|| {
            spacecraft
                .par_iter()
                .filter_map(|sc| {
                    let mut prop = setup.with(*sc, self.almanac.clone());
                    if trajectory {
                        match prop.until_epoch_with_traj(epoch) {
                            Ok((state, traj)) => Some(PropagationResult {
                                state,
                                trajectory: Some(PyTrajectory { inner: traj }),
                            }),
                            Err(e) => {
                                eprintln!("{e}");
                                None
                            }
                        }
                    } else {
                        match prop.until_epoch(epoch) {
                            Ok(state) => Some(PropagationResult {
                                state,
                                trajectory: None,
                            }),
                            Err(e) => {
                                eprintln!("{e}");
                                None
                            }
                        }
                    }
                })
                .collect::<Vec<PropagationResult>>()
        }))
    }

    /// Propagates the initialization state for the desired duration, optionally not building the trajectory
    #[pyo3(signature = (spacecraft, duration, trajectory=true))]
    fn many_for_duration(
        &self,
        py: Python,
        spacecraft: Vec<Spacecraft>,
        duration: Duration,
        trajectory: bool,
    ) -> Result<Vec<PropagationResult>, PropagationError> {
        let setup = Arc::new(CorePropagator::new(
            self.dynamics
                .build(self.almanac.clone())
                .map_err(|msg| PropagationError::PropGenericError { msg })?,
            self.method,
            self.options,
        ));

        Ok(py.detach(|| {
            spacecraft
                .par_iter()
                .filter_map(|sc| {
                    let mut prop = setup.with(*sc, self.almanac.clone());
                    if trajectory {
                        match prop.for_duration_with_traj(duration) {
                            Ok((state, traj)) => Some(PropagationResult {
                                state,
                                trajectory: Some(PyTrajectory { inner: traj }),
                            }),
                            Err(e) => {
                                eprintln!("{e}");
                                None
                            }
                        }
                    } else {
                        match prop.for_duration(duration) {
                            Ok(state) => Some(PropagationResult {
                                state,
                                trajectory: None,
                            }),
                            Err(e) => {
                                eprintln!("{e}");
                                None
                            }
                        }
                    }
                })
                .collect::<Vec<PropagationResult>>()
        }))
    }

    #[pyo3(signature = (spacecraft, event, max_duration, trigger=1, event_frame=None, trajectory=true))]
    fn many_until_event(
        &self,
        py: Python,
        spacecraft: Vec<Spacecraft>,
        event: &Event,
        max_duration: Duration,
        trigger: usize,
        event_frame: Option<Frame>,
        trajectory: bool,
    ) -> Result<Vec<PropagationResult>, PropagationError> {
        let setup = Arc::new(CorePropagator::new(
            self.dynamics
                .build(self.almanac.clone())
                .map_err(|msg| PropagationError::PropGenericError { msg })?,
            self.method,
            self.options,
        ));

        Ok(py.detach(|| {
            spacecraft
                .par_iter()
                .filter_map(|sc| {
                    let mut prop = setup.with(*sc, self.almanac.clone());
                    match prop.until_nth_event(max_duration, event, event_frame, trigger) {
                        Ok((state, traj)) => Some(PropagationResult {
                            state,
                            trajectory: if trajectory {
                                Some(PyTrajectory { inner: traj })
                            } else {
                                None
                            },
                        }),
                        Err(e) => {
                            eprintln!("{e}");
                            None
                        }
                    }
                })
                .collect::<Vec<PropagationResult>>()
        }))
    }
}

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

    /// Add another state to this trajectory.
    fn push(&mut self, spacecraft: Spacecraft) {
        self.inner.states.push(spacecraft);
        self.inner.finalize();
    }

    /// Append many spacecraft to this trajectory.
    fn append(&mut self, many_spacecraft: Vec<Spacecraft>) {
        for sc in many_spacecraft {
            self.inner.states.push(sc);
        }
        self.inner.finalize();
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

    fn __repr__(&self) -> String {
        format!("{} @ {self:p}", self.inner)
    }
}
