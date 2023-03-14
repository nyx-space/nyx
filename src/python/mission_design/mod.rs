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

use crate::io::trajectory_data::DynamicTrajectory;
use crate::io::ConfigError;
use crate::md::ui::{PropOpts, Propagator, SpacecraftDynamics};
use crate::md::{Event, StateParameter};

use crate::{NyxError, Spacecraft};
use hifitime::{Duration, Epoch, Unit};
use pyo3::{prelude::*, py_run};

use self::trajectory::Traj;

mod events;
pub mod spacecraft;
mod trajectory;

/// Mission design
pub(crate) fn register_md(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "nyx_space.mission_design")?;

    sm.add_class::<DynamicTrajectory>()?;
    sm.add_class::<SpacecraftDynamics>()?;
    sm.add_class::<StateParameter>()?;
    sm.add_class::<Event>()?;
    sm.add_function(wrap_pyfunction!(propagate, sm)?)?;

    py_run!(
        py,
        sm,
        "import sys; sys.modules['nyx_space.mission_design'] = sm"
    );
    parent_module.add_submodule(sm)?;
    Ok(())
}

/// Propagates the provided spacecraft with the provided dynamics until the provided stopping condition (duration, epoch, or event [and optionally the count]).
#[pyfunction]
#[pyo3(
    text_signature = "(spacecraft, dynamics, duration=None, epoch=None, event=None, event_count=None, min_step=None, max_step=None, fixed_step=None, tolerance=None)"
)]
fn propagate(
    spacecraft: Spacecraft,
    dynamics: SpacecraftDynamics,
    duration: Option<Duration>,
    epoch: Option<Epoch>,
    event: Option<Event>,
    event_count: Option<usize>,
    min_step: Option<Duration>,
    max_step: Option<Duration>,
    fixed_step: Option<Duration>,
    tolerance: Option<f64>,
) -> Result<(Spacecraft, Traj), NyxError> {
    let opts = match fixed_step {
        Some(step) => PropOpts::with_fixed_step(step),
        None => {
            let mut opts = PropOpts::default();
            if let Some(step) = min_step {
                opts.set_min_step(step);
            }
            if let Some(step) = max_step {
                opts.set_max_step(step);
            }
            if let Some(tol) = tolerance {
                opts.tolerance = tol;
            }
            opts
        }
    };
    info!("Propagator options: {opts}");

    if let Some(event) = event {
        let max_duration = match duration {
            Some(duration) => duration,
            None => {
                warn!("No maximum duration provided to search, setting to 30 days");
                30 * Unit::Day
            }
        };

        let (sc, traj) = match event_count {
            Some(count) => Propagator::rk89(dynamics, opts)
                .with(spacecraft)
                .until_nth_event(max_duration, &event, count)?,
            None => Propagator::rk89(dynamics, opts)
                .with(spacecraft)
                .until_event(max_duration, &event)?,
        };

        Ok((sc, Traj { inner: traj }))
    } else if let Some(duration) = duration {
        let (sc, traj) = Propagator::rk89(dynamics, opts)
            .with(spacecraft)
            .for_duration_with_traj(duration)?;

        Ok((sc, Traj { inner: traj }))
    } else if let Some(epoch) = epoch {
        let (sc, traj) = Propagator::rk89(dynamics, opts)
            .with(spacecraft)
            .until_epoch_with_traj(epoch)?;

        Ok((sc, Traj { inner: traj }))
    } else {
        Err(NyxError::ConfigError(ConfigError::InvalidConfig(
            "Either duration or epoch must be provided for a propagation to happen".to_string(),
        )))
    }
}
