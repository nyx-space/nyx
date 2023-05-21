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
use crate::io::{ConfigError, ExportCfg};
use crate::md::prelude::{PropOpts, Propagator, SpacecraftDynamics};
use crate::md::{Event, StateParameter};

use crate::propagators::{
    CashKarp45, Dormand45, Dormand78, Fehlberg45, RK2Fixed, RK4Fixed, Verner56,
};
use crate::{NyxError, Spacecraft};
use hifitime::{Duration, Epoch, Unit};
use pyo3::{prelude::*, py_run};

use self::trajectory::Traj;

mod events;
pub mod spacecraft;
mod trajectory;

/// Mission design
pub(crate) fn register_md(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "_nyx_space.mission_design")?;

    sm.add_class::<DynamicTrajectory>()?;
    sm.add_class::<SpacecraftDynamics>()?;
    sm.add_class::<StateParameter>()?;
    sm.add_class::<Event>()?;
    sm.add_class::<ExportCfg>()?;
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
///
/// Available methods: rk89, dormand78, dormand45, rk45 (or fehlberg45), cashkarp45, verner56, rk4, rk2
#[pyfunction]
#[pyo3(
    text_signature = "(spacecraft, dynamics, duration=None, epoch=None, event=None, event_count=None, min_step=None, max_step=None, fixed_step=None, tolerance=None, method='rk89')"
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
    method: Option<String>,
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

    let prop_setup = match method {
        Some(value) => match value.to_lowercase().as_str() {
            "rk89" => Propagator::rk89(dynamics, opts),
            "dormand78" => Propagator::new::<Dormand78>(dynamics, opts),
            "dormand45" => Propagator::new::<Dormand45>(dynamics, opts),
            "rk45" | "fehlberg45" => Propagator::new::<Fehlberg45>(dynamics, opts),
            "cashkarp45" => Propagator::new::<CashKarp45>(dynamics, opts),
            "verner56" => Propagator::new::<Verner56>(dynamics, opts),
            "rk4" => Propagator::new::<RK4Fixed>(dynamics, opts),
            "rk2" => Propagator::new::<RK2Fixed>(dynamics, opts),
            _ => {
                return Err(NyxError::ConfigError(ConfigError::InvalidConfig(format!(
                    "Unknown propagation method: {}",
                    value
                ))))
            }
        },
        None => Propagator::rk89(dynamics, opts),
    };

    if let Some(event) = event {
        let max_duration = match duration {
            Some(duration) => duration,
            None => {
                warn!("No maximum duration provided to search, setting to 30 days");
                30 * Unit::Day
            }
        };

        let (sc, traj) = match event_count {
            Some(count) => {
                prop_setup
                    .with(spacecraft)
                    .until_nth_event(max_duration, &event, count)?
            }
            None => prop_setup
                .with(spacecraft)
                .until_event(max_duration, &event)?,
        };

        Ok((sc, Traj { inner: traj }))
    } else if let Some(duration) = duration {
        let (sc, traj) = prop_setup
            .with(spacecraft)
            .for_duration_with_traj(duration)?;

        Ok((sc, Traj { inner: traj }))
    } else if let Some(epoch) = epoch {
        let (sc, traj) = prop_setup.with(spacecraft).until_epoch_with_traj(epoch)?;

        Ok((sc, Traj { inner: traj }))
    } else {
        Err(NyxError::ConfigError(ConfigError::InvalidConfig(
            "Either duration or epoch must be provided for a propagation to happen".to_string(),
        )))
    }
}
