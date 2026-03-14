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

use std::collections::BTreeMap;
use std::sync::Arc;

use anise::prelude::Almanac;
use hifitime::Epoch;
use indexmap::IndexMap;
use log::info;
use snafu::ResultExt;

use crate::dynamics::{guidance::Thruster, SpacecraftDynamics};
use crate::dynamics::{GravityField, OrbitalDynamics};
use crate::errors::FromPropSnafu;
use crate::io::gravity::GravityFieldData;
use crate::md::Trajectory;
use crate::propagators::Propagator;
use crate::{NyxError, Spacecraft, State};

pub mod config;
pub mod discrete_event;

use config::*;
use discrete_event::*;

#[derive(Clone, Debug)]
pub struct SpacecraftPlan {
    pub timeline: BTreeMap<Epoch, TimelinePhase>,
    pub thruster_sets: IndexMap<String, Thruster>,
    pub propagators: IndexMap<String, PropagatorConfig>,
    prop_setups: IndexMap<String, Propagator<SpacecraftDynamics>>,
}

impl SpacecraftPlan {
    pub fn validate(&self) -> Result<(), String> {
        // Check that the last statement is a terminate
        if let Some((_, TimelinePhase::Phase { .. })) = self.timeline.iter().last() {
            return Err("final phase must be a Terminate".into());
        }

        // Check that all of the thruster set indexes reference an available thruster
        for (epoch, phase) in &self.timeline {
            if let TimelinePhase::Phase {
                name: _,
                propagator,
                guidance,
                on_entry: _,
                disabled: _,
            } = phase
            {
                // Check that the propagator exists
                if self.propagators.get(propagator).is_none() {
                    return Err(format!("{epoch}: no propagator named `{propagator}`"));
                }
                if let Some(guidance) = guidance {
                    let thruster = guidance.thruster_model();
                    if self.thruster_sets.get(thruster).is_none() {
                        return Err(format!("{epoch}: no thruster set named {thruster}"));
                    }
                }
            }
        }

        Ok(())
    }

    /// Set up the propagators that are used in the timeline
    pub fn setup(&mut self, almanac: Arc<Almanac>) -> Result<(), String> {
        // Don't set up anything if this is not a valid timeline
        self.validate()?;

        for phase in self.timeline.values() {
            if let TimelinePhase::Phase {
                name: _,
                propagator,
                guidance: _,
                on_entry: _,
                disabled,
            } = phase
            {
                if !disabled && self.prop_setups.get(propagator).is_none() {
                    // Set up the propagator -- fetch the config first
                    // We know the config exists because validate would catch missing names.
                    let cfg = &self.propagators[propagator];
                    // Build the orbital dynamics
                    let mut orbital_dyn = OrbitalDynamics::two_body();
                    if let Some(point_masses) = &cfg.accel_models.point_masses {
                        orbital_dyn
                            .accel_models
                            .push(Arc::new(point_masses.clone()));
                    }
                    if let Some((gravity_cfg, frame_uid)) = &cfg.accel_models.gravity_field {
                        let grav_data = GravityFieldData::from_config(gravity_cfg.clone())
                            .map_err(|e| e.to_string())?;
                        let compute_frame =
                            almanac.frame_info(*frame_uid).map_err(|e| e.to_string())?;
                        let gravity_field = GravityField::from_stor(compute_frame, grav_data);
                        orbital_dyn.accel_models.push(gravity_field);
                    }
                    // Build the spacecraft dynamics
                    let mut sc_dyn = SpacecraftDynamics::new(orbital_dyn);

                    if let Some(srp) = &cfg.force_models.solar_pressure {
                        sc_dyn.force_models.push(Arc::new(srp.clone()));
                    }

                    if let Some(drag) = &cfg.force_models.drag {
                        sc_dyn.force_models.push(Arc::new(drag.clone()));
                    }

                    // And set it all up!
                    let setup = Propagator::new(sc_dyn, cfg.method, cfg.options);

                    self.prop_setups.insert(propagator.clone(), setup);
                }
            }
        }

        Ok(())
    }

    /// Propagate this plan starting the relevant phase for the probided state, and propagating until the end of the plan
    /// or until the provided phase name.
    pub fn propagate(
        &self,
        mut state: Spacecraft,
        until_phase: Option<String>,
        almanac: Arc<Almanac>,
    ) -> Result<Trajectory, NyxError> {
        let mut phase_iterator = self.timeline.range(state.epoch()..).peekable();

        let mut traj = Trajectory::new();

        while let Some((epoch, phase)) = phase_iterator.next() {
            match phase {
                TimelinePhase::Terminate => {
                    info!("[{epoch}] plan completed");
                    return Ok(traj);
                }
                TimelinePhase::Phase {
                    name,
                    propagator,
                    guidance,
                    on_entry,
                    disabled,
                } => {
                    // Check stop condition
                    if let Some(ref target) = until_phase {
                        if target == name {
                            return Ok(traj);
                        }
                    }

                    if !disabled {
                        info!("[{epoch}] executing {name}");
                        // TODO: Apply the on_entry condition

                        let end_time = phase_iterator
                            .peek()
                            .expect("validate did not catch missing terminate")
                            .0;

                        // Include the guidance if any is available
                        let (next_state, phase_traj) = if let Some(guid_cfg) = guidance {
                            // Clone the propagator to add the dynamics
                            let mut setup = self.prop_setups[propagator].clone();
                            match &**guid_cfg {
                                GuidanceConfig::FiniteBurn {
                                    maneuver,
                                    thruster_model: _,
                                } => {
                                    setup.dynamics.guid_law = Some(Arc::new(*maneuver));
                                }
                                _ => unimplemented!(),
                            }
                            setup
                                .with(state, almanac.clone())
                                .until_epoch_with_traj(*end_time)
                                .context(FromPropSnafu)?
                        } else {
                            self.prop_setups[propagator]
                                .with(state, almanac.clone())
                                .until_epoch_with_traj(*end_time)
                                .context(FromPropSnafu)?
                        };
                        state = next_state;
                        // Append all the states
                        phase_traj
                            .states
                            .iter()
                            .copied()
                            .for_each(|sc_state| traj.states.push(sc_state));
                        traj.finalize();
                    }
                }
            }
        }
        unreachable!("spacecraft plan never finished?!")
    }
}
