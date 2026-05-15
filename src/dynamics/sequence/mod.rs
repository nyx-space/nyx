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
use hifitime::{Epoch, Unit};
use indexmap::IndexMap;
use log::{debug, info};
use serde::{Deserialize, Serialize};
#[cfg(feature = "python")]
use pyo3::prelude::*;
use serde_dhall::{SimpleType, StaticType};
use snafu::ResultExt;
use std::collections::HashMap;

use crate::dynamics::guidance::{Kluever, Ruggiero};
use crate::dynamics::{guidance::Thruster, SpacecraftDynamics};
use crate::dynamics::{GravityField, OrbitalDynamics};
use crate::errors::{FromAlmanacSnafu, FromPropSnafu};
use crate::io::gravity::GravityFieldData;
use crate::md::Trajectory;
use crate::propagators::Propagator;
use crate::{NyxError, Spacecraft, State};

mod config;
mod discrete_event;

pub use config::*;
pub use discrete_event::*;

#[derive(Clone, Debug, Default, Serialize, Deserialize)]
#[cfg_attr(feature = "python", pyclass)]
pub struct SpacecraftSequence {
    #[serde(serialize_with = "map_as_pairs", deserialize_with = "pairs_as_map")]
    pub seq: BTreeMap<Epoch, Phase>,
    #[serde(with = "indexmap::map::serde_seq")]
    pub thruster_sets: IndexMap<String, Thruster>,
    #[serde(with = "indexmap::map::serde_seq")]
    pub propagators: IndexMap<String, PropagatorConfig>,
    #[serde(skip)]
    prop_setups: IndexMap<String, Propagator<SpacecraftDynamics>>,
}

impl SpacecraftSequence {
    pub fn validate(&self) -> Result<(), String> {
        // Check that the last statement is a terminate
        if let Some((_, Phase::Activity { .. })) = self.seq.iter().last() {
            return Err("final phase must be a Terminate".into());
        }

        // Check that all of the thruster set indexes reference an available thruster
        for (epoch, phase) in &self.seq {
            if let Phase::Activity {
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
                    let thruster = &guidance.thruster_model;
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

        for phase in self.seq.values() {
            if let Phase::Activity {
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
                        orbital_dyn.accel_models.push(point_masses.clone());
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
                        sc_dyn.force_models.push(srp.clone());
                    }

                    if let Some(drag) = &cfg.force_models.drag {
                        sc_dyn.force_models.push(drag.clone());
                    }

                    // And set it all up!
                    let setup = Propagator::new(sc_dyn, cfg.method, cfg.options);

                    self.prop_setups.insert(propagator.clone(), setup);
                    debug!("built `{propagator}`");
                }
            }
        }

        Ok(())
    }

    /// Propagate this plan starting the relevant phase for the probided state, and propagating until the end of the plan
    /// or until the provided phase name. Returns the trajectory for each phase, allowing for each phase to be in its own central body.
    pub fn propagate(
        &self,
        mut state: Spacecraft,
        until_phase: Option<String>,
        almanac: Arc<Almanac>,
    ) -> Result<Vec<Trajectory>, NyxError> {
        let tick = Epoch::now().unwrap();
        let mut phase_iterator = self.seq.range(state.epoch()..).peekable();

        let mut trajs = Vec::with_capacity(self.seq.len());

        while let Some((epoch, phase)) = phase_iterator.next() {
            match phase {
                Phase::Terminate => {
                    let tock = (Epoch::now().unwrap() - tick).round(Unit::Millisecond * 1);
                    info!("[{epoch}] plan completed in {tock}");
                    return Ok(trajs);
                }
                Phase::Activity {
                    name,
                    propagator,
                    guidance,
                    on_entry,
                    disabled,
                } => {
                    // Check stop condition
                    if let Some(ref target) = until_phase {
                        if target == name {
                            return Ok(trajs);
                        }
                    }

                    if *disabled {
                        info!("[{epoch}] skipping disabled {name}");
                    } else {
                        info!("[{epoch}] executing {name}");
                        if let Some(discrete_event) = on_entry {
                            match &**discrete_event {
                                DiscreteEvent::FrameSwap { new_frame } => {
                                    if !new_frame.orient_origin_match(state.orbit.frame)
                                        || !new_frame.ephem_origin_match(state.orbit.frame)
                                    {
                                        state = state.with_orbit(
                                            almanac
                                                .translate_to(state.orbit, *new_frame, None)
                                                .map_err(|source| {
                                                    anise::errors::AlmanacError::Ephemeris {
                                                        action: "central body swap",
                                                        source: Box::new(source),
                                                    }
                                                })
                                                .context(FromAlmanacSnafu {
                                                    action: "central body swap",
                                                })?,
                                        );
                                        info!("[{epoch}] central body swapped to {new_frame}");
                                    }
                                }
                                DiscreteEvent::Staging {
                                    impulsive_maneuver,
                                    decrement_properties,
                                } => {
                                    if let Some(mnvr) = impulsive_maneuver {
                                        info!("[{epoch}] staging, with maneuver {mnvr}");
                                        state = state
                                            .with_orbit(state.orbit.with_dv_km_s(mnvr.dv_km_s));
                                    }
                                    if let Some(decr) = decrement_properties {
                                        if let Some(mass) = decr.mass {
                                            state.mass.dry_mass_kg -= mass.dry_mass_kg;
                                            state.mass.prop_mass_kg -= mass.prop_mass_kg;
                                            state.mass.extra_mass_kg -= mass.extra_mass_kg;
                                        }
                                        if let Some(srp) = decr.srp {
                                            state.srp.area_m2 -= srp.area_m2;
                                            state.srp.coeff_reflectivity -= srp.coeff_reflectivity;
                                        }
                                        if let Some(drag) = decr.drag {
                                            state.drag.area_m2 -= drag.area_m2;
                                            state.drag.coeff_drag -= drag.coeff_drag;
                                        }
                                    }
                                }
                                DiscreteEvent::Docking {
                                    impulsive_maneuver,
                                    increment_properties,
                                } => {
                                    if let Some(mnvr) = impulsive_maneuver {
                                        info!("[{epoch}] docking, with maneuver {mnvr}");
                                        state = state
                                            .with_orbit(state.orbit.with_dv_km_s(mnvr.dv_km_s));
                                    }
                                    if let Some(incr) = increment_properties {
                                        if let Some(mass) = incr.mass {
                                            state.mass.dry_mass_kg += mass.dry_mass_kg;
                                            state.mass.prop_mass_kg += mass.prop_mass_kg;
                                            state.mass.extra_mass_kg += mass.extra_mass_kg;
                                        }
                                        if let Some(srp) = incr.srp {
                                            state.srp.area_m2 += srp.area_m2;
                                            state.srp.coeff_reflectivity += srp.coeff_reflectivity;
                                        }
                                        if let Some(drag) = incr.drag {
                                            state.drag.area_m2 += drag.area_m2;
                                            state.drag.coeff_drag += drag.coeff_drag;
                                        }
                                    }
                                }
                            }
                        }

                        let end_time = phase_iterator
                            .peek()
                            .expect("validate did not catch missing terminate")
                            .0;

                        // Include the guidance if any is available
                        let (next_state, mut phase_traj) = if let Some(guid_cfg) = guidance {
                            // Clone the propagator to add the dynamics
                            let mut setup = self.prop_setups[propagator].clone();
                            setup.dynamics.decrement_mass = !guid_cfg.disable_prop_mass;
                            state.thruster = Some(self.thruster_sets[&guid_cfg.thruster_model]);
                            match &guid_cfg.law {
                                SteeringLaw::FiniteBurn(maneuver) => {
                                    setup.dynamics.guid_law = Some(Arc::new(*maneuver));
                                }
                                SteeringLaw::Ruggiero {
                                    objectives,
                                    max_eclipse_prct,
                                } => {
                                    let guid = Ruggiero {
                                        objectives: objectives.clone(),
                                        max_eclipse_prct: *max_eclipse_prct,
                                        init_state: state,
                                    };
                                    setup.dynamics.guid_law = Some(Arc::new(guid));
                                }
                                SteeringLaw::Kluever {
                                    objectives,
                                    max_eclipse_prct,
                                } => {
                                    let guid = Kluever {
                                        objectives: objectives.clone(),
                                        max_eclipse_prct: *max_eclipse_prct,
                                    };
                                    setup.dynamics.guid_law = Some(Arc::new(guid));
                                }
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
                        info!("[{epoch}] {name} completed: {next_state:x}");
                        state = next_state;
                        phase_traj.name = Some(name.clone());
                        trajs.push(phase_traj);
                    }
                }
            }
        }
        unreachable!("spacecraft plan never finished?!")
    }
}

impl StaticType for SpacecraftSequence {
    fn static_type() -> serde_dhall::SimpleType {
        let mut repr = HashMap::new();

        // seq maps to "sequence" in Dhall
        // Serialized as List { _1: Text, _2: Phase }
        let mut seq_entry = HashMap::new();
        seq_entry.insert("_1".to_string(), SimpleType::Text); // Epoch serializes to Text
        seq_entry.insert("_2".to_string(), Phase::static_type());

        repr.insert(
            "seq".to_string(),
            SimpleType::List(Box::new(SimpleType::Record(seq_entry))),
        );

        // thruster_sets maps to "thruster_set" in Dhall (matches your serialization name)
        // Serialized as List { _1: Text, _2: Thruster }
        let mut thruster_sets = HashMap::new();
        thruster_sets.insert("_1".to_string(), SimpleType::Text);
        thruster_sets.insert("_2".to_string(), Thruster::static_type());

        repr.insert(
            "thruster_sets".to_string(), // Keep as "thruster_set" if that's your Dhall preference
            SimpleType::List(Box::new(SimpleType::Record(thruster_sets))),
        );

        let mut propagators = HashMap::new();
        propagators.insert("_1".to_string(), SimpleType::Text);
        propagators.insert("_2".to_string(), PropagatorConfig::static_type());

        repr.insert(
            "propagators".to_string(), // Keep as "thruster_set" if that's your Dhall preference
            SimpleType::List(Box::new(SimpleType::Record(propagators))),
        );

        SimpleType::Record(repr)
    }
}

/* serialization helper functions */

fn map_as_pairs<S, K, V>(map: &BTreeMap<K, V>, serializer: S) -> Result<S::Ok, S::Error>
where
    S: serde::Serializer,
    K: Serialize + Clone,
    V: Serialize + Clone,
{
    // This turns the map into a sequence of (K, V) which Serde-Dhall sees as {_1, _2}
    serializer.collect_seq(map.iter())
}

// You'll need a symmetric deserializer if you plan to read these back
fn pairs_as_map<'de, D, K, V>(deserializer: D) -> Result<BTreeMap<K, V>, D::Error>
where
    D: serde::Deserializer<'de>,
    K: Deserialize<'de> + Ord,
    V: Deserialize<'de>,
{
    let pairs: Vec<(K, V)> = Vec::deserialize(deserializer)?;
    Ok(pairs.into_iter().collect())
}
pub mod python;
