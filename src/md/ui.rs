/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

extern crate csv;
extern crate rayon;

use self::rayon::prelude::*;
use super::MdHdlr;
pub use crate::celestia::*;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, U6};
pub use crate::dynamics::orbital::OrbitalDynamics;
use crate::dynamics::solarpressure::SolarPressure;
use crate::dynamics::spacecraft::Spacecraft;
use crate::dynamics::sph_harmonics::Harmonics;
pub use crate::dynamics::{Dynamics, NyxError};
use crate::io::formatter::*;
use crate::io::quantity::ParsingError;
use crate::io::scenario::ConditionSerde;
use crate::io::scenario::ScenarioSerde;
use crate::propagators::Propagator;
use crate::time::{Duration, Epoch, SECONDS_PER_DAY};
use crate::SpacecraftState;
use std::convert::TryFrom;
use std::str::FromStr;
use std::sync::{mpsc::channel, Arc};
use std::time::Instant;

/// An MDProcess allows the creation and propagation of a spacecraft subjected to some dynamics
pub struct MDProcess<'a>
where
    DefaultAllocator: Allocator<f64, U6>,
{
    pub sc_dyn: Arc<Spacecraft<'a>>,
    pub init_state: SpacecraftState,
    pub formatter: Option<StateFormatter>,
    pub prop_time: Option<Duration>,
    pub prop_event: Option<ConditionSerde>,
    pub prop_tol: f64,
    pub name: String,
}

impl<'a> MDProcess<'a>
where
    DefaultAllocator: Allocator<f64, U6>,
{
    pub fn try_from_scenario(
        scen: &ScenarioSerde,
        prop_name: String,
        stm_flag: bool,
        cosm: Arc<Cosm>,
    ) -> Result<(Self, Option<StateFormatter>), ParsingError> {
        match scen.propagator.get(&prop_name.to_lowercase()) {
            None => Err(ParsingError::PropagatorNotFound(prop_name)),
            Some(prop) => {
                #[allow(unused_assignments)]
                let mut sc_dyn: Spacecraft;
                let init_sc;

                // Validate the output
                let formatter = if let Some(output) = &prop.output {
                    match scen.output.get(&output.to_lowercase()) {
                        None => {
                            return Err(ParsingError::MD(format!(
                                "propagator `{}` refers to undefined output `{}`",
                                prop_name, output
                            )))
                        }
                        Some(out) => Some(out.to_state_formatter(cosm.clone())),
                    }
                } else {
                    None
                };
                let spacecraft = match scen.spacecraft.get(&prop.dynamics.to_lowercase()) {
                    None => {
                        return Err(ParsingError::MD(format!(
                            "propagator `{}` refers to undefined spacecraft `{}`",
                            prop_name, prop.dynamics
                        )))
                    }
                    Some(spacecraft) => spacecraft,
                };

                // Get and validate the orbital dynamics
                let dynamics = match scen
                    .orbital_dynamics
                    .get(&spacecraft.orbital_dynamics.to_lowercase())
                {
                    None => {
                        return Err(ParsingError::MD(format!(
                            "spacecraft `{}` refers to undefined dynamics `{}`",
                            &prop.dynamics, spacecraft.orbital_dynamics
                        )))
                    }
                    Some(dynamics) => dynamics,
                };

                let state_name = &dynamics.initial_state.to_lowercase();
                let init_state = match scen.state.get(state_name) {
                    None => {
                        // Let's see if there is a delta state
                        if scen.delta_state.is_some() {
                            match scen
                                .delta_state
                                .as_ref()
                                .unwrap()
                                .get(&dynamics.initial_state.to_lowercase())
                            {
                                None => {
                                    return Err(ParsingError::MD(format!(
                                        "dynamics `{}` refers to unknown state `{}`",
                                        prop.dynamics, dynamics.initial_state
                                    )))
                                }
                                Some(delta_state) => {
                                    // Rebuilt the state
                                    let inherited = match scen.state.get(&delta_state.inherit) {
                                        None => {
                                            return Err(ParsingError::MD(format!(
                                                "delta state `{}` refers to unknown state `{}`",
                                                state_name, delta_state.inherit
                                            )))
                                        }
                                        Some(base) => {
                                            let state_frame = &cosm.frame(base.frame.as_str());
                                            base.as_state(*state_frame)?
                                        }
                                    };
                                    delta_state.as_state(inherited)?
                                }
                            }
                        } else {
                            return Err(ParsingError::MD(format!(
                                "dynamics `{}` refers to unknown state `{}`",
                                prop.dynamics, dynamics.initial_state
                            )));
                        }
                    }
                    Some(init_state_sd) => {
                        let state_frame = &cosm.frame(init_state_sd.frame.as_str());
                        let mut init = init_state_sd.as_state(*state_frame)?;
                        if stm_flag {
                            // Specify that we want to compute the STM.
                            init.enable_stm();
                        }
                        init
                    }
                };

                // Create the dynamics
                if let Some(pts_masses) = &dynamics.point_masses {
                    // Get the object IDs from name
                    let mut bodies = Vec::with_capacity(10);
                    for obj in pts_masses {
                        match Bodies::try_from(obj.to_string()) {
                            Ok(b) => bodies.push(b),
                            Err(e) => {
                                return Err(ParsingError::LoadingError(format!("{:?}", e)));
                            }
                        }
                    }
                    // Remove bodies which are part of the state
                    if let Some(pos) = bodies
                        .iter()
                        .position(|x| x.exb_id() == init_state.frame.exb_id())
                    {
                        bodies.remove(pos);
                    }

                    sc_dyn = Spacecraft::new_raw(OrbitalDynamics::point_masses(
                        init_state.frame,
                        &bodies,
                        cosm.clone(),
                    ));

                    init_sc = SpacecraftState::new(
                        init_state,
                        spacecraft.dry_mass,
                        spacecraft.fuel_mass.unwrap_or(0.0),
                    );
                } else {
                    sc_dyn = Spacecraft::new_raw(OrbitalDynamics::two_body());

                    init_sc = SpacecraftState::new(
                        init_state,
                        spacecraft.dry_mass,
                        spacecraft.fuel_mass.unwrap_or(0.0),
                    );

                    // Add the force models
                    if let Some(force_models) = &spacecraft.force_models {
                        if scen.force_models.as_ref().is_none() {
                            return Err(ParsingError::MD(format!(
                                "spacecraft `{}` refers to force models but none are defined",
                                prop.dynamics
                            )));
                        }
                        for mdl in force_models {
                            match scen.force_models.as_ref().unwrap().get(&mdl.to_lowercase()) {
                                None => {
                                    return Err(ParsingError::MD(format!(
                                        "spacecraft `{}` refers to unknown force model `{}`",
                                        prop.dynamics, mdl
                                    )))
                                }
                                Some(amdl) => {
                                    let eme2k = &cosm.frame("EME2000");
                                    let luna = &cosm.frame("Luna");
                                    for smdl in amdl.srp.values() {
                                        // Note that an Arc is immutable, but we want to specify everything
                                        // so we create the SRP without the wrapper
                                        let mut srp = SolarPressure::default_raw(
                                            smdl.sc_area,
                                            vec![*eme2k, *luna],
                                            cosm.clone(),
                                        );
                                        srp.phi = smdl.phi;
                                        srp.cr = smdl.cr;
                                        sc_dyn.add_model(Arc::new(srp));
                                    }
                                }
                            }
                        }
                    }
                }

                // Add the acceleration models if applicable
                if let Some(accel_models) = &dynamics.accel_models {
                    // In this case, we'll need to recreate the orbital dynamics because it's behind an immutable Arc.
                    let mut orbital_dyn = OrbitalDynamics::new_raw(vec![]);
                    for mdl in accel_models {
                        match scen.accel_models.as_ref().unwrap().get(&mdl.to_lowercase()) {
                            None => {
                                return Err(ParsingError::MD(format!(
                                    "dynamics `{}` refers to unknown state `{}`",
                                    prop.dynamics, dynamics.initial_state
                                )))
                            }
                            Some(amdl) => {
                                for hmdl in amdl.harmonics.values() {
                                    let in_mem = hmdl.load();
                                    let compute_frame = &cosm.frame(hmdl.frame.as_str());

                                    let hh =
                                        Harmonics::from_stor(*compute_frame, in_mem, cosm.clone());
                                    orbital_dyn.add_model(hh);
                                }
                            }
                        }
                    }
                    // And set these into the spacecraft dynamics
                    sc_dyn.orbital_dyn = Arc::new(orbital_dyn);
                }

                // Validate the stopping condition
                // Check if it's a stopping condition
                let prop_event = if let Some(conditions) = &scen.conditions {
                    match conditions.get(&prop.stop_cond) {
                        Some(c) => Some(c.clone()),
                        None => None,
                    }
                } else {
                    None
                };
                // Let's see if it's a relative time
                let prop_time = if prop_event.is_none() {
                    match Duration::from_str(prop.stop_cond.as_str()) {
                        Ok(d) => Some(d),
                        Err(_) => match Epoch::from_str(prop.stop_cond.as_str()) {
                            Err(e) => return Err(ParsingError::IllDefined(format!("{}", e))),
                            Ok(epoch) => Some(epoch - init_state.dt),
                        },
                    }
                } else {
                    None
                };

                // Let's see if it's a relative time
                let prop_tol = prop.tolerance.unwrap_or(1e-12);

                Ok((
                    Self {
                        sc_dyn: Arc::new(sc_dyn),
                        init_state: init_sc,
                        formatter: None,
                        prop_time,
                        prop_event,
                        prop_tol,
                        name: prop_name.clone(),
                    },
                    formatter,
                ))
            }
        }
    }

    // #[inline]
    // pub fn propagator(&'a self) -> PropInstance<'a, Spacecraft<'a>, RSSStepPV> {
    //     let mut p = Propagator::default(&self.sc_dyn);
    //     p.set_tolerance(self.prop_tol);
    //     p.with(self.init_state)
    // }

    pub fn execute(mut self) -> Result<(), NyxError> {
        self.execute_with(vec![])
    }

    /// Execute the MD with the provided handlers. Note that you must initialize your own CSV output if that's desired.
    pub fn execute_with(
        &mut self,
        mut hdlrs: Vec<Box<dyn MdHdlr<SpacecraftState>>>,
    ) -> Result<(), NyxError> {
        // Get the prop time before the mutable ref of the propagator
        let maybe_prop_time = self.prop_time;
        let maybe_prop_event = self.prop_event.clone();

        // Build the propagator
        let mut prop_setup = Propagator::default(self.sc_dyn.clone());
        prop_setup.set_tolerance(self.prop_tol);
        let mut prop = prop_setup.with(self.init_state);
        // let mut prop = self.propagator();
        // Set up the channels
        let (tx, rx) = channel();
        prop.tx_chan = Some(tx);

        let mut initial_state = Some(prop.state);

        info!("Initial state: {}", prop.state);
        let start = Instant::now();

        // Run
        match maybe_prop_event {
            Some(prop_event) => {
                let stop_cond = prop_event.to_condition(initial_state.unwrap().orbit.dt);
                info!("Propagating until event {:?}", stop_cond.event);
                let rslt = prop.until_event(stop_cond);
                if rslt.is_err() {
                    panic!("{:?}", rslt.err());
                }
            }
            None => {
                // Elapsed seconds propagation
                let prop_time = maybe_prop_time.unwrap();
                info!(
                    "Propagating for {} seconds (~ {:.3} days)",
                    prop_time,
                    prop_time / SECONDS_PER_DAY
                );
                prop.for_duration(prop_time)?;
            }
        }

        info!(
            "Final state:   {} (computed in {:.3} seconds)",
            prop.state,
            (Instant::now() - start).as_secs_f64()
        );

        while let Ok(prop_state) = rx.try_recv() {
            // Provide to the handler
            hdlrs.par_iter_mut().for_each(|hdlr| {
                // let mut hdlr = hdlr_am.lock().unwrap();
                if let Some(first_state) = initial_state {
                    hdlr.handle(&first_state);
                }
                hdlr.handle(&prop_state);
            });

            // We've used the initial state (if it was there)
            if initial_state.is_some() {
                initial_state = None;
            }
        }

        Ok(())
    }
}
