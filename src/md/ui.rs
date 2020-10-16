extern crate csv;
extern crate rayon;

use self::rayon::prelude::*;
use super::MdHdlr;
pub use crate::celestia::*;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, U6};
pub use crate::dynamics::orbital::{OrbitalDynamics, OrbitalDynamicsStm, OrbitalDynamicsT};
use crate::dynamics::solarpressure::SolarPressure;
use crate::dynamics::spacecraft::{Spacecraft, SpacecraftState};
use crate::dynamics::sph_harmonics::{Harmonics, HarmonicsDiff};
pub use crate::dynamics::{Dynamics, NyxError};
use crate::io::formatter::*;
use crate::io::quantity::{parse_duration, ParsingError};
use crate::io::scenario::ConditionSerde;
use crate::io::scenario::ScenarioSerde;
use crate::propagators::error_ctrl::RSSStepPV;
use crate::propagators::{PropOpts, Propagator};
use crate::time::{Epoch, SECONDS_PER_DAY};
use std::str::FromStr;
use std::sync::mpsc::channel;
use std::time::Instant;

pub enum StmState<N, M> {
    Without(N),
    With(M),
    Unknown,
}

impl<N, M> StmState<N, M> {
    pub fn with(&self) -> bool {
        matches!(self, StmState::With(_))
    }
}

pub type StmStateFlag = StmState<(), ()>;

/// An MDProcess allows the creation and propagation of a spacecraft subjected to some dynamics
pub struct MDProcess<'a>
where
    DefaultAllocator: Allocator<f64, U6>,
{
    sc_dyn: StmState<Spacecraft<'a, OrbitalDynamics<'a>>, Spacecraft<'a, OrbitalDynamicsStm<'a>>>,
    pub formatter: Option<StateFormatter<'a>>,
    pub prop_time_s: Option<f64>,
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
        stm_flag: StmStateFlag,
        cosm: &'a Cosm,
    ) -> Result<(Self, Option<StateFormatter<'a>>), ParsingError> {
        match scen.propagator.get(&prop_name.to_lowercase()) {
            None => Err(ParsingError::PropagatorNotFound(prop_name)),
            Some(prop) => {
                #[allow(unused_assignments)]
                let mut sc_dyn_flagged = StmState::Unknown;

                // Validate the output
                let formatter = if let Some(output) = &prop.output {
                    match scen.output.get(&output.to_lowercase()) {
                        None => {
                            return Err(ParsingError::MD(format!(
                                "propagator `{}` refers to undefined output `{}`",
                                prop_name, output
                            )))
                        }
                        Some(out) => Some(out.to_state_formatter(cosm)),
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
                                            let state_frame = cosm.frame(base.frame.as_str());
                                            base.as_state(state_frame)?
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
                        let state_frame = cosm.frame(init_state_sd.frame.as_str());
                        init_state_sd.as_state(state_frame)?
                    }
                };

                // Create the dynamics
                if let Some(pts_masses) = &dynamics.point_masses {
                    // Get the object IDs from name
                    let mut bodies = Vec::with_capacity(10);
                    for obj in pts_masses {
                        match cosm.try_frame(obj) {
                            Ok(frame) => bodies.push(frame.exb_id()),
                            Err(_) => {
                                // Let's try with "j2000" appended
                                match cosm.try_frame(format!("{} j2000", obj).as_str()) {
                                    Ok(frame) => bodies.push(frame.exb_id()),
                                    Err(_) => {
                                        bodies.push(
                                            cosm.frame(
                                                format!("{} barycenter j2000", obj).as_str(),
                                            )
                                            .exb_id(),
                                        );
                                    }
                                }
                            }
                        }
                    }
                    // Remove bodies which are part of the state
                    if let Some(pos) = bodies.iter().position(|x| *x == init_state.frame.exb_id()) {
                        bodies.remove(pos);
                    }
                    if stm_flag.with() {
                        let mut sc_dyn_stm = Spacecraft::with_stm(
                            OrbitalDynamicsStm::point_masses(init_state, bodies, cosm),
                            spacecraft.dry_mass,
                        );
                        if let Some(fuel_mass) = spacecraft.fuel_mass {
                            sc_dyn_stm.fuel_mass = fuel_mass;
                        }

                        sc_dyn_flagged = StmState::With(sc_dyn_stm);
                    } else {
                        let mut sc_dyn = Spacecraft::new(
                            OrbitalDynamics::point_masses(init_state, bodies, cosm),
                            spacecraft.dry_mass,
                        );
                        if let Some(fuel_mass) = spacecraft.fuel_mass {
                            sc_dyn.fuel_mass = fuel_mass;
                        }

                        sc_dyn_flagged = StmState::Without(sc_dyn);
                    }
                } else if stm_flag.with() {
                    let mut sc_dyn_stm = Spacecraft::with_stm(
                        OrbitalDynamicsStm::two_body(init_state),
                        spacecraft.dry_mass,
                    );
                    if let Some(fuel_mass) = spacecraft.fuel_mass {
                        sc_dyn_stm.fuel_mass = fuel_mass;
                    }

                    sc_dyn_flagged = StmState::With(sc_dyn_stm);
                } else {
                    let mut sc_dyn =
                        Spacecraft::new(OrbitalDynamics::two_body(init_state), spacecraft.dry_mass);
                    if let Some(fuel_mass) = spacecraft.fuel_mass {
                        sc_dyn.fuel_mass = fuel_mass;
                    }
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
                                    for smdl in amdl.srp.values() {
                                        let mut srp = SolarPressure::default(
                                            smdl.sc_area,
                                            vec![cosm.frame("EME2000"), cosm.frame("Luna")],
                                            &cosm,
                                        );
                                        srp.phi = smdl.phi;
                                        srp.cr = smdl.cr;
                                        match sc_dyn_flagged {
                                            StmState::With(ref mut sc_dyn_stm) => {
                                                sc_dyn_stm.add_model(Box::new(srp));
                                            }
                                            StmState::Without(ref mut sc_dyn) => {
                                                sc_dyn.add_model(Box::new(srp));
                                            }
                                            _ => panic!("should not happen"),
                                        };
                                    }
                                }
                            }
                        }
                    }

                    sc_dyn_flagged = StmState::Without(sc_dyn);
                }

                // Add the acceleration models if applicable
                if let Some(accel_models) = &dynamics.accel_models {
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
                                    let compute_frame = cosm.frame(hmdl.frame.as_str());

                                    if stm_flag.with() {
                                        let hh =
                                            HarmonicsDiff::from_stor(compute_frame, in_mem, &cosm);
                                        match sc_dyn_flagged {
                                            StmState::With(ref mut sc_dyn_stm) => {
                                                sc_dyn_stm.orbital_dyn.add_model(Box::new(hh));
                                            }
                                            _ => panic!("should not happen"),
                                        };
                                    } else {
                                        let hh = Harmonics::from_stor(compute_frame, in_mem, &cosm);
                                        match sc_dyn_flagged {
                                            StmState::Without(ref mut sc_dyn) => {
                                                sc_dyn.orbital_dyn.add_model(Box::new(hh));
                                            }
                                            _ => panic!("should not happen"),
                                        };
                                    }
                                }
                            }
                        }
                    }
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
                let prop_time_s = if prop_event.is_none() {
                    match parse_duration(prop.stop_cond.as_str()) {
                        Ok(duration) => Some(duration.v()),
                        Err(e) => {
                            // Check to see if it's an Epoch
                            match Epoch::from_str(prop.stop_cond.as_str()) {
                                Err(_) => return Err(e),
                                Ok(epoch) => Some(epoch - init_state.dt),
                            }
                        }
                    }
                } else {
                    None
                };

                // Let's see if it's a relative time
                let prop_tol = match prop.tolerance {
                    Some(tol) => tol,
                    None => 1e-12,
                };

                Ok((
                    Self {
                        sc_dyn: sc_dyn_flagged,
                        formatter: None,
                        prop_time_s,
                        prop_event,
                        prop_tol,
                        name: prop_name.clone(),
                    },
                    formatter,
                ))
            }
        }
    }

    pub fn propagator(&mut self) -> Propagator<Spacecraft<'a, OrbitalDynamics<'a>>, RSSStepPV> {
        match self.sc_dyn {
            StmState::Without(ref mut sc_dyn) => {
                let mut p = Propagator::default(sc_dyn, &PropOpts::default());
                p.set_tolerance(self.prop_tol);
                p
            }
            _ => panic!("these dynamics are defined with stm. Use stm_propagator instead"),
        }
    }

    pub fn stm_propagator(
        &mut self,
    ) -> Propagator<Spacecraft<'a, OrbitalDynamicsStm<'a>>, RSSStepPV> {
        match self.sc_dyn {
            StmState::With(ref mut sc_dyn) => {
                let mut p = Propagator::default(sc_dyn, &PropOpts::default());
                p.set_tolerance(self.prop_tol);
                p
            }
            _ => panic!("these dynamics are defined without stm. Use propagator instead"),
        }
    }

    pub fn state(&self) -> SpacecraftState {
        match &self.sc_dyn {
            StmState::With(sc_dyn) => sc_dyn.state(),
            StmState::Without(sc_dyn) => sc_dyn.state(),
            _ => panic!("only call state() on an initialized MDProcess"),
        }
    }

    pub fn execute(mut self) -> Result<(), NyxError> {
        self.execute_with(vec![])
    }

    /// Execute the MD with the provided handlers. Note that you must initialize your own CSV output if that's desired.
    pub fn execute_with(
        &mut self,
        mut hdlrs: Vec<Box<dyn MdHdlr<SpacecraftState>>>,
    ) -> Result<(), NyxError> {
        // Get the prop time before the mutable ref of the propagator
        let maybe_prop_time = self.prop_time_s;
        let maybe_prop_event = self.prop_event.clone();

        // Build the propagator
        let mut prop = self.propagator();
        // Set up the channels
        let (tx, rx) = channel();
        prop.tx_chan = Some(tx);

        let mut initial_state = Some(prop.state());

        info!("Initial state: {}", prop.state());
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
                prop.until_time_elapsed(prop_time)?;
            }
        }

        info!(
            "Final state:   {} (computed in {:.3} seconds)",
            prop.state(),
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
