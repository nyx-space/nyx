extern crate regex;

use self::regex::Regex;
pub use crate::celestia::*;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, U6};
pub use crate::dynamics::orbital::OrbitalDynamics;
use crate::dynamics::sph_harmonics::Harmonics;
use crate::io::scenario::ScenarioSerde;
use crate::propagators::{PropOpts, Propagator};
use crate::time::Epoch;
use std::str::FromStr;
use std::sync::mpsc::channel;

#[derive(Debug)]
pub enum MDError {
    ParsingError(String),
}

pub struct MDProcess<'a>
where
    DefaultAllocator: Allocator<f64, U6>,
{
    /// Propagator used for the mission design
    // prop: &'a mut Propagator<'a, OrbitalDynamics<'a>, RSSStepPV>,
    orbit_dyn: OrbitalDynamics<'a>,
    /// Vector of estimates available after a pass
    pub output: Vec<State>,
    pub prop_time_s: Option<f64>,
}

impl<'a> MDProcess<'a>
where
    DefaultAllocator: Allocator<f64, U6>,
{
    pub fn new(orbit_dyn: OrbitalDynamics<'a>) -> Self {
        Self {
            orbit_dyn,
            output: Vec::with_capacity(65_535),
            prop_time_s: None,
        }
    }

    pub fn try_from_scenario(scen: ScenarioSerde, cosm: &'a Cosm) -> Result<Vec<Self>, MDError> {
        let mut seq = Vec::with_capacity(10);
        for seq_name in &scen.sequence {
            match scen.propagator.get(seq_name) {
                None => {
                    return Err(MDError::ParsingError(format!(
                        "sequence refers to undefined propagator `{}` ",
                        seq_name,
                    )))
                }
                Some(prop) => {
                    // let mut spacecraft_dynamics;
                    let mut orbital_dynamics;
                    let mut init_state;
                    // Validate the orbital dynamics
                    match scen.orbital_dynamics.get(&prop.dynamics) {
                        None => {
                            return Err(MDError::ParsingError(format!(
                                "propagator `{}` refers to undefined dynamics `{}`",
                                seq_name, prop.dynamics
                            )))
                        }
                        Some(dynamics) => match scen.state.get(&dynamics.initial_state) {
                            None => {
                                return Err(MDError::ParsingError(format!(
                                    "dynamics `{}` refers to unknown state `{}`",
                                    prop.dynamics, dynamics.initial_state
                                )))
                            }
                            Some(init_state_sd) => {
                                // Let's check the frames
                                let state_frame = cosm.frame(init_state_sd.frame.as_str());
                                init_state = init_state_sd.as_state(state_frame);
                                if let Some(integ_frame_name) = &dynamics.integration_frame {
                                    let integ_frame = cosm.frame(integ_frame_name);
                                    init_state = cosm.frame_chg(&init_state, integ_frame);
                                }
                                // Create the dynamics
                                if let Some(pts_masses) = &dynamics.point_masses {
                                    // Get the object IDs from name
                                    let mut bodies = Vec::with_capacity(10);
                                    for obj in pts_masses {
                                        match cosm.try_frame(obj) {
                                            Ok(frame) => bodies.push(frame.exb_id()),
                                            Err(_) => {
                                                // Let's try with "j2000" appended
                                                match cosm
                                                    .try_frame(format!("{} j2000", obj).as_str())
                                                {
                                                    Ok(frame) => bodies.push(frame.exb_id()),
                                                    Err(_) => {
                                                        bodies.push(
                                                            cosm.frame(
                                                                format!("{} barycenter j2000", obj)
                                                                    .as_str(),
                                                            )
                                                            .exb_id(),
                                                        );
                                                    }
                                                }
                                            }
                                        }
                                    }
                                    // Remove bodies which are part of the state
                                    if let Some(pos) =
                                        bodies.iter().position(|x| *x == state_frame.exb_id())
                                    {
                                        bodies.remove(pos);
                                    }
                                    orbital_dynamics =
                                        OrbitalDynamics::point_masses(init_state, bodies, cosm);
                                } else {
                                    orbital_dynamics = OrbitalDynamics::two_body(init_state);
                                }

                                // Add the acceleration models if applicable
                                if let Some(accel_models) = &dynamics.accel_models {
                                    for mdl in accel_models {
                                        match scen.accel_models.get(mdl) {
                                            None => {
                                                return Err(MDError::ParsingError(format!(
                                                    "dynamics `{}` refers to unknown state `{}`",
                                                    prop.dynamics, dynamics.initial_state
                                                )))
                                            }
                                            Some(amdl) => {
                                                for hmdl in amdl.harmonics.values() {
                                                    let in_mem = hmdl.load();
                                                    let compute_frame =
                                                        cosm.frame(hmdl.frame.as_str());
                                                    let hh = Harmonics::from_stor(
                                                        compute_frame,
                                                        in_mem,
                                                        &cosm,
                                                    );
                                                    orbital_dynamics.add_model(Box::new(hh));
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        },
                    }
                    // Validate the stopping condition
                    // Let's see if it's a relative time
                    let reg = Regex::new(r"^(\d+\.?\d*)\W*(\w+)$").unwrap();
                    let prop_time_s = match reg.captures(prop.stop_cond.as_str()) {
                        Some(cap) => {
                            let mut prop_time_s = cap[1].to_owned().parse::<f64>().unwrap();
                            match cap[2].to_owned().to_lowercase().as_str() {
                                "days" | "day" => prop_time_s *= 84_600.0,
                                "hours" | "hour" => prop_time_s *= 3_600.0,
                                "min" | "mins" | "minute" | "minutes" => prop_time_s *= 60.0,
                                _ => {
                                    return Err(MDError::ParsingError(format!(
                                        "unknown duration unit in `{}`",
                                        prop.stop_cond
                                    )))
                                }
                            }
                            prop_time_s
                        }
                        None => {
                            // Check to see if it's an Epoch
                            match Epoch::from_str(prop.stop_cond.as_str()) {
                                Err(_) => {
                                    return Err(MDError::ParsingError(format!(
                                        "Could not parse stopping condition: `{}`",
                                        prop.stop_cond
                                    )))
                                }
                                Ok(epoch) => epoch - init_state.dt,
                            }
                        }
                    };
                    let mut me = Self::new(orbital_dynamics);
                    me.prop_time_s = Some(prop_time_s);
                    seq.push(me);
                }
            }
        }
        Ok(seq)
    }

    pub fn execute(&mut self) {
        // Build the propagator
        let mut prop = Propagator::default(&mut self.orbit_dyn, &PropOpts::default());
        // Set up the channels
        let (tx, rx) = channel();
        prop.tx_chan = Some(&tx);
        // Run
        prop.until_time_elapsed(self.prop_time_s.unwrap());

        while let Ok(prop_state) = rx.try_recv() {
            self.output.push(prop_state);
        }
    }
}
