extern crate csv;
extern crate regex;

use self::regex::Regex;
pub use crate::celestia::*;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, U6};
pub use crate::dynamics::orbital::{OrbitalDynamics, OrbitalDynamicsStm};
use crate::dynamics::sph_harmonics::{Harmonics, HarmonicsDiff};
pub use crate::dynamics::Dynamics;
use crate::io::formatter::*;
use crate::io::scenario::ScenarioSerde;
use crate::io::ParsingError;
use crate::propagators::{PropOpts, Propagator};
use crate::time::{Epoch, SECONDS_PER_DAY, SECONDS_PER_HOUR, SECONDS_PER_MINUTE};
use std::str::FromStr;
use std::sync::mpsc::channel;
use std::time::Instant;

pub struct MDProcess<'a>
where
    DefaultAllocator: Allocator<f64, U6>,
{
    orbital_dyn: Option<OrbitalDynamics<'a>>,
    orbital_dyn_stm: Option<OrbitalDynamicsStm<'a>>,
    formatter: Option<StateFormatter<'a>>,
    pub output: Vec<State>,
    pub prop_time_s: Option<f64>,
}

impl<'a> MDProcess<'a>
where
    DefaultAllocator: Allocator<f64, U6>,
{
    pub fn try_from_scenario(
        scen: ScenarioSerde,
        cosm: &'a Cosm,
    ) -> Result<Vec<Self>, ParsingError> {
        let mut seq = Vec::with_capacity(10);
        for seq_name in &scen.sequence {
            match scen.propagator.get(seq_name) {
                None => {
                    return Err(ParsingError::MD(format!(
                        "sequence refers to undefined propagator `{}` ",
                        seq_name,
                    )))
                }
                Some(prop) => {
                    // let mut spacecraft_dynamics;
                    let mut orbital_dyn = None;
                    let mut orbital_dyn_stm = None;
                    let mut init_state;
                    // Validate the output
                    let formatter = if let Some(output) = &prop.output {
                        match scen.output.get(output) {
                            None => {
                                return Err(ParsingError::MD(format!(
                                    "propagator `{}` refers to undefined output `{}`",
                                    seq_name, output
                                )))
                            }
                            Some(out) => Some(out.to_state_formatter(cosm)),
                        }
                    } else {
                        None
                    };
                    // Validate the orbital dynamics
                    match scen.orbital_dynamics.get(&prop.dynamics) {
                        None => {
                            return Err(ParsingError::MD(format!(
                                "propagator `{}` refers to undefined dynamics `{}`",
                                seq_name, prop.dynamics
                            )))
                        }
                        Some(dynamics) => match scen.state.get(&dynamics.initial_state) {
                            None => {
                                return Err(ParsingError::MD(format!(
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
                                    if dynamics.with_stm() {
                                        orbital_dyn_stm = Some(OrbitalDynamicsStm::point_masses(
                                            init_state, bodies, cosm,
                                        ));
                                    } else {
                                        orbital_dyn = Some(OrbitalDynamics::point_masses(
                                            init_state, bodies, cosm,
                                        ));
                                    }
                                } else if dynamics.with_stm() {
                                    orbital_dyn_stm =
                                        Some(OrbitalDynamicsStm::two_body(init_state));
                                } else {
                                    orbital_dyn = Some(OrbitalDynamics::two_body(init_state));
                                }

                                // Add the acceleration models if applicable
                                if let Some(accel_models) = &dynamics.accel_models {
                                    for mdl in accel_models {
                                        match scen.accel_models.get(mdl) {
                                            None => {
                                                return Err(ParsingError::MD(format!(
                                                    "dynamics `{}` refers to unknown state `{}`",
                                                    prop.dynamics, dynamics.initial_state
                                                )))
                                            }
                                            Some(amdl) => {
                                                for hmdl in amdl.harmonics.values() {
                                                    let in_mem = hmdl.load();
                                                    let compute_frame =
                                                        cosm.frame(hmdl.frame.as_str());
                                                    if dynamics.with_stm() {
                                                        let hh = HarmonicsDiff::from_stor(
                                                            compute_frame,
                                                            in_mem,
                                                            &cosm,
                                                        );
                                                        let mut new_orbital_dyn_stm =
                                                            orbital_dyn_stm.take().unwrap();
                                                        new_orbital_dyn_stm.add_model(Box::new(hh));
                                                        orbital_dyn_stm = Some(new_orbital_dyn_stm);
                                                    } else {
                                                        let hh = Harmonics::from_stor(
                                                            compute_frame,
                                                            in_mem,
                                                            &cosm,
                                                        );
                                                        let mut new_orbital_dyn =
                                                            orbital_dyn.take().unwrap();
                                                        new_orbital_dyn.add_model(Box::new(hh));
                                                        orbital_dyn = Some(new_orbital_dyn);
                                                    }
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
                                "days" | "day" => prop_time_s *= SECONDS_PER_DAY,
                                "hours" | "hour" => prop_time_s *= SECONDS_PER_HOUR,
                                "min" | "mins" | "minute" | "minutes" => {
                                    prop_time_s *= SECONDS_PER_MINUTE
                                }
                                _ => {
                                    return Err(ParsingError::MD(format!(
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
                                    return Err(ParsingError::MD(format!(
                                        "Could not parse stopping condition: `{}`",
                                        prop.stop_cond
                                    )))
                                }
                                Ok(epoch) => epoch - init_state.dt,
                            }
                        }
                    };

                    let me = Self {
                        orbital_dyn,
                        orbital_dyn_stm,
                        formatter,
                        output: Vec::with_capacity(65_535),
                        prop_time_s: Some(prop_time_s),
                    };

                    seq.push(me);
                }
            }
        }
        Ok(seq)
    }

    pub fn execute(&mut self) {
        // Create the output file
        let mut maybe_wtr = match &self.formatter {
            Some(fmtr) => {
                let mut wtr =
                    csv::Writer::from_path(fmtr.filename.clone()).expect("could not create file");
                wtr.serialize(&fmtr.headers)
                    .expect("could not write headers");
                info!("Saving output to {}", fmtr.filename);
                Some(wtr)
            }
            None => None,
        };

        // Build the propagator
        if let Some(mut orbital_dyn) = self.orbital_dyn.take() {
            let mut prop = Propagator::default(&mut orbital_dyn, &PropOpts::default());
            // Set up the channels
            let (tx, rx) = channel();
            prop.tx_chan = Some(&tx);
            // Run
            let prop_time = self.prop_time_s.unwrap();
            info!(
                "Propagating for {} seconds (~ {:.3} days)",
                prop_time,
                prop_time / SECONDS_PER_DAY
            );
            let start = Instant::now();
            prop.until_time_elapsed(prop_time);
            info!(
                "Done in {:.3} seconds",
                (Instant::now() - start).as_secs_f64()
            );

            while let Ok(prop_state) = rx.try_recv() {
                self.output.push(prop_state);
                if let Some(wtr) = &mut maybe_wtr {
                    wtr.serialize(self.formatter.as_ref().unwrap().format(&prop_state))
                        .expect("could not format state");
                }
            }
        } else if let Some(mut orbital_dyn_stm) = self.orbital_dyn_stm.take() {
            let mut prop = Propagator::default(&mut orbital_dyn_stm, &PropOpts::default());
            // Set up the channels
            let (tx, rx) = channel();
            prop.tx_chan = Some(&tx);
            // Run
            let prop_time = self.prop_time_s.unwrap();
            info!(
                "Propagating for {} seconds (~ {:.3} days)",
                prop_time,
                prop_time / SECONDS_PER_DAY
            );
            let start = Instant::now();
            prop.until_time_elapsed(prop_time);
            info!(
                "Done in {:.3} seconds",
                (Instant::now() - start).as_secs_f64()
            );

            while let Ok(prop_state) = rx.try_recv() {
                self.output.push(prop_state);
                if let Some(wtr) = &mut maybe_wtr {
                    wtr.serialize(self.formatter.as_ref().unwrap().format(&prop_state))
                        .expect("could not format state");
                }
            }
        }
    }
}
