extern crate regex;
use self::regex::Regex;
use crate::celestia::*;
use crate::dynamics::*;
use crate::io::scenario::{PropagatorKind, ScenarioSerde};
use crate::log::warn;
use crate::od::ui::*;
use crate::propagators::*;
use crate::time::Epoch;
use std::str::FromStr;

pub fn maybe_run(scen: ScenarioSerde, cosm: &Cosm, for_real: bool) {
    if !for_real {
        warn!("DRY-RUN ENABLED");
    }

    // Parse the scenario and create the run at the same time

    for seq_name in &scen.sequence {
        match scen.propagator.get(seq_name) {
            None => panic!("sequence refers to undefined propagator `{}` ", seq_name),
            Some(prop) => {
                let mut orbital_dynamics;
                let mut init_state;
                // Validate the orbital dynamics
                match scen.orbital_dynamics.get(&prop.dynamics) {
                    None => panic!(
                        "propagator `{}` refers to undefined dynamics `{}`",
                        seq_name, prop.dynamics
                    ),
                    Some(dynamics) => match scen.state.get(&dynamics.initial_state) {
                        None => panic!(
                            "dynamics `{}` refers to unknown state `{}`",
                            prop.dynamics, dynamics.initial_state
                        ),
                        Some(init_state_sd) => {
                            // Let's check the frames
                            let state_frame = cosm.frame(init_state_sd.frame.as_str());
                            init_state = init_state_sd.as_state(state_frame);
                            if let Some(integ_frame_name) = &dynamics.integration_frame {
                                let integ_frame = cosm.frame(integ_frame_name);
                                init_state = cosm.frame_chg(&init_state, integ_frame);
                            }

                            // Create the dynamics
                            orbital_dynamics = orbital::OrbitalDynamics::two_body(init_state);
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
                            _ => panic!("unknown duration unit in `{}`", prop.stop_cond),
                        }
                        prop_time_s
                    }
                    None => {
                        // Check to see if it's an Epoch
                        match Epoch::from_str(prop.stop_cond.as_str()) {
                            Err(e) => panic!("{}", e),
                            Ok(epoch) => epoch - init_state.dt,
                        }
                    }
                };

                // Build the propagator
                let mut propagator = match &prop.kind {
                    Some(kind) => match kind {
                        PropagatorKind::Rk89 => {
                            Propagator::default(&mut orbital_dynamics, &PropOpts::default())
                        }
                        _ => panic!("only RK89 is supported for now"),
                    },
                    None => Propagator::default(&mut orbital_dynamics, &PropOpts::default()),
                };

                // And run!

                propagator.until_time_elapsed(prop_time_s);
                println!("{}", orbital_dynamics.state());
            }
        }
    }
}

pub fn dry_run(scenario: ScenarioSerde, cosm: &Cosm) {
    maybe_run(scenario, cosm, false)
}

pub fn run(scenario: ScenarioSerde, cosm: &Cosm) {
    maybe_run(scenario, cosm, true)
}
