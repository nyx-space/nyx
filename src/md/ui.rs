extern crate csv;
extern crate regex;

use self::regex::Regex;
pub use crate::celestia::*;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, U6};
pub use crate::dynamics::orbital::{OrbitalDynamics, OrbitalDynamicsStm};
use crate::dynamics::spacecraft::Spacecraft;
use crate::dynamics::sph_harmonics::Harmonics;
pub use crate::dynamics::Dynamics;
use crate::io::formatter::*;
use crate::io::scenario::ScenarioSerde;
use crate::io::ParsingError;
use crate::propagators::{PropOpts, Propagator};
use crate::time::{Epoch, SECONDS_PER_DAY, SECONDS_PER_HOUR, SECONDS_PER_MINUTE};
use std::str::FromStr;
use std::sync::mpsc::channel;
use std::time::Instant;

/// An MDProcess allows the creation and propagation of a spacecraft subjected to some dynamics
pub struct MDProcess<'a>
where
    DefaultAllocator: Allocator<f64, U6>,
{
    sc_dyn: Spacecraft<'a, OrbitalDynamics<'a>>,
    formatter: Option<StateFormatter<'a>>,
    pub output: Vec<State>,
    pub prop_time_s: Option<f64>,
    pub name: String,
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
                    let mut sc_dyn;
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
                    match scen.spacecraft.get(&prop.dynamics) {
                        None => {
                            return Err(ParsingError::MD(format!(
                                "propagator `{}` refers to undefined spacecraft `{}`",
                                seq_name, prop.dynamics
                            )))
                        }
                        Some(spacecraft) =>
                        // Validate the orbital dynamics
                        {
                            match scen.orbital_dynamics.get(&spacecraft.orbital_dynamics) {
                                None => {
                                    return Err(ParsingError::MD(format!(
                                        "spacecraft `{}` refers to undefined dynamics `{}`",
                                        &prop.dynamics, spacecraft.orbital_dynamics
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
                                        if let Some(integ_frame_name) = &dynamics.integration_frame
                                        {
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
                                                        match cosm.try_frame(
                                                            format!("{} j2000", obj).as_str(),
                                                        ) {
                                                            Ok(frame) => {
                                                                bodies.push(frame.exb_id())
                                                            }
                                                            Err(_) => {
                                                                bodies.push(
                                                                    cosm.frame(
                                                                        format!(
                                                                            "{} barycenter j2000",
                                                                            obj
                                                                        )
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
                                            if let Some(pos) = bodies
                                                .iter()
                                                .position(|x| *x == state_frame.exb_id())
                                            {
                                                bodies.remove(pos);
                                            }
                                            if dynamics.with_stm() {
                                                return Err(ParsingError::MD(format!(
                                                    "dynamics `{}` require STM, but this is a mission design scenario",
                                                    prop.dynamics
                                                )));
                                            } else {
                                                sc_dyn = Spacecraft::new(
                                                    OrbitalDynamics::point_masses(
                                                        init_state, bodies, cosm,
                                                    ),
                                                    spacecraft.dry_mass,
                                                );
                                                if let Some(fuel_mass) = spacecraft.fuel_mass {
                                                    sc_dyn.fuel_mass = fuel_mass;
                                                }
                                            }
                                        } else if dynamics.with_stm() {
                                            return Err(ParsingError::MD(format!(
                                                "dynamics `{}` require STM, but this is a mission design scenario",
                                                prop.dynamics
                                            )));
                                        } else {
                                            sc_dyn = Spacecraft::new(
                                                OrbitalDynamics::two_body(init_state),
                                                spacecraft.dry_mass,
                                            );
                                            if let Some(fuel_mass) = spacecraft.fuel_mass {
                                                sc_dyn.fuel_mass = fuel_mass;
                                            }
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

                                                            let hh = Harmonics::from_stor(
                                                                compute_frame,
                                                                in_mem,
                                                                &cosm,
                                                            );
                                                            sc_dyn
                                                                .orbital_dyn
                                                                .add_model(Box::new(hh));
                                                        }
                                                    }
                                                }
                                            }
                                        }
                                    }
                                },
                            }
                        }
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
                        sc_dyn,
                        formatter,
                        output: Vec::with_capacity(65_535),
                        prop_time_s: Some(prop_time_s),
                        name: seq_name.clone(),
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
        let mut prop = Propagator::default(&mut self.sc_dyn, &PropOpts::default());
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
            self.output.push(prop_state.orbit);
            if let Some(wtr) = &mut maybe_wtr {
                wtr.serialize(self.formatter.as_ref().unwrap().fmt(&prop_state.orbit))
                    .expect("could not format state");
            }
        }
    }
}
