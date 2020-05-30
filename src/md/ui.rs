extern crate csv;

pub use crate::celestia::*;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, U6};
pub use crate::dynamics::orbital::{OrbitalDynamics, OrbitalDynamicsStm, OrbitalDynamicsT};
use crate::dynamics::spacecraft::{Spacecraft, SpacecraftState};
use crate::dynamics::sph_harmonics::{Harmonics, HarmonicsDiff};
pub use crate::dynamics::Dynamics;
use crate::io::formatter::*;
use crate::io::scenario::ScenarioSerde;
use crate::io::{parse_duration, ParsingError};
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
        match self {
            StmState::With(_) => true,
            _ => false,
        }
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
    pub output: Vec<State>,
    pub prop_time_s: Option<f64>,
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
    ) -> Result<Self, ParsingError> {
        match scen.propagator.get(&prop_name) {
            None => Err(ParsingError::PropagatorNotFound(prop_name)),
            Some(prop) => {
                #[allow(unused_assignments)]
                let mut sc_dyn_flagged = StmState::Unknown;
                let mut init_state;
                // Validate the output
                let formatter = if let Some(output) = &prop.output {
                    match scen.output.get(output) {
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
                match scen.spacecraft.get(&prop.dynamics) {
                    None => {
                        return Err(ParsingError::MD(format!(
                            "propagator `{}` refers to undefined spacecraft `{}`",
                            prop_name, prop.dynamics
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
                                                    match cosm.try_frame(
                                                        format!("{} j2000", obj).as_str(),
                                                    ) {
                                                        Ok(frame) => bodies.push(frame.exb_id()),
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
                                        if let Some(pos) =
                                            bodies.iter().position(|x| *x == state_frame.exb_id())
                                        {
                                            bodies.remove(pos);
                                        }
                                        if stm_flag.with() {
                                            let mut sc_dyn_stm = Spacecraft::with_stm(
                                                OrbitalDynamicsStm::point_masses(
                                                    init_state, bodies, cosm,
                                                ),
                                                spacecraft.dry_mass,
                                            );
                                            if let Some(fuel_mass) = spacecraft.fuel_mass {
                                                sc_dyn_stm.fuel_mass = fuel_mass;
                                            }

                                            sc_dyn_flagged = StmState::With(sc_dyn_stm);
                                        } else {
                                            let mut sc_dyn = Spacecraft::new(
                                                OrbitalDynamics::point_masses(
                                                    init_state, bodies, cosm,
                                                ),
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
                                        let mut sc_dyn = Spacecraft::new(
                                            OrbitalDynamics::two_body(init_state),
                                            spacecraft.dry_mass,
                                        );
                                        if let Some(fuel_mass) = spacecraft.fuel_mass {
                                            sc_dyn.fuel_mass = fuel_mass;
                                        }

                                        sc_dyn_flagged = StmState::Without(sc_dyn);
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

                                                        if stm_flag.with() {
                                                            let hh = HarmonicsDiff::from_stor(
                                                                compute_frame,
                                                                in_mem,
                                                                &cosm,
                                                            );
                                                            match sc_dyn_flagged {
                                                                StmState::With(
                                                                    ref mut sc_dyn_stm,
                                                                ) => {
                                                                    sc_dyn_stm
                                                                        .orbital_dyn
                                                                        .add_model(Box::new(hh));
                                                                }
                                                                _ => panic!("should not happen"),
                                                            };
                                                        } else {
                                                            let hh = Harmonics::from_stor(
                                                                compute_frame,
                                                                in_mem,
                                                                &cosm,
                                                            );
                                                            match sc_dyn_flagged {
                                                                StmState::Without(
                                                                    ref mut sc_dyn,
                                                                ) => {
                                                                    sc_dyn
                                                                        .orbital_dyn
                                                                        .add_model(Box::new(hh));
                                                                }
                                                                _ => panic!("should not happen"),
                                                            };
                                                        }
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
                let prop_time_s = match parse_duration(prop.stop_cond.as_str()) {
                    Ok(duration) => duration,
                    Err(e) => {
                        // Check to see if it's an Epoch
                        match Epoch::from_str(prop.stop_cond.as_str()) {
                            Err(_) => return Err(e),
                            Ok(epoch) => epoch - init_state.dt,
                        }
                    }
                };

                Ok(Self {
                    sc_dyn: sc_dyn_flagged,
                    formatter,
                    output: Vec::with_capacity(65_535),
                    prop_time_s: Some(prop_time_s),
                    name: prop_name.clone(),
                })
            }
        }
    }

    pub fn propagator(&mut self) -> Propagator<Spacecraft<'a, OrbitalDynamics<'a>>, RSSStepPV> {
        match self.sc_dyn {
            StmState::Without(ref mut sc_dyn) => Propagator::default(sc_dyn, &PropOpts::default()),
            _ => panic!("these dynamics are defined with stm. Use stm_propagator instead"),
        }
    }

    pub fn stm_propagator(
        &mut self,
    ) -> Propagator<Spacecraft<'a, OrbitalDynamicsStm<'a>>, RSSStepPV> {
        match self.sc_dyn {
            StmState::With(ref mut sc_dyn) => Propagator::default(sc_dyn, &PropOpts::default()),
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

        // Get the prop time before the mutable ref
        let prop_time = self.prop_time_s.unwrap();

        // Build the propagator
        let mut prop = self.propagator();
        // Set up the channels
        let (tx, rx) = channel();
        prop.tx_chan = Some(tx);
        // Run
        info!(
            "Propagating for {} seconds (~ {:.3} days)",
            prop_time,
            prop_time / SECONDS_PER_DAY
        );
        info!("Initial state: {}", prop.state());
        let start = Instant::now();
        prop.until_time_elapsed(prop_time);
        info!(
            "Final state:   {} (computed in {:.3} seconds)",
            prop.state(),
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
