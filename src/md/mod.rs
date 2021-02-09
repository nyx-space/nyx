extern crate bacon_sci;
extern crate csv;

use crate::io::formatter::StateFormatter;
use crate::{Orbit, SpacecraftState};
use std::fs::File;
pub mod ui;

pub mod trajectory;

pub type ScTraj = trajectory::Traj<SpacecraftState>;
pub type Ephemeris = trajectory::Traj<Orbit>;

/// A Mission Design handler
pub trait MdHdlr<StateType: Copy>: Send + Sync {
    fn handle(&mut self, state: &StateType);
}

pub struct OrbitStateOutput {
    csv_out: csv::Writer<File>,
    fmtr: StateFormatter,
}

impl OrbitStateOutput {
    pub fn new(fmtr: StateFormatter) -> Self {
        let mut wtr = csv::Writer::from_path(fmtr.filename.clone()).expect("could not create file");
        wtr.serialize(&fmtr.headers)
            .expect("could not write headers");
        info!("Saving output to {}", fmtr.filename);

        Self { csv_out: wtr, fmtr }
    }
}

impl MdHdlr<SpacecraftState> for OrbitStateOutput {
    fn handle(&mut self, state: &SpacecraftState) {
        self.csv_out
            .serialize(self.fmtr.fmt(&state.orbit))
            .expect("could not format state");
    }
}
