extern crate csv;

use crate::dynamics::spacecraft::SpacecraftState;
use crate::io::formatter::StateFormatter;
use std::fs::File;
pub mod ui;

/// A Mission Design handler
pub trait MdHdlr<StateType: Copy>: Send + Sync {
    fn handle(&mut self, state: &StateType);
}

pub struct OrbitStateOutput<'a> {
    csv_out: csv::Writer<File>,
    fmtr: StateFormatter<'a>,
}

impl<'a> OrbitStateOutput<'a> {
    pub fn new(fmtr: StateFormatter<'a>) -> Self {
        let mut wtr = csv::Writer::from_path(fmtr.filename.clone()).expect("could not create file");
        wtr.serialize(&fmtr.headers)
            .expect("could not write headers");
        info!("Saving output to {}", fmtr.filename);

        Self { csv_out: wtr, fmtr }
    }
}

impl<'a> MdHdlr<SpacecraftState> for OrbitStateOutput<'a> {
    fn handle(&mut self, state: &SpacecraftState) {
        self.csv_out
            .serialize(self.fmtr.fmt(&state.orbit))
            .expect("could not format state");
    }
}
