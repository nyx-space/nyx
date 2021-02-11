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
