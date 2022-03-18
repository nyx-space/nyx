/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use self::csv::{QuoteStyle, Writer, WriterBuilder};
use crate::cosmic::Orbit;
use std::fs::File;

/// Exports to the XYZV data type used in Cosmographia
pub struct Cosmographia {
    wtr: Writer<File>,
}

impl Cosmographia {
    pub fn from_path(path: String) -> Cosmographia {
        Cosmographia {
            wtr: WriterBuilder::new()
                .delimiter(b' ')
                .quote_style(QuoteStyle::Never)
                .has_headers(false)
                .from_path(path)
                .expect("could not create file"),
        }
    }

    pub fn append(&mut self, s: Orbit) {
        self.wtr.serialize(s).expect("could not write to XYZV file");
    }
}
