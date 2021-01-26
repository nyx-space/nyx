extern crate csv;

use self::csv::{QuoteStyle, Writer, WriterBuilder};
use crate::celestia::Orbit;
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
