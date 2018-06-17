extern crate csv;

use self::csv::{QuoteStyle, Writer, WriterBuilder};
use celestia::{CoordinateFrame, State};
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

    pub fn append<F: CoordinateFrame>(&mut self, s: State<F>) {
        self.wtr.serialize(s).expect("could not write to XYZV file");
    }
}
