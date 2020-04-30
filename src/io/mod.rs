extern crate flate2;
extern crate serde_derive;

/// Handles loading of gravity models using files of NASA PDS and GMAT COF. Several gunzipped files are provided with nyx.
pub mod gravity;

/// Handles writing to an XYZV file
pub mod cosmo;

/// Handles reading from frames defined in input files
pub mod frame_serde;

/// Handles reading random variables
pub mod rv;

pub mod scenario;
