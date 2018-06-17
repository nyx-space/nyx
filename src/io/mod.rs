extern crate flate2;

/// Handles loading of gravity models using files of NASA PDS and GMAT COF. Several gunzipped files are provided with nyx.
pub mod gravity;

/// Handles writing to an XYZV file
pub mod cosmo;
