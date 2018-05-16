use super::flate2::read::GzDecoder;
use std::io::prelude::*;
use std::fs::File;
use std::str::FromStr;
use std::collections::HashMap;

pub trait GravityPotentialStor
where
    Self: Sized,
{
    /// Returns the maximum order of this gravity potential storage
    fn max_order(&self) -> u16;
    /// Returns the maximum degree of this gravity potential storage
    fn max_degree(&self) -> u16;
    /// Returns the C_nm and S_nm for the provided order and degree.
    ///
    /// WARNING: It's up to the caller to ensure that no degree or order greater than stored
    /// in this `GravityPotentialStor` is requested. Depending on the implementor, this call might `panic!`.
    fn cs_nm(&self, order: u16, degree: u16) -> (f64, f64);
}

/// `MemoryBackend` loads the requested gravity potential files and stores them in memory (in a HashMap).
///
/// WARNING: This memory backend may require a lot of RAM (e.g. EMG2008 2190x2190 requires nearly 400 MB of RAM).
#[derive(Clone)]
pub struct MemoryBackend {
    order: u16,
    degree: u16,
    // data is (degree, order) -> (C_nm, S_nm)
    data: HashMap<(u16, u16), (f64, f64)>,
}

impl MemoryBackend {
    /// Initialize `MemoryBackend` as an EARTH J<sub>2</sub> only using the JGM3 model (available in GMAT)
    ///
    /// Use the embedded Earth parameter. If others are needed, load from `from_shadr` or `from_egm`.
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn j2_jgm3() -> MemoryBackend {
        let mut data = HashMap::new();
        data.insert((0, 0), (0.0, 0.0));
        data.insert((1, 0), (0.0, 0.0));
        data.insert((2, 0), (-4.84165374886470e-04, 0.0));
        MemoryBackend {
            order: 2,
            degree: 0,
            data: data,
        }
    }

    /// Initialize `MemoryBackend` as an EARTH J<sub>2</sub> only using the JGM2 model (available in GMAT)
    ///
    /// Use the embedded Earth parameter. If others are needed, load from `from_shadr` or `from_egm`.
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn j2_jgm2() -> MemoryBackend {
        let mut data = HashMap::new();
        data.insert((0, 0), (0.0, 0.0));
        data.insert((1, 0), (0.0, 0.0));
        data.insert((2, 0), (-4.8416539e-04, 0.0));
        MemoryBackend {
            order: 2,
            degree: 0,
            data: data,
        }
    }

    /// Initialize `MemoryBackend` as J<sub>2</sub> only using the EGM2008 model (from the GRACE mission, best model as of 2018)
    ///
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn j2_egm2008() -> MemoryBackend {
        let mut data = HashMap::new();
        data.insert((0, 0), (0.0, 0.0));
        data.insert((1, 0), (0.0, 0.0));
        data.insert((2, 0), (-0.484165143790815e-03, 0.0));
        MemoryBackend {
            order: 2,
            degree: 0,
            data: data,
        }
    }

    /// Initialize `MemoryBackend` from the file path (must be a gunzipped file)
    ///
    /// Gravity models provided by `nyx`: TODO: add github links and examples
    /// + EMG2008 to 2190 for Earth (tide free)
    /// + Moon to 1500 (from SHADR file)
    /// + Mars to 120 (from SHADR file)
    /// + Venus to 150 (from SHADR file)
    pub fn from_shadr(filepath: String, degree: u16, order: u16, gunzipped: bool) -> MemoryBackend {
        let mut f = File::open(filepath.clone()).expect("could not open file");
        let mut buffer = vec![0; 0];
        if gunzipped {
            let mut d = GzDecoder::new(f);
            d.read_to_end(&mut buffer)
                .expect("could not read the full file");
        } else {
            f.read_to_end(&mut buffer).expect("to end");
        }
        Self::load(
            true, //SHADR has a header which we ignore
            degree,
            order,
            String::from_utf8(buffer).expect("could not decode file contents as utf8"),
            filepath,
        )
    }

    pub fn from_egm(filepath: String, degree: u16, order: u16, gunzipped: bool) -> MemoryBackend {
        let mut f = File::open(filepath.clone()).expect("could not open file");
        let mut buffer = vec![0; 0];
        if gunzipped {
            let mut d = GzDecoder::new(f);
            d.read_to_end(&mut buffer)
                .expect("could not read the full file");
        } else {
            f.read_to_end(&mut buffer).expect("to end");
        }
        Self::load(
            false,
            degree,
            order,
            String::from_utf8(buffer).expect("could not decode file contents as utf8"),
            filepath,
        )
    }

    /// `load` handles the actual loading in memory.
    fn load(
        skip_first_line: bool,
        degree: u16,
        order: u16,
        data_as_str: String,
        filepath: String,
    ) -> MemoryBackend {
        let mut data: HashMap<(u16, u16), (f64, f64)>;
        data = HashMap::new();
        let mut max_order: u16 = 0;
        let mut max_degree: u16 = 0;
        for (lno, line) in data_as_str.split("\n").enumerate() {
            if lno == 0 && skip_first_line {
                continue;
            }
            // These variables need to be declared as mutable because rustc does not know
            // we won't match each ino more than once.
            let mut cur_order: u16 = 0;
            let mut cur_degree: u16 = 0;
            let mut c_nm: f64 = 0.0;
            let mut s_nm: f64 = 0.0;
            for (ino, item) in line.split_whitespace().enumerate() {
                match ino {
                    0 => match u16::from_str(item) {
                        Ok(val) => cur_order = val,
                        Err(_) => {
                            println!("could not parse order on line {} -- ignoring line", lno);
                            break;
                        }
                    },
                    1 => match u16::from_str(item) {
                        Ok(val) => cur_degree = val,
                        Err(_) => {
                            println!("could not parse degree on line {} -- ignoring line", lno);
                            break;
                        }
                    },
                    2 => match f64::from_str(&item.replace("D", "E")) {
                        Ok(val) => c_nm = val,
                        Err(_) => {
                            println!(
                                "could not parse C_nm `{}` on line {} -- ignoring line",
                                item, lno
                            );
                            break;
                        }
                    },
                    3 => match f64::from_str(&item.replace("D", "E")) {
                        Ok(val) => s_nm = val,
                        Err(_) => {
                            println!("could not parse S_nm on line {} -- ignoring line", lno);
                            break;
                        }
                    },
                    _ => break, // We aren't storing the covariance of these harmonics
                }
            }

            if cur_order > order {
                // The file is organized by order, so once we've passed the maximum order we want,
                // we can safely stop reading the file.
                break;
            }

            // Only insert this data into the hashmap if it's within the required degree as well
            if cur_degree <= degree {
                data.insert((cur_order, cur_degree), (c_nm, s_nm));
            }
            // This serves as a warning.
            max_order = if order > max_order { order } else { max_order };
            max_degree = if degree > max_degree {
                degree
            } else {
                max_degree
            };
        }
        if max_degree < degree || max_order < order {
            warn!(
                "{} only contained (degree, order) of ({}, {}) instead of requested ({}, {})",
                filepath, max_degree, max_order, degree, order
            );
        } else {
            info!(
                "{} loaded with (degree, order) = ({}, {})",
                filepath, degree, order
            );
        }
        MemoryBackend {
            order: max_order,
            degree: max_degree,
            data: data,
        }
    }
}

impl GravityPotentialStor for MemoryBackend {
    fn max_order(&self) -> u16 {
        self.order
    }

    fn max_degree(&self) -> u16 {
        self.degree
    }

    fn cs_nm(&self, order: u16, degree: u16) -> (f64, f64) {
        let &(c, s) = self.data.get(&(order, degree)).unwrap();
        (c, s)
    }
}
