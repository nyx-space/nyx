use super::flate2::read::GzDecoder;
use std::collections::HashMap;
use std::fs::File;
use std::io::prelude::*;
use std::str::FromStr;

/// All gravity potential storage backends must implement this trait in order to be used in the provided dynamics.
/// Currently, only a HashMap based storage is provided. However, the use of this trait enables any application
/// from storing the gravity potential in another way, such as a remote database.
pub trait GravityPotentialStor
where
    Self: Sized,
{
    /// Returns the maximum degree of this gravity potential storage (Jn=J2,J3...)
    fn max_degree(&self) -> u16;
    /// Returns the maximum order of this gravity potential storage (Jnm=Jn2,Jn3...)
    fn max_order(&self) -> u16;
    /// Returns the C_nm and S_nm for the provided order and degree.
    ///
    /// WARNING: It's up to the caller to ensure that no degree or order greater than stored
    /// in this `GravityPotentialStor` is requested. Depending on the implementor, this call might `panic!`.
    fn cs_nm(&self, degree: u16, order: u16) -> (f64, f64);
}

/// `MemoryBackend` loads the requested gravity potential files and stores them in memory (in a HashMap).
///
/// WARNING: This memory backend may require a lot of RAM (e.g. EMG2008 2190x2190 requires nearly 400 MB of RAM).
#[derive(Clone)]
pub struct MemoryBackend {
    degree: u16,
    order: u16,
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
        data.insert((1, 1), (0.0, 0.0));
        data.insert((2, 0), (-4.841_653_748_864_70e-04, 0.0));
        data.insert((2, 1), (0.0, 0.0));
        data.insert((2, 2), (0.0, 0.0));

        MemoryBackend {
            degree: 2,
            order: 0,
            data,
        }
    }

    /// Initialize `MemoryBackend` as an EARTH J<sub>2</sub> only using the JGM2 model (available in GMAT)
    ///
    /// Use the embedded Earth parameter. If others are needed, load from `from_shadr` or `from_egm`.
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn j2_jgm2() -> MemoryBackend {
        let mut data = HashMap::new();
        data.insert((1, 0), (0.0, 0.0));
        data.insert((2, 0), (-4.841_653_9e-04, 0.0));
        MemoryBackend {
            degree: 2,
            order: 0,
            data,
        }
    }

    /// Initialize `MemoryBackend` as J<sub>2</sub> only using the EGM2008 model (from the GRACE mission, best model as of 2018)
    ///
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn j2_egm2008() -> MemoryBackend {
        let mut data = HashMap::new();
        data.insert((1, 0), (0.0, 0.0));
        data.insert((2, 0), (-0.484_165_143_790_815e-03, 0.0));
        MemoryBackend {
            degree: 2,
            order: 0,
            data,
        }
    }

    /// Initialize `MemoryBackend` from the file path (must be a gunzipped file)
    ///
    /// Gravity models provided by `nyx`: TODO: add github links and examples
    /// + EMG2008 to 2190 for Earth (tide free)
    /// + Moon to 1500 (from SHADR file)
    /// + Mars to 120 (from SHADR file)
    /// + Venus to 150 (from SHADR file)
    pub fn from_shadr(filepath: &str, degree: u16, order: u16, gunzipped: bool) -> MemoryBackend {
        let mut f = File::open(filepath).expect("could not open file");
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

    pub fn from_egm(filepath: &str, degree: u16, order: u16, gunzipped: bool) -> MemoryBackend {
        let mut f = File::open(filepath).expect("could not open file");
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

    pub fn from_cof(filepath: &str, degree: u16, order: u16, gunzipped: bool) -> MemoryBackend {
        let mut f = File::open(filepath).expect("could not open file");
        let mut buffer = vec![0; 0];
        if gunzipped {
            let mut d = GzDecoder::new(f);
            d.read_to_end(&mut buffer)
                .expect("could not read the full file");
        } else {
            f.read_to_end(&mut buffer).expect("to end");
        }

        // Since the COF files are so specific, we just code everything up in here.

        let mut data: HashMap<(u16, u16), (f64, f64)>;
        data = HashMap::new();
        // Immediately add data which will be requested but may not exist (will be overwritten if it does)
        data.insert((1, 0), (0.0, 0.0));
        data.insert((1, 1), (0.0, 0.0));
        let mut max_order: u16 = 0;
        let mut max_degree: u16 = 0;
        for (lno, line) in String::from_utf8(buffer)
            .expect("error decoding utf8")
            .split('\n')
            .enumerate()
        {
            if line.is_empty() || line.chars().next().unwrap() != 'R' {
                continue; // This is either a comment, a header or "END"
            }
            // These variables need to be declared as mutable because rustc does not know
            // we nwon't match each ino more than once.
            let mut cur_degree: u16 = 0;
            let mut cur_order: u16 = 0;
            let mut c_nm: f64 = 0.0;
            let mut s_nm: f64 = 0.0;
            for (ino, item) in line.split_whitespace().enumerate() {
                match ino {
                    0 => continue, // We need this so we don't break at every first item
                    1 => match u16::from_str(item) {
                        Ok(val) => cur_degree = val,
                        Err(_) => {
                            println!("could not parse degree on line {} -- ignoring line", lno);
                            break;
                        }
                    },
                    2 => match u16::from_str(item) {
                        Ok(val) => cur_order = val,
                        Err(_) => {
                            println!("could not parse order on line {} -- ignoring line", lno);
                            break;
                        }
                    },
                    3 => {
                        // If we are at degree zero, then there is only one item, so we can parse that and
                        // set the S_nm to zero.
                        if degree == 0 {
                            s_nm = 0.0;
                            match f64::from_str(item) {
                                Ok(val) => c_nm = val,
                                Err(_) => {
                                    println!(
                                        "could not parse C_nm `{}` on line {} -- ignoring line",
                                        item, lno
                                    );
                                    break;
                                }
                            }
                        } else {
                            // There is a space as a delimiting character between the C_nm and S_nm only if the S_nm
                            // is a positive number, otherwise, they are continuous (what a great format).
                            if (item.matches('-').count() == 3
                                && item.chars().next().unwrap() != '-')
                                || item.matches('-').count() == 4
                            {
                                // Now we have two items concatenated into one... great
                                let parts: Vec<&str> = item.split('-').collect();
                                if parts.len() == 5 {
                                    // That mean we have five minus signs, so both the C and S are negative.
                                    let c_nm_str = "-".to_owned() + parts[1] + "-" + parts[2];
                                    match f64::from_str(&c_nm_str) {
                                        Ok(val) => c_nm = val,
                                        Err(_) => {
                                            println!("could not parse C_nm `{}` on line {} -- ignoring line", item, lno);
                                            break;
                                        }
                                    }
                                    // That mean we have five minus signs, so both the C and S are negative.
                                    let s_nm_str = "-".to_owned() + parts[3] + "-" + parts[4];
                                    match f64::from_str(&s_nm_str) {
                                        Ok(val) => s_nm = val,
                                        Err(_) => {
                                            println!("could not parse S_nm `{}` on line {} -- ignoring line", item, lno);
                                            break;
                                        }
                                    }
                                } else {
                                    // That mean we have fouve minus signs, and since both values are concatenated, C_nm is positive and S_nm is negative
                                    let c_nm_str = parts[0].to_owned() + "-" + parts[1];
                                    match f64::from_str(&c_nm_str) {
                                        Ok(val) => c_nm = val,
                                        Err(_) => {
                                            println!("could not parse C_nm `{}` on line {} -- ignoring line", item, lno);
                                            break;
                                        }
                                    }
                                    // That mean we have five minus signs, so both the C and S are negative.
                                    let s_nm_str = "-".to_owned() + parts[2] + "-" + parts[3];
                                    match f64::from_str(&s_nm_str) {
                                        Ok(val) => s_nm = val,
                                        Err(_) => {
                                            println!("could not parse S_nm `{}` on line {} -- ignoring line", item, lno);
                                            break;
                                        }
                                    }
                                }
                            } else {
                                // We only have the first item, and that's the C_nm
                                match f64::from_str(item) {
                                    Ok(val) => c_nm = val,
                                    Err(_) => {
                                        println!(
                                            "could not parse C_nm `{}` on line {} -- ignoring line",
                                            item, lno
                                        );
                                        break;
                                    }
                                }
                            }
                        }
                    }
                    4 => match f64::from_str(item) {
                        // If this exists, then the S_nm is positive.
                        Ok(val) => s_nm = val,
                        Err(_) => {
                            println!("could not parse S_nm on line {} -- ignoring line", lno);
                            break;
                        }
                    },
                    _ => break, // We aren't storing the covariance of these harmonics
                }
            }

            if cur_degree > degree {
                // The file is organized by degree, so once we've passed the maximum degree we want,
                // we can safely stop reading the file.
                break;
            }

            // Only insert this data into the hashmap if it's within the required order as well
            if cur_order <= order {
                data.insert((cur_degree, cur_order), (c_nm, s_nm));
            }
            if max_degree == 0 {
                // Let's populate with zeros.
                for n in 0..degree {
                    for m in 0..order {
                        data.insert((n, m), (0.0, 0.0));
                    }
                }
            }

            // This serves as a warning.
            max_order = if cur_order > max_order {
                cur_order
            } else {
                max_order
            };
            max_degree = if cur_degree > max_degree {
                cur_degree
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
            degree: max_degree,
            order: max_order,
            data,
        }
    }

    /// `load` handles the actual loading in memory.
    fn load(
        skip_first_line: bool,
        degree: u16,
        order: u16,
        data_as_str: String,
        filepath: &str,
    ) -> MemoryBackend {
        let mut data: HashMap<(u16, u16), (f64, f64)>;
        data = HashMap::new();
        // Immediately add data which will be requested but may not exist (will be overwritten if it does)
        data.insert((1, 0), (0.0, 0.0));
        data.insert((1, 1), (0.0, 0.0));
        let mut max_degree: u16 = 0;
        let mut max_order: u16 = 0;
        for (lno, line) in data_as_str.split('\n').enumerate() {
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
                        Ok(val) => cur_degree = val,
                        Err(_) => {
                            println!("could not parse degree on line {} -- ignoring line", lno);
                            break;
                        }
                    },
                    1 => match u16::from_str(item) {
                        Ok(val) => cur_order = val,
                        Err(_) => {
                            println!("could not parse order on line {} -- ignoring line", lno);
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

            if cur_degree > degree {
                // The file is organized by degree, so once we've passed the maximum degree we want,
                // we can safely stop reading the file.
                break;
            }

            // Only insert this data into the hashmap if it's within the required order as well
            if cur_order <= order {
                data.insert((cur_degree, cur_order), (c_nm, s_nm));
            }
            // This serves as a warning.
            max_order = if cur_order > max_order {
                cur_order
            } else {
                max_order
            };
            max_degree = if cur_degree > max_degree {
                cur_degree
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
            data,
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

    fn cs_nm(&self, degree: u16, order: u16) -> (f64, f64) {
        let &(c, s) = &self.data[&(degree, order)];
        (c, s)
    }
}
