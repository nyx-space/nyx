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

use super::flate2::read::GzDecoder;
use crate::linalg::DMatrix;
use crate::NyxError;
use std::fs::File;
use std::io::prelude::*;
use std::str::FromStr;

/// All gravity potential storage backends must implement this trait in order to be used in the provided dynamics.
/// Currently, only a HashMap based storage is provided. However, the use of this trait enables any application
/// from storing the gravity potential in another way, such as a remote database.
pub trait GravityPotentialStor
where
    Self: Clone + Sized + Sync,
{
    /// Returns the maximum degree of this gravity potential storage (Jn=J2,J3...)
    fn max_degree_n(&self) -> usize;
    /// Returns the maximum order of this gravity potential storage (Jnm=Jn2,Jn3...)
    fn max_order_m(&self) -> usize;
    /// Returns the C_nm and S_nm for the provided order and degree.
    ///
    /// WARNING: It's up to the caller to ensure that no degree or order greater than stored
    /// in this `GravityPotentialStor` is requested. Depending on the implementor, this call might `panic!`.
    fn cs_nm(&self, degree: usize, order: usize) -> (f64, f64);
}

/// `HarmonicsMem` loads the requested gravity potential files and stores them in memory (in a HashMap).
///
/// WARNING: This memory backend may require a lot of RAM (e.g. EMG2008 2190x2190 requires nearly 400 MB of RAM).
#[derive(Clone)]
pub struct HarmonicsMem {
    degree: usize,
    order: usize,
    c_nm: DMatrix<f64>,
    s_nm: DMatrix<f64>,
}

impl HarmonicsMem {
    /// Initialize `HarmonicsMem` with a custom J2 value
    pub fn from_j2(j2: f64) -> HarmonicsMem {
        let mut c_nm = DMatrix::from_element(3, 3, 0.0);
        c_nm[(2, 0)] = j2;

        HarmonicsMem {
            degree: 2 + 1,
            order: 0,
            c_nm,
            s_nm: DMatrix::from_element(3, 3, 0.0),
        }
    }

    /// Initialize `HarmonicsMem` as an EARTH J<sub>2</sub> only using the JGM3 model (available in GMAT)
    ///
    /// Use the embedded Earth parameter. If others are needed, load from `from_shadr` or `from_egm`.
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn j2_jgm3() -> HarmonicsMem {
        Self::from_j2(-4.841_653_748_864_70e-04)
    }

    /// Initialize `HarmonicsMem` as an EARTH J<sub>2</sub> only using the JGM2 model (available in GMAT)
    ///
    /// Use the embedded Earth parameter. If others are needed, load from `from_shadr` or `from_egm`.
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn j2_jgm2() -> HarmonicsMem {
        Self::from_j2(-4.841_653_9e-04)
    }

    /// Initialize `HarmonicsMem` as J<sub>2</sub> only using the EGM2008 model (from the GRACE mission, best model as of 2018)
    ///
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn j2_egm2008() -> HarmonicsMem {
        Self::from_j2(-0.484_165_143_790_815e-03)
    }

    /// Initialize `HarmonicsMem` from the file path (must be a gunzipped file)
    ///
    /// Gravity models provided by `nyx`:
    /// + EMG2008 to 2190 for Earth (tide free)
    /// + Moon to 1500 (from SHADR file)
    /// + Mars to 120 (from SHADR file)
    /// + Venus to 150 (from SHADR file)
    pub fn from_shadr(
        filepath: &str,
        degree: usize,
        order: usize,
        gunzipped: bool,
    ) -> Result<HarmonicsMem, NyxError> {
        Self::load(
            gunzipped, true, //SHADR has a header which we ignore
            degree, order, filepath,
        )
    }

    pub fn from_egm(
        filepath: &str,
        degree: usize,
        order: usize,
        gunzipped: bool,
    ) -> Result<HarmonicsMem, NyxError> {
        Self::load(gunzipped, false, degree, order, filepath)
    }

    pub fn from_cof(
        filepath: &str,
        degree: usize,
        order: usize,
        gunzipped: bool,
    ) -> Result<HarmonicsMem, NyxError> {
        let mut f = File::open(filepath)
            .map_err(|_| NyxError::FileUnreadable(format!("File not found: {}", filepath)))?;
        let mut buffer = vec![0; 0];
        if gunzipped {
            let mut d = GzDecoder::new(f);
            d.read_to_end(&mut buffer).map_err(|_| {
                NyxError::FileUnreadable("could not read file as gunzip".to_string())
            })?;
        } else {
            f.read_to_end(&mut buffer)
                .map_err(|_| NyxError::FileUnreadable("could not read file to end".to_string()))?;
        }

        let data_as_str = String::from_utf8(buffer).map_err(|_| {
            NyxError::FileUnreadable("could not decode file contents as utf8".to_string())
        })?;

        // Since the COF files are so specific, we just code everything up in here.

        let mut c_nm_mat = DMatrix::from_element(degree + 1, degree + 1, 0.0);
        let mut s_nm_mat = DMatrix::from_element(degree + 1, degree + 1, 0.0);
        let mut max_order: usize = 0;
        let mut max_degree: usize = 0;
        for (lno, line) in data_as_str.split('\n').enumerate() {
            if line.is_empty() || !line.starts_with('R') {
                continue; // This is either a comment, a header or "END"
            }
            // These variables need to be declared as mutable because rustc does not know
            // we nwon't match each ino more than once.
            let mut cur_degree: usize = 0;
            let mut cur_order: usize = 0;
            let mut c_nm: f64 = 0.0;
            let mut s_nm: f64 = 0.0;
            for (ino, item) in line.split_whitespace().enumerate() {
                match ino {
                    0 => continue, // We need this so we don't break at every first item
                    1 => match usize::from_str(item) {
                        Ok(val) => cur_degree = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable(format!(
                                "Harmonics file: 
                                could not parse degree `{}` on line {}",
                                item, lno
                            )));
                        }
                    },
                    2 => match usize::from_str(item) {
                        Ok(val) => cur_order = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable(format!(
                                "Harmonics file: 
                                could not parse order `{}` on line {}",
                                item, lno
                            )));
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
                                    return Err(NyxError::FileUnreadable(format!(
                                        "Harmonics file: 
                                        could not parse C_nm `{}` on line {}",
                                        item, lno
                                    )));
                                }
                            }
                        } else {
                            // There is a space as a delimiting character between the C_nm and S_nm only if the S_nm
                            // is a positive number, otherwise, they are continuous (what a great format).
                            if (item.matches('-').count() == 3 && !item.starts_with('-'))
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
                                            return Err(NyxError::FileUnreadable(format!(
                                                "Harmonics file: 
                                                could not parse C_nm `{}` on line {}",
                                                item, lno
                                            )));
                                        }
                                    }
                                    // That mean we have five minus signs, so both the C and S are negative.
                                    let s_nm_str = "-".to_owned() + parts[3] + "-" + parts[4];
                                    match f64::from_str(&s_nm_str) {
                                        Ok(val) => s_nm = val,
                                        Err(_) => {
                                            return Err(NyxError::FileUnreadable(format!(
                                                "Harmonics file: 
                                                could not parse S_nm `{}` on line {}",
                                                item, lno
                                            )));
                                        }
                                    }
                                } else {
                                    // That mean we have fouve minus signs, and since both values are concatenated, C_nm is positive and S_nm is negative
                                    let c_nm_str = parts[0].to_owned() + "-" + parts[1];
                                    match f64::from_str(&c_nm_str) {
                                        Ok(val) => c_nm = val,
                                        Err(_) => {
                                            return Err(NyxError::FileUnreadable(format!(
                                                "Harmonics file: 
                                                could not parse C_nm `{}` on line {}",
                                                item, lno
                                            )));
                                        }
                                    }
                                    // That mean we have five minus signs, so both the C and S are negative.
                                    let s_nm_str = "-".to_owned() + parts[2] + "-" + parts[3];
                                    match f64::from_str(&s_nm_str) {
                                        Ok(val) => s_nm = val,
                                        Err(_) => {
                                            return Err(NyxError::FileUnreadable(format!(
                                                "Harmonics file: 
                                                could not parse S_nm `{}` on line {}",
                                                item, lno
                                            )));
                                        }
                                    }
                                }
                            } else {
                                // We only have the first item, and that's the C_nm
                                match f64::from_str(item) {
                                    Ok(val) => c_nm = val,
                                    Err(_) => {
                                        return Err(NyxError::FileUnreadable(format!(
                                            "Harmonics file: 
                                            could not parse C_nm `{}` on line {}",
                                            item, lno
                                        )));
                                    }
                                }
                            }
                        }
                    }
                    4 => match f64::from_str(item) {
                        // If this exists, then the S_nm is positive.
                        Ok(val) => s_nm = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable(format!(
                                "Harmonics file: 
                                could not parse S_nm `{}` on line {}",
                                item, lno
                            )));
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
                c_nm_mat[(cur_degree, cur_order)] = c_nm;
                s_nm_mat[(cur_degree, cur_order)] = s_nm;
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
        Ok(HarmonicsMem {
            degree: max_degree,
            order: max_order,
            c_nm: c_nm_mat,
            s_nm: s_nm_mat,
        })
    }

    /// `load` handles the actual loading in memory.
    fn load(
        gunzipped: bool,
        skip_first_line: bool,
        degree: usize,
        order: usize,
        filepath: &str,
    ) -> Result<HarmonicsMem, NyxError> {
        let mut f = File::open(filepath)
            .map_err(|_| NyxError::FileUnreadable(format!("File not found: {}", filepath)))?;
        let mut buffer = vec![0; 0];
        if gunzipped {
            let mut d = GzDecoder::new(f);
            d.read_to_end(&mut buffer).map_err(|_| {
                NyxError::FileUnreadable("could not read file as gunzip".to_string())
            })?;
        } else {
            f.read_to_end(&mut buffer)
                .map_err(|_| NyxError::FileUnreadable("could not read file to end".to_string()))?;
        }

        let data_as_str = String::from_utf8(buffer).map_err(|_| {
            NyxError::FileUnreadable("could not decode file contents as utf8".to_string())
        })?;

        let mut c_nm_mat = DMatrix::from_element(degree + 1, degree + 1, 0.0);
        let mut s_nm_mat = DMatrix::from_element(degree + 1, degree + 1, 0.0);

        let mut max_degree: usize = 0;
        let mut max_order: usize = 0;
        for (lno, line) in data_as_str.split('\n').enumerate() {
            if lno == 0 && skip_first_line {
                continue;
            }
            // These variables need to be declared as mutable because rustc does not know
            // we won't match each ino more than once.
            let mut cur_order: usize = 0;
            let mut cur_degree: usize = 0;
            let mut c_nm: f64 = 0.0;
            let mut s_nm: f64 = 0.0;
            for (ino, item) in line.replace(",", " ").split_whitespace().enumerate() {
                match ino {
                    0 => match usize::from_str(item) {
                        Ok(val) => cur_degree = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable(format!(
                                "Harmonics file: 
                                could not parse degree on line {} (`{}`)",
                                lno, item
                            )));
                        }
                    },
                    1 => match usize::from_str(item) {
                        Ok(val) => cur_order = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable(format!(
                                "Harmonics file: 
                                could not parse order on line {} (`{}`)",
                                lno, item
                            )));
                        }
                    },
                    2 => match f64::from_str(&item.replace("D", "E")) {
                        Ok(val) => c_nm = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable(format!(
                                "Harmonics file: 
                                could not parse C_nm `{}` on line {}",
                                item, lno
                            )));
                        }
                    },
                    3 => match f64::from_str(&item.replace("D", "E")) {
                        Ok(val) => s_nm = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable(format!(
                                "Harmonics file: 
                                could not parse S_nm `{}` on line {}",
                                item, lno
                            )));
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
                c_nm_mat[(cur_degree, cur_order)] = c_nm;
                s_nm_mat[(cur_degree, cur_order)] = s_nm;
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
        Ok(HarmonicsMem {
            order: max_order,
            degree: max_degree,
            c_nm: c_nm_mat,
            s_nm: s_nm_mat,
        })
    }
}

impl GravityPotentialStor for HarmonicsMem {
    fn max_order_m(&self) -> usize {
        self.order
    }

    fn max_degree_n(&self) -> usize {
        self.degree
    }

    fn cs_nm(&self, degree: usize, order: usize) -> (f64, f64) {
        (self.c_nm[(degree, order)], self.s_nm[(degree, order)])
    }
}

#[test]
fn test_load_harmonic_files() {
    HarmonicsMem::from_cof("data/JGM3.cof.gz", 50, 50, true).expect("could not load JGM3");

    HarmonicsMem::from_egm("data/EGM2008_to2190_TideFree.gz", 120, 120, true)
        .expect("could not load EGM2008");

    HarmonicsMem::from_shadr("data/Luna_jggrx_1500e_sha.tab.gz", 1500, 1500, true)
        .expect("could not load jggrx");
}
