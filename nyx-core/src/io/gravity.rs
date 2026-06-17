/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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

use crate::NyxError;
use crate::linalg::DMatrix;
use anise::errors::AlmanacError;
use anise::frames::{Frame, FrameUid};
use anise::prelude::Almanac;
use flate2::read::GzDecoder;
use log::{info, warn};
use serde::{Deserialize, Serialize};
use serde_dhall::{SimpleType, StaticType};
use std::collections::HashMap;
use std::fmt::Debug;
use std::fs::File;
use std::io::prelude::*;
use std::path::{Path, PathBuf};
use std::str::FromStr;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Configuration holder for gravity field.
///
/// Data is first loaded as a SHADR, if that fails, Nyx will try to load it as a COF file.
#[derive(Clone, Serialize, Deserialize, Debug)]
#[cfg_attr(feature = "python", pyclass(from_py_object, get_all, set_all))]
pub struct GravityFieldConfig {
    /// Desired degree
    pub degree: usize,
    /// Desired order
    pub order: usize,
    /// Path to the file, relative to the current working directory
    pub filepath: PathBuf,
    /// Set to true if the data is gunzipped
    pub gunzipped: bool,
    /// The frame in which to compute this gravity field
    pub frame: FrameUid,
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl GravityFieldConfig {
    #[pyo3(signature=(degree, order, filepath, frame, gunzipped=true))]
    #[new]
    fn py_new(
        degree: usize,
        order: usize,
        filepath: PathBuf,
        frame: FrameUid,
        gunzipped: bool,
    ) -> Self {
        Self {
            filepath,
            gunzipped,
            degree,
            order,
            frame,
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}

/// `GravityFieldData` loads the requested gravity potential files and stores them in memory (in a HashMap).
///
/// WARNING: This memory backend may require a lot of RAM (e.g. EMG2008 2190x2190 requires nearly 400 MB of RAM).
#[derive(Clone)]
pub struct GravityFieldData {
    degree: usize,
    order: usize,
    c_nm: DMatrix<f64>,
    s_nm: DMatrix<f64>,
    pub frame: Frame,
}

impl GravityFieldData {
    pub fn from_config(cfg: GravityFieldConfig, almanac: &Almanac) -> Result<Self, NyxError> {
        let frame = almanac
            .frame_info(cfg.frame)
            .map_err(|e| NyxError::FromAlmanacError {
                source: Box::new(AlmanacError::GenericError { err: e.to_string() }),
                action: "fetching gravity field frame",
            })?;

        if !cfg.gunzipped && cfg.filepath.ends_with(".cof")
            || cfg.gunzipped && cfg.filepath.ends_with(".cof.gz")
        {
            Self::from_cof(cfg.filepath, cfg.degree, cfg.order, cfg.gunzipped, frame)
        } else {
            Self::from_shadr(cfg.filepath, cfg.degree, cfg.order, cfg.gunzipped, frame)
        }
    }

    /// Initialize `GravityFieldData` with a custom J2 value
    pub fn from_j2(j2: f64, frame: Frame) -> GravityFieldData {
        let mut c_nm = DMatrix::from_element(3, 3, 0.0);
        c_nm[(2, 0)] = j2;

        GravityFieldData {
            degree: 2,
            order: 0,
            c_nm,
            s_nm: DMatrix::from_element(3, 3, 0.0),
            frame,
        }
    }

    /// Initialize `GravityFieldData` from the file path (must be a gunzipped file)
    ///
    /// Gravity models provided by `nyx`:
    /// + EMG2008 to 2190 for Earth (tide free)
    /// + Moon to 1500 (from SHADR file)
    /// + Mars to 120 (from SHADR file)
    /// + Venus to 150 (from SHADR file)
    pub fn from_shadr<P: AsRef<Path> + Debug>(
        filepath: P,
        degree: usize,
        order: usize,
        gunzipped: bool,
        frame: Frame,
    ) -> Result<GravityFieldData, NyxError> {
        Self::load(
            filepath, gunzipped, true, //SHADR has a header which we ignore
            degree, order, frame,
        )
    }

    pub fn from_cof<P: AsRef<Path> + Debug>(
        filepath: P,
        degree: usize,
        order: usize,
        gunzipped: bool,
        frame: Frame,
    ) -> Result<GravityFieldData, NyxError> {
        let mut f = File::open(&filepath).map_err(|_| NyxError::FileUnreadable {
            msg: format!("File not found: {filepath:?}"),
        })?;
        let mut buffer = vec![0; 0];
        if gunzipped {
            let mut d = GzDecoder::new(f);
            d.read_to_end(&mut buffer)
                .map_err(|_| NyxError::FileUnreadable {
                    msg: "could not read file as gunzip".to_string(),
                })?;
        } else {
            f.read_to_end(&mut buffer)
                .map_err(|_| NyxError::FileUnreadable {
                    msg: "could not read file to end".to_string(),
                })?;
        }

        let data_as_str = String::from_utf8(buffer).map_err(|_| NyxError::FileUnreadable {
            msg: "could not decode file contents as utf8".to_string(),
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
                            return Err(NyxError::FileUnreadable {
                                msg: format!(
                                    "Harmonics file:
                                could not parse degree `{item}` on line {lno}"
                                ),
                            });
                        }
                    },
                    2 => match usize::from_str(item) {
                        Ok(val) => cur_order = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable {
                                msg: format!(
                                    "Harmonics file:
                                could not parse order `{item}` on line {lno}"
                                ),
                            });
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
                                    return Err(NyxError::FileUnreadable {
                                        msg: format!(
                                            "Harmonics file:
                                        could not parse C_nm `{item}` on line {lno}"
                                        ),
                                    });
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
                                            return Err(NyxError::FileUnreadable {
                                                msg: format!(
                                                    "Harmonics file:
                                                could not parse C_nm `{item}` on line {lno}"
                                                ),
                                            });
                                        }
                                    }
                                    // That mean we have five minus signs, so both the C and S are negative.
                                    let s_nm_str = "-".to_owned() + parts[3] + "-" + parts[4];
                                    match f64::from_str(&s_nm_str) {
                                        Ok(val) => s_nm = val,
                                        Err(_) => {
                                            return Err(NyxError::FileUnreadable {
                                                msg: format!(
                                                    "Harmonics file:
                                                could not parse S_nm `{item}` on line {lno}"
                                                ),
                                            });
                                        }
                                    }
                                } else {
                                    // That mean we have fouve minus signs, and since both values are concatenated, C_nm is positive and S_nm is negative
                                    let c_nm_str = parts[0].to_owned() + "-" + parts[1];
                                    match f64::from_str(&c_nm_str) {
                                        Ok(val) => c_nm = val,
                                        Err(_) => {
                                            return Err(NyxError::FileUnreadable {
                                                msg: format!(
                                                    "Harmonics file:
                                                could not parse C_nm `{item}` on line {lno}"
                                                ),
                                            });
                                        }
                                    }
                                    // That mean we have five minus signs, so both the C and S are negative.
                                    let s_nm_str = "-".to_owned() + parts[2] + "-" + parts[3];
                                    match f64::from_str(&s_nm_str) {
                                        Ok(val) => s_nm = val,
                                        Err(_) => {
                                            return Err(NyxError::FileUnreadable {
                                                msg: format!(
                                                    "Harmonics file:
                                                could not parse S_nm `{item}` on line {lno}"
                                                ),
                                            });
                                        }
                                    }
                                }
                            } else {
                                // We only have the first item, and that's the C_nm
                                match f64::from_str(item) {
                                    Ok(val) => c_nm = val,
                                    Err(_) => {
                                        return Err(NyxError::FileUnreadable {
                                            msg: format!(
                                                "Harmonics file:
                                            could not parse C_nm `{item}` on line {lno}"
                                            ),
                                        });
                                    }
                                }
                            }
                        }
                    }
                    4 => match f64::from_str(item) {
                        // If this exists, then the S_nm is positive.
                        Ok(val) => s_nm = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable {
                                msg: format!(
                                    "Harmonics file:
                                could not parse S_nm `{item}` on line {lno}"
                                ),
                            });
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
                "{filepath:?} only contained (degree, order) of ({max_degree}, {max_order}) instead of requested ({degree}, {order})"
            );
        } else {
            info!("{filepath:?} loaded with (degree, order) = ({degree}, {order})");
        }
        Ok(GravityFieldData {
            degree: max_degree,
            order: max_order,
            c_nm: c_nm_mat,
            s_nm: s_nm_mat,
            frame,
        })
    }

    /// `load` handles the actual loading in memory.
    fn load<P: AsRef<Path> + Debug>(
        filepath: P,
        gunzipped: bool,
        skip_first_line: bool,
        degree: usize,
        order: usize,
        frame: Frame,
    ) -> Result<GravityFieldData, NyxError> {
        let mut f = File::open(&filepath).map_err(|_| NyxError::FileUnreadable {
            msg: format!("File not found: {filepath:?}"),
        })?;
        let mut buffer = vec![0; 0];
        if gunzipped {
            let mut d = GzDecoder::new(f);
            d.read_to_end(&mut buffer)
                .map_err(|_| NyxError::FileUnreadable {
                    msg: "could not read file as gunzip".to_string(),
                })?;
        } else {
            f.read_to_end(&mut buffer)
                .map_err(|_| NyxError::FileUnreadable {
                    msg: "could not read file to end".to_string(),
                })?;
        }

        let data_as_str = String::from_utf8(buffer).map_err(|_| NyxError::FileUnreadable {
            msg: "could not decode file contents as utf8".to_string(),
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
            for (ino, item) in line.replace(',', " ").split_whitespace().enumerate() {
                match ino {
                    0 => match usize::from_str(item) {
                        Ok(val) => cur_degree = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable {
                                msg: format!(
                                    "Harmonics file:
                                could not parse degree on line {lno} (`{item}`)",
                                ),
                            });
                        }
                    },
                    1 => match usize::from_str(item) {
                        Ok(val) => cur_order = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable {
                                msg: format!(
                                    "Harmonics file:
                                could not parse order on line {lno} (`{item}`)"
                                ),
                            });
                        }
                    },
                    2 => match f64::from_str(&item.replace('D', "E")) {
                        Ok(val) => c_nm = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable {
                                msg: format!(
                                    "Harmonics file:
                                could not parse C_nm `{item}` on line {lno}"
                                ),
                            });
                        }
                    },
                    3 => match f64::from_str(&item.replace('D', "E")) {
                        Ok(val) => s_nm = val,
                        Err(_) => {
                            return Err(NyxError::FileUnreadable {
                                msg: format!(
                                    "Harmonics file:
                                could not parse S_nm `{item}` on line {lno}"
                                ),
                            });
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
                "{filepath:?} only contained (degree, order) of ({max_degree}, {max_order}) instead of requested ({degree}, {order})",
            );
        } else {
            info!("{filepath:?} loaded with (degree, order) = ({degree}, {order})");
        }
        Ok(GravityFieldData {
            order: max_order,
            degree: max_degree,
            c_nm: c_nm_mat,
            s_nm: s_nm_mat,
            frame,
        })
    }

    /// Returns the maximum order of this gravity potential storage (Jnm=Jn2,Jn3...)
    pub fn max_order_m(&self) -> usize {
        self.order
    }

    /// Returns the maximum degree of this gravity potential storage (Jn=J2,J3...)
    pub fn max_degree_n(&self) -> usize {
        self.degree
    }

    /// Returns the C_nm and S_nm for the provided order and degree.
    pub fn cs_nm(&self, degree: usize, order: usize) -> (f64, f64) {
        (self.c_nm[(degree, order)], self.s_nm[(degree, order)])
    }
}

impl StaticType for GravityFieldConfig {
    fn static_type() -> SimpleType {
        let mut fields = HashMap::new();

        fields.insert("filepath".to_string(), String::static_type());
        fields.insert("gunzipped".to_string(), bool::static_type());
        fields.insert("degree".to_string(), usize::static_type());
        fields.insert("order".to_string(), usize::static_type());

        SimpleType::Record(fields)
    }
}

#[cfg(test)]
#[test]
fn test_load_harmonic_files() {
    use anise::constants::frames::IAU_EARTH_FRAME;

    let data_folder: PathBuf = [env!("CARGO_MANIFEST_DIR"), "../data/01_planetary"]
        .iter()
        .collect();

    GravityFieldData::from_cof(
        data_folder.join("JGM3.cof.gz"),
        50,
        50,
        true,
        IAU_EARTH_FRAME,
    )
    .expect("could not load JGM3");

    GravityFieldData::from_shadr(
        data_folder.join("EGM2008_to2190_TideFree.gz"),
        120,
        120,
        true,
        IAU_EARTH_FRAME,
    )
    .expect("could not load EGM2008");

    GravityFieldData::from_shadr(
        data_folder.join("Luna_jggrx_1500e_sha.tab.gz"),
        1500,
        1500,
        true,
        IAU_EARTH_FRAME,
    )
    .expect("could not load jggrx");
}
