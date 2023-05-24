/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::TrajError;
use super::{ExportCfg, Traj};
use crate::cosmic::{Cosm, Frame, Orbit};
use crate::errors::NyxError;
use crate::md::prelude::StateParameter;
use crate::md::EventEvaluator;
use crate::time::Epoch;
use crate::time::TimeUnits;
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufReader, Read};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;
use std::time::Instant;

impl Traj<Orbit> {
    /// Allows converting the source trajectory into the (almost) equivalent trajectory in another frame.
    /// This simply converts each state into the other frame and may lead to aliasing due to the Nyquist–Shannon sampling theorem.
    #[allow(clippy::map_clone)]
    pub fn to_frame(&self, new_frame: Frame, cosm: Arc<Cosm>) -> Result<Self, NyxError> {
        if self.states.is_empty() {
            return Err(NyxError::Trajectory(TrajError::CreationError(
                "No trajectory to convert".to_string(),
            )));
        }
        let start_instant = Instant::now();
        let mut traj = Self::new();
        for state in &self.states {
            traj.states.push(cosm.frame_chg(state, new_frame));
        }
        traj.finalize();

        info!(
            "Converted trajectory from {} to {} in {} ms: {traj}",
            self.first().frame,
            new_frame,
            (Instant::now() - start_instant).as_millis()
        );
        Ok(traj)
    }

    /// Exports this trajectory to the provided filename in parquet format with only the epoch, the geodetic latitude, longitude, and height at one state per minute.
    /// Must provide a body fixed frame to correctly compute the latitude and longitude.
    #[allow(clippy::identity_op)]
    pub fn to_groundtrack_parquet<P: AsRef<Path>>(
        &self,
        path: P,
        body_fixed_frame: Frame,
        events: Option<Vec<&dyn EventEvaluator<Orbit>>>,
        metadata: Option<HashMap<String, String>>,
        cosm: Arc<Cosm>,
    ) -> Result<PathBuf, Box<dyn Error>> {
        let traj = self.to_frame(body_fixed_frame, cosm)?;

        let mut cfg = ExportCfg::default();
        cfg.append_field(StateParameter::GeodeticLatitude);
        cfg.append_field(StateParameter::GeodeticLongitude);
        cfg.append_field(StateParameter::GeodeticHeight);
        cfg.append_field(StateParameter::Rmag);
        cfg.set_step(1.minutes());
        cfg.metadata = metadata;

        traj.to_parquet(path, events, cfg)
    }

    /// Initialize a new orbit trajectory from the path to a CCSDS OEM file.
    ///
    /// # Limitations
    /// 1. Only text versions of the OEM format are supported
    /// 2. The covariance information, if present, is ignored
    /// 3. Only one spacecraft per OEM file is supported.
    ///
    /// # Thanks
    /// GPT-4 because I didn't want to spend too much time coding this up since it'll be a feature in ANISE.
    pub fn from_oem_file<P: AsRef<Path>>(path: P) -> Result<Self, NyxError> {
        let cosm = Cosm::de438();
        // Open the file
        let file =
            File::open(path).map_err(|e| NyxError::CCSDS(format!("File opening error: {e}")))?;
        let mut reader = BufReader::new(file);

        // Read the file contents into a buffer
        let mut buffer = String::new();
        reader
            .read_to_string(&mut buffer)
            .map_err(|e| NyxError::CCSDS(format!("File read error: {e}")))?;

        // Split the file contents into lines
        let lines: Vec<&str> = buffer.lines().collect();

        // Parse the Orbit Element messages
        let mut frame: Option<Frame> = None;
        let mut time_system = String::new();

        let ignored_tokens = ["CCSDS_OMM_VERS", "CREATION_DATE", "ORIGINATOR"];

        let mut traj = Self::default();

        let mut start = false;

        'lines: for (lno, line) in lines.iter().enumerate() {
            for tok in &ignored_tokens {
                if line.starts_with(tok) {
                    // Ignore this token
                    continue 'lines;
                }
            }
            if line.starts_with("OBJECT_NAME") {
                // Extract the object ID from the line
                let parts: Vec<&str> = line.split("=").collect();
                let name = parts[1].trim().to_string();
                debug!("[line: {lno}] Found object {name}");
                traj.name = Some(name);
            } else if line.starts_with("REF_FRAME") {
                let parts: Vec<&str> = line.split("=").collect();
                let mut ref_frame = parts[1].trim();
                if ref_frame == "ICRF" {
                    ref_frame = "EME2000";
                }
                frame = Some(cosm.try_frame(ref_frame)?);
            } else if line.starts_with("TIME_SYSTEM") {
                let parts: Vec<&str> = line.split("=").collect();
                time_system = parts[1].trim().to_string();
                debug!("[line: {lno}] Found time system `{time_system}`");
            } else if line.starts_with("META_STOP") {
                // We can start parsing now
                start = true;
            } else if !line.is_empty() && start {
                // Split the line into components
                let parts: Vec<&str> = line.split_whitespace().collect();

                // Extract the values
                let epoch_str = format!("{} {time_system}", parts[0].to_string());
                match parts[1].parse::<f64>() {
                    Ok(x_km) => {
                        // Look good!
                        let y_km = parts[2].parse::<f64>().unwrap();
                        let z_km = parts[3].parse::<f64>().unwrap();
                        let vx_km_s = parts[4].parse::<f64>().unwrap();
                        let vy_km_s = parts[5].parse::<f64>().unwrap();
                        let vz_km_s = parts[6].parse::<f64>().unwrap();

                        debug!("[line: {lno}] Parsing epoch `{epoch_str}`");
                        let orbit = Orbit {
                            epoch: Epoch::from_str(epoch_str.trim()).map_err(|e| {
                                NyxError::CCSDS(format!("Parsing epoch error: {e}"))
                            })?,
                            x_km,
                            y_km,
                            z_km,
                            vx_km_s,
                            vy_km_s,
                            vz_km_s,
                            frame: frame.unwrap(),
                            stm: None,
                        };

                        traj.states.push(orbit);
                    }
                    Err(_) => {
                        // Probably a comment
                        continue;
                    }
                };
            }
        }

        Ok(traj)
    }
}

#[test]
fn test_load_oem_leo() {
    use pretty_env_logger;
    use std::env;

    // All three samples were taken from https://github.com/bradsease/oem/blob/main/tests/samples/real/
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "tests",
        "ccsds",
        "oem",
        "LEO_10s.oem",
    ]
    .iter()
    .collect();

    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let leo_traj = Traj::from_oem_file(path).unwrap();

    assert_eq!(leo_traj.states.len(), 361);
    assert_eq!(leo_traj.name.unwrap(), "TEST_OBJ".to_string());
}

#[test]
fn test_load_oem_meo() {
    use pretty_env_logger;
    use std::env;

    // All three samples were taken from https://github.com/bradsease/oem/blob/main/tests/samples/real/
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "tests",
        "ccsds",
        "oem",
        "MEO_60s.oem",
    ]
    .iter()
    .collect();

    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let leo_traj = Traj::from_oem_file(path).unwrap();

    assert_eq!(leo_traj.states.len(), 61);
    assert_eq!(leo_traj.name.unwrap(), "TEST_OBJ".to_string());
}

#[test]
fn test_load_oem_geo() {
    use pretty_env_logger;
    use std::env;

    // All three samples were taken from https://github.com/bradsease/oem/blob/main/tests/samples/real/
    let path: PathBuf = [
        env!("CARGO_MANIFEST_DIR"),
        "data",
        "tests",
        "ccsds",
        "oem",
        "GEO_20s.oem",
    ]
    .iter()
    .collect();

    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    let leo_traj = Traj::from_oem_file(path).unwrap();

    assert_eq!(leo_traj.states.len(), 181);
    assert_eq!(leo_traj.name.unwrap(), "TEST_OBJ".to_string());
}
