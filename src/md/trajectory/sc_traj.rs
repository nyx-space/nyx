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

use anise::prelude::{Almanac, Frame, Orbit};
use hifitime::Epoch;
use snafu::ResultExt;

use super::TrajError;
use super::{ExportCfg, Traj};
use crate::cosmic::Spacecraft;
use crate::errors::{FromAlmanacSnafu, NyxError};
use crate::md::prelude::StateParameter;
use crate::md::EventEvaluator;
use crate::time::{Duration, TimeUnits};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;
#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

impl Traj<Spacecraft> {
    /// Allows converting the source trajectory into the (almost) equivalent trajectory in another frame
    #[allow(clippy::map_clone)]
    pub fn to_frame(&self, new_frame: Frame, almanac: Arc<Almanac>) -> Result<Self, NyxError> {
        if self.states.is_empty() {
            return Err(NyxError::Trajectory {
                source: TrajError::CreationError {
                    msg: "No trajectory to convert".to_string(),
                },
            });
        }

        #[cfg(not(target_arch = "wasm32"))]
        let start_instant = Instant::now();
        let mut traj = Self::new();
        for state in &self.states {
            let new_orbit =
                almanac
                    .transform_to(state.orbit, new_frame, None)
                    .context(FromAlmanacSnafu {
                        action: "transforming trajectory into new frame",
                    })?;
            traj.states.push(state.with_orbit(new_orbit));
        }
        traj.finalize();

        #[cfg(not(target_arch = "wasm32"))]
        info!(
            "Converted trajectory from {} to {} in {} ms: {traj}",
            self.first().orbit.frame,
            new_frame,
            (Instant::now() - start_instant).as_millis()
        );

        #[cfg(target_arch = "wasm32")]
        info!(
            "Converted trajectory from {} to {}: {traj}",
            self.first().orbit.frame,
            new_frame,
        );

        Ok(traj)
    }

    /// A shortcut to `to_parquet_with_cfg`
    pub fn to_parquet_with_step<P: AsRef<Path>>(
        &self,
        path: P,
        step: Duration,
        almanac: Arc<Almanac>,
    ) -> Result<(), Box<dyn Error>> {
        self.to_parquet_with_cfg(
            path,
            ExportCfg {
                step: Some(step),
                ..Default::default()
            },
            almanac,
        )?;

        Ok(())
    }

    /// Exports this trajectory to the provided filename in parquet format with only the epoch, the geodetic latitude, longitude, and height at one state per minute.
    /// Must provide a body fixed frame to correctly compute the latitude and longitude.
    #[allow(clippy::identity_op)]
    pub fn to_groundtrack_parquet<P: AsRef<Path>>(
        &self,
        path: P,
        body_fixed_frame: Frame,
        events: Option<Vec<&dyn EventEvaluator<Spacecraft>>>,
        metadata: Option<HashMap<String, String>>,
        almanac: Arc<Almanac>,
    ) -> Result<PathBuf, Box<dyn Error>> {
        let traj = self.to_frame(body_fixed_frame, almanac.clone())?;

        let mut cfg = ExportCfg::builder()
            .step(1.minutes())
            .fields(vec![
                StateParameter::GeodeticLatitude,
                StateParameter::GeodeticLongitude,
                StateParameter::GeodeticHeight,
                StateParameter::Rmag,
            ])
            .build();
        cfg.metadata = metadata;

        traj.to_parquet(path, events, cfg, almanac)
    }

    /// Convert this spacecraft trajectory into an Orbit trajectory, loosing all references to the spacecraft
    pub fn downcast(&self) -> Self {
        unimplemented!()
    }

    /// Initialize a new spacecraft trajectory from the path to a CCSDS OEM file.
    ///
    /// CCSDS OEM only contains the orbit information, so you must provide a template spacecraft since we'll upcast the orbit trajectory into a spacecraft trajectory.
    pub fn from_oem_file<P: AsRef<Path>>(path: P, template: Spacecraft) -> Result<Self, NyxError> {
        // Open the file
        let file = File::open(path).map_err(|e| NyxError::CCSDS {
            msg: format!("File opening error: {e}"),
        })?;
        let reader = BufReader::new(file);

        // Parse the Orbit Element messages
        let mut time_system = String::new();

        let ignored_tokens: HashSet<_> = [
            "CCSDS_OMM_VERS".to_string(),
            "CREATION_DATE".to_string(),
            "ORIGINATOR".to_string(),
        ]
        .into();

        let mut traj = Self::default();

        let mut parse = false;

        let mut center_name = None;
        let mut orient_name = None;

        'lines: for (lno, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| NyxError::CCSDS {
                msg: format!("File read error: {e}"),
            })?;
            let line = line.trim();
            if line.is_empty() {
                continue;
            }

            if ignored_tokens.iter().any(|t| line.starts_with(t)) {
                continue 'lines;
            }
            if line.starts_with("OBJECT_NAME") {
                // Extract the object ID from the line
                let parts: Vec<&str> = line.split('=').collect();
                let name = parts[1].trim().to_string();
                debug!("[line: {}] Found object {name}", lno + 1);
                traj.name = Some(name);
            } else if line.starts_with("CENTER_NAME") {
                let parts: Vec<&str> = line.split('=').collect();
                center_name = Some(parts[1].trim().to_owned());
            } else if line.starts_with("REF_FRAME") {
                let parts: Vec<&str> = line.split('=').collect();
                orient_name = Some(parts[1].trim().to_owned());
            } else if line.starts_with("TIME_SYSTEM") {
                let parts: Vec<&str> = line.split('=').collect();
                time_system = parts[1].trim().to_string();
                debug!("[line: {}] Found time system `{time_system}`", lno + 1);
            } else if line.starts_with("META_STOP") {
                // We can start parsing now
                parse = true;
            } else if line.starts_with("META_START") {
                // Stop the parsing
                parse = false;
            } else if line.starts_with("COVARIANCE_START") {
                // Stop the parsing
                warn!("[line: {}] Skipping covariance in OEM parsing", lno + 1);
                parse = false;
            } else if parse {
                let frame = Frame::from_name(
                    center_name.clone().unwrap().as_str(),
                    orient_name.clone().unwrap().as_str(),
                )
                .map_err(|e| NyxError::CCSDS {
                    msg: format!("frame error `{center_name:?} {orient_name:?}`: {e}"),
                })?;
                // Split the line into components
                let parts: Vec<&str> = line.split_whitespace().collect();

                if parts.len() < 7 {
                    debug!("[line: {}] Could not understand `{parts:?}`", lno + 1);
                } else {
                    // Extract the values
                    let epoch_str = format!("{} {time_system}", parts[0]);
                    match parts[1].parse::<f64>() {
                        Ok(x_km) => {
                            // Look good!
                            let y_km = parts[2].parse::<f64>().unwrap();
                            let z_km = parts[3].parse::<f64>().unwrap();
                            let vx_km_s = parts[4].parse::<f64>().unwrap();
                            let vy_km_s = parts[5].parse::<f64>().unwrap();
                            let vz_km_s = parts[6].parse::<f64>().unwrap();

                            let epoch =
                                Epoch::from_str(epoch_str.trim()).map_err(|e| NyxError::CCSDS {
                                    msg: format!("Parsing epoch error: {e}"),
                                })?;

                            let orbit = Orbit::new(
                                x_km, y_km, z_km, vx_km_s, vy_km_s, vz_km_s, epoch, frame,
                            );

                            traj.states.push(template.with_orbit(orbit));
                        }
                        Err(_) => {
                            // Probably a comment
                            debug!("[line: {}] Could not parse `{parts:?}`", lno + 1);
                            continue;
                        }
                    };
                }
            }
        }

        traj.finalize();

        Ok(traj)
    }
}
