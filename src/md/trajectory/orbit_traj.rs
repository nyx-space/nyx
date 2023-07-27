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
use crate::io::watermark::prj_name_ver;
use crate::md::prelude::StateParameter;
use crate::md::EventEvaluator;
use crate::time::{Epoch, Format, Formatter, TimeUnits};
use crate::{Spacecraft, State};
use std::collections::{HashMap, HashSet};
use std::error::Error;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::str::FromStr;
use std::sync::Arc;
#[cfg(not(target_arch = "wasm32"))]
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

        #[cfg(not(target_arch = "wasm32"))]
        let start_instant = Instant::now();
        let mut traj = Self::new();
        for state in &self.states {
            traj.states.push(cosm.frame_chg(state, new_frame));
        }
        traj.finalize();

        #[cfg(not(target_arch = "wasm32"))]
        info!(
            "Converted trajectory from {} to {} in {} ms: {traj}",
            self.first().frame,
            new_frame,
            (Instant::now() - start_instant).as_millis()
        );

        #[cfg(target_arch = "wasm32")]
        info!(
            "Converted trajectory from {} to {}: {traj}",
            self.first().frame,
            new_frame,
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

        traj.to_parquet(path, events, cfg)
    }

    /// Convert this orbit trajectory into a spacecraft trajectory by copying the provided template and setting its orbit state to that of each state of the trajectory
    pub fn upcast(&self, template: Spacecraft) -> Traj<Spacecraft> {
        let mut out = Traj::new();
        for orbit in &self.states {
            out.states.push(template.with_orbit(*orbit));
        }
        out
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
        let reader = BufReader::new(file);

        // Parse the Orbit Element messages
        let mut frame: Option<Frame> = None;
        let mut time_system = String::new();

        let ignored_tokens: HashSet<_> = [
            "CCSDS_OMM_VERS".to_string(),
            "CREATION_DATE".to_string(),
            "ORIGINATOR".to_string(),
        ]
        .into();

        let mut traj = Self::default();

        let mut parse = false;

        'lines: for (lno, line) in reader.lines().enumerate() {
            let line = line.map_err(|e| NyxError::CCSDS(format!("File read error: {e}")))?;
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
                let center_name = parts[1].trim();
                frame = Some(cosm.try_frame(&format!("{center_name} J2000"))?);
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

    pub fn to_oem_file<P: AsRef<Path>>(
        &self,
        path: P,
        cfg: ExportCfg,
    ) -> Result<PathBuf, NyxError> {
        if self.states.is_empty() {
            return Err(NyxError::CCSDS(
                "Cannot export an empty trajectory to OEM".to_string(),
            ));
        }
        let tick = Epoch::now().unwrap();
        info!("Exporting trajectory to CCSDS OEM file...");

        // Grab the path here before we move stuff.
        let path_buf = cfg.actual_path(path);

        let metadata = cfg.metadata.unwrap_or(HashMap::new());

        let file = File::create(&path_buf)
            .map_err(|e| NyxError::CCSDS(format!("File creation error: {e}")))?;
        let mut writer = BufWriter::new(file);

        let err_hdlr = |e| NyxError::CCSDS(format!("Could not write: {e}"));

        // Build the states iterator -- this does require copying the current states but I can't either get a reference or a copy of all the states.
        let states = if cfg.start_epoch.is_some() || cfg.end_epoch.is_some() || cfg.step.is_some() {
            // Must interpolate the data!
            let start = cfg.start_epoch.unwrap_or_else(|| self.first().epoch());
            let end = cfg.end_epoch.unwrap_or_else(|| self.last().epoch());
            let step = cfg.step.unwrap_or_else(|| 1.minutes());
            self.every_between(step, start, end).collect()
        } else {
            self.states.to_vec()
        };

        // Epoch formmatter.
        let iso8601_no_ts = Format::from_str("%Y-%m-%dT%H:%M:%S.%f").unwrap();

        // Write mandatory metadata
        writeln!(writer, "CCSDS_OMM_VERS = 2.0").map_err(err_hdlr)?;
        writeln!(
            writer,
            "CREATION_DATE = {}",
            Formatter::new(Epoch::now().unwrap(), iso8601_no_ts)
        )
        .map_err(err_hdlr)?;
        writeln!(
            writer,
            "ORIGINATOR = {}\n",
            metadata
                .get("originator")
                .unwrap_or(&"Nyx Space".to_string())
        )
        .map_err(err_hdlr)?;

        writeln!(writer, "META_START").map_err(err_hdlr)?;
        // Write optional metadata
        if let Some(object_name) = metadata.get("object_name") {
            writeln!(writer, "OBJECT_NAME = {}", object_name).map_err(err_hdlr)?;
        } else if let Some(object_name) = &self.name {
            writeln!(writer, "OBJECT_NAME = {}", object_name).map_err(err_hdlr)?;
        }

        let frame_str = states[0].frame.to_string();
        let splt: Vec<&str> = frame_str.split(' ').collect();
        let center = splt[0];
        let ref_frame = frame_str.replace(center, " ");
        writeln!(
            writer,
            "REF_FRAME = {}",
            match ref_frame.trim() {
                "J2000" => "ICRF",
                _ => ref_frame.trim(),
            }
        )
        .map_err(err_hdlr)?;

        writeln!(writer, "CENTER_NAME = {center}",).map_err(err_hdlr)?;

        writeln!(writer, "TIME_SYSTEM = {}", states[0].epoch.time_scale).map_err(err_hdlr)?;

        writeln!(
            writer,
            "START_TIME = {}",
            Formatter::new(states[0].epoch(), iso8601_no_ts)
        )
        .map_err(err_hdlr)?;
        writeln!(
            writer,
            "USEABLE_START_TIME = {}",
            Formatter::new(states[0].epoch(), iso8601_no_ts)
        )
        .map_err(err_hdlr)?;
        writeln!(
            writer,
            "USEABLE_STOP_TIME = {}",
            Formatter::new(states[states.len() - 1].epoch(), iso8601_no_ts)
        )
        .map_err(err_hdlr)?;
        writeln!(
            writer,
            "STOP_TIME = {}",
            Formatter::new(states[states.len() - 1].epoch(), iso8601_no_ts)
        )
        .map_err(err_hdlr)?;

        writeln!(writer, "META_STOP\n").map_err(err_hdlr)?;

        writeln!(
            writer,
            "COMMENT Generated by {} provided in AGPLv3 license -- https://nyxspace.com/\n",
            prj_name_ver()
        )
        .map_err(err_hdlr)?;

        for state in &states {
            writeln!(
                writer,
                "{} {:E} {:E} {:E} {:E} {:E} {:E}",
                Formatter::new(state.epoch, iso8601_no_ts),
                state.x_km,
                state.y_km,
                state.z_km,
                state.vx_km_s,
                state.vy_km_s,
                state.vz_km_s
            )
            .map_err(err_hdlr)?;
        }

        #[allow(clippy::writeln_empty_string)]
        writeln!(writer, "").map_err(err_hdlr)?;

        // Return the path this was written to
        let tock_time = Epoch::now().unwrap() - tick;
        info!(
            "Trajectory written to {} in {tock_time}",
            path_buf.display()
        );
        Ok(path_buf)
    }
}

#[cfg(test)]
mod ut_ccsds_oem {

    use crate::md::prelude::{Cosm, OrbitalDynamics, Propagator};
    use crate::time::{Epoch, TimeUnits};
    use crate::{io::ExportCfg, md::prelude::Traj, Orbit};
    use pretty_env_logger;
    use std::env;
    use std::str::FromStr;
    use std::{collections::HashMap, path::PathBuf};

    #[test]
    fn test_load_oem_leo() {
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

        let _ = pretty_env_logger::try_init();

        let traj: Traj<Orbit> = Traj::<Orbit>::from_oem_file(path).unwrap();

        // This trajectory has two duplicate epochs, which should be removed by the call to finalize()
        assert_eq!(traj.states.len(), 361);
        assert_eq!(traj.name.unwrap(), "TEST_OBJ".to_string());
    }

    #[test]
    fn test_load_oem_meo() {
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

        let _ = pretty_env_logger::try_init();

        let traj = Traj::<Orbit>::from_oem_file(path).unwrap();

        assert_eq!(traj.states.len(), 61);
        assert_eq!(traj.name.unwrap(), "TEST_OBJ".to_string());
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

        let _ = pretty_env_logger::try_init();

        let traj: Traj<Orbit> = Traj::<Orbit>::from_oem_file(path).unwrap();

        assert_eq!(traj.states.len(), 181);
        assert_eq!(traj.name.as_ref().unwrap(), &"TEST_OBJ".to_string());

        // Reexport this to CCSDS.
        let cfg = ExportCfg::builder()
            .timestamp(true)
            .metadata(HashMap::from([
                ("originator".to_string(), "Test suite".to_string()),
                ("object_name".to_string(), "TEST_OBJ".to_string()),
            ]))
            .build();

        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "output_data",
            "GEO_20s_rebuilt.oem",
        ]
        .iter()
        .collect();

        let out_path = traj.to_oem_file(path.clone(), cfg).unwrap();
        // And reload, make sure we have the same data.
        let traj_reloaded: Traj<Orbit> = Traj::<Orbit>::from_oem_file(out_path).unwrap();

        assert_eq!(traj_reloaded, traj);

        // Now export after trimming one state on either end
        let cfg = ExportCfg::builder()
            .timestamp(true)
            .metadata(HashMap::from([
                ("originator".to_string(), "Test suite".to_string()),
                ("object_name".to_string(), "TEST_OBJ".to_string()),
            ]))
            .step(20.seconds())
            .start_epoch(traj.first().epoch + 1.seconds())
            .end_epoch(traj.last().epoch - 1.seconds())
            .build();
        let out_path = traj.to_oem_file(path, cfg).unwrap();
        // And reload, make sure we have the same data.
        let traj_reloaded: Traj<Orbit> = Traj::<Orbit>::from_oem_file(out_path).unwrap();

        // Note that the number of states has changed because we interpolated with a step similar to the original one but
        // we started with a different time.
        assert_eq!(traj_reloaded.states.len(), traj.states.len() - 1);
        assert_eq!(
            traj_reloaded.first().epoch,
            traj.first().epoch + 1.seconds()
        );
        // Note: because we used a fixed step, the last epoch is actually an offset of step size - end offset
        // from the original trajectory
        assert_eq!(traj_reloaded.last().epoch, traj.last().epoch - 19.seconds());
    }

    #[test]
    fn test_moon_frame_long_prop() {
        let cosm = Cosm::de438();
        let epoch = Epoch::from_str("2022-06-13T12:00:00").unwrap();
        let orbit = Orbit::keplerian_altitude(
            350.0,
            0.02,
            30.0,
            45.0,
            85.0,
            0.0,
            epoch,
            cosm.frame("Moon J2000"),
        );

        let mut traj = Propagator::default_dp78(OrbitalDynamics::two_body())
            .with(orbit)
            .for_duration_with_traj(45.days())
            .unwrap()
            .1;
        // Set the name of this object
        traj.name = Some("TEST_MOON_OBJ".to_string());

        // Export CCSDS OEM file
        let path: PathBuf = [env!("CARGO_MANIFEST_DIR"), "output_data", "moon_45days.oem"]
            .iter()
            .collect();

        let out_path = traj.to_oem_file(path, ExportCfg::default()).unwrap();

        // And reload
        let traj_reloaded: Traj<Orbit> = Traj::<Orbit>::from_oem_file(out_path).unwrap();

        assert_eq!(traj, traj_reloaded);
    }
}
