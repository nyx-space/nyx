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

use super::TrajError;
use super::{ExportCfg, Traj};
use crate::cosmic::Spacecraft;
use crate::errors::{FromAlmanacSnafu, NyxError};
use crate::io::{InputOutputError, MissingDataSnafu, ParquetSnafu, StdIOSnafu};
use crate::md::prelude::{Interpolatable, StateParameter};
use crate::time::{Duration, Epoch, TimeUnits};
use crate::State;
use anise::analysis::prelude::OrbitalElement;
use anise::astro::Aberration;
use anise::ephemerides::ephemeris::Ephemeris;
use anise::ephemerides::EphemerisError;
use anise::errors::AlmanacError;
use anise::prelude::{Almanac, Frame};
use arrow::array::RecordBatchReader;
use arrow::array::{Float64Array, StringArray};
use hifitime::TimeSeries;
use log::info;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use snafu::{ensure, ResultExt};
use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::path::{Path, PathBuf};
use std::sync::Arc;
#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant;

impl Traj<Spacecraft> {
    /// Builds a new trajectory built from the SPICE BSP (SPK) file loaded in the provided Almanac, provided the start and stop epochs.
    ///
    /// If the start and stop epochs are not provided, then the full domain of the trajectory will be used.
    pub fn from_bsp(
        target_frame: Frame,
        observer_frame: Frame,
        almanac: Arc<Almanac>,
        sc_template: Spacecraft,
        step: Duration,
        start_epoch: Option<Epoch>,
        end_epoch: Option<Epoch>,
        ab_corr: Option<Aberration>,
        name: Option<String>,
    ) -> Result<Self, AlmanacError> {
        let (domain_start, domain_end) =
            almanac
                .spk_domain(target_frame.ephemeris_id)
                .map_err(|e| AlmanacError::Ephemeris {
                    action: "could not fetch domain",
                    source: Box::new(e),
                })?;

        let start_epoch = start_epoch.unwrap_or(domain_start);
        let end_epoch = end_epoch.unwrap_or(domain_end);

        let time_series = TimeSeries::inclusive(start_epoch, end_epoch, step);
        let mut states = Vec::with_capacity(time_series.len());
        for epoch in time_series {
            let orbit = almanac.transform(target_frame, observer_frame, epoch, ab_corr)?;

            states.push(sc_template.with_orbit(orbit));
        }

        Ok(Self { name, states })
    }
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

    /// Exports this trajectory to the provided filename in parquet format with only the epoch, the geodetic latitude, longitude, and height at one state per minute.
    /// Must provide a body fixed frame to correctly compute the latitude and longitude.
    #[allow(clippy::identity_op)]
    pub fn to_groundtrack_parquet<P: AsRef<Path>>(
        &self,
        path: P,
        body_fixed_frame: Frame,
        metadata: Option<HashMap<String, String>>,
        almanac: Arc<Almanac>,
    ) -> Result<PathBuf, Box<dyn Error>> {
        let traj = self.to_frame(body_fixed_frame, almanac)?;

        let mut cfg = ExportCfg::builder()
            .step(1.minutes())
            .fields(vec![
                StateParameter::Element(OrbitalElement::Latitude),
                StateParameter::Element(OrbitalElement::Longitude),
                StateParameter::Element(OrbitalElement::Height),
                StateParameter::Element(OrbitalElement::Rmag),
            ])
            .build();
        cfg.metadata = metadata;

        traj.to_parquet(path, cfg)
    }

    /// Export this spacecraft trajectory estimate to an ANISE Ephemeris
    pub fn to_ephemeris(&self, object_id: String, cfg: ExportCfg) -> Ephemeris {
        let mut ephem = Ephemeris::new(object_id);

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

        for sc_state in &states {
            ephem.insert_orbit(sc_state.orbit());
        }

        ephem
    }

    /// Initialize a new spacecraft trajectory from the path to a CCSDS OEM file.
    ///
    /// CCSDS OEM only contains the orbit information but Nyx builds spacecraft trajectories.
    /// If not spacecraft template is provided, then a default massless spacecraft will be built.
    pub fn from_oem_file<P: AsRef<Path>>(
        path: P,
        tpl_option: Option<Spacecraft>,
    ) -> Result<Self, EphemerisError> {
        // Read the ephemeris
        let ephem = Ephemeris::from_ccsds_oem_file(path)?;
        // Rebuild a trajectory by applying the template
        let template = tpl_option.unwrap_or_default();
        let mut traj = Self::default();
        for record in &ephem {
            traj.states.push(template.with_orbit(record.orbit));
        }
        traj.name = Some(ephem.object_id);

        Ok(traj)
    }

    pub fn to_oem_file<P: AsRef<Path>>(
        &self,
        path: P,
        object_id: String,
        originator: Option<String>,
        object_name: Option<String>,
        cfg: ExportCfg,
    ) -> Result<(), EphemerisError> {
        let ephem = self.to_ephemeris(object_id, cfg);
        ephem.write_ccsds_oem(path, originator, object_name)
    }

    pub fn from_parquet<P: AsRef<Path>>(path: P) -> Result<Self, InputOutputError> {
        let file = File::open(&path).context(StdIOSnafu {
            action: "opening trajectory file",
        })?;

        let builder = ParquetRecordBatchReaderBuilder::try_new(file).unwrap();

        let mut metadata = HashMap::new();
        // Build the custom metadata
        if let Some(file_metadata) = builder.metadata().file_metadata().key_value_metadata() {
            for key_value in file_metadata {
                if !key_value.key.starts_with("ARROW:") {
                    metadata.insert(
                        key_value.key.clone(),
                        key_value.value.clone().unwrap_or("[unset]".to_string()),
                    );
                }
            }
        }

        // Check the schema
        let mut has_epoch = false; // Required
        let mut frame = None;

        let mut found_fields = vec![
            (StateParameter::Element(OrbitalElement::X), false),
            (StateParameter::Element(OrbitalElement::Y), false),
            (StateParameter::Element(OrbitalElement::Z), false),
            (StateParameter::Element(OrbitalElement::VX), false),
            (StateParameter::Element(OrbitalElement::VY), false),
            (StateParameter::Element(OrbitalElement::VZ), false),
            (StateParameter::PropMass(), false),
        ];

        let reader = builder.build().context(ParquetSnafu {
            action: "building output trajectory file",
        })?;

        for field in &reader.schema().fields {
            if field.name().as_str() == "Epoch (UTC)" {
                has_epoch = true;
            } else {
                for potential_field in &mut found_fields {
                    if field.name() == potential_field.0.to_field(None).name() {
                        potential_field.1 = true;
                        if potential_field.0 != StateParameter::PropMass() {
                            if let Some(frame_info) = field.metadata().get("Frame") {
                                // Frame is expected to be serialized as Dhall.
                                match serde_dhall::from_str(frame_info).parse::<Frame>() {
                                    Err(e) => {
                                        return Err(InputOutputError::ParseDhall {
                                            data: frame_info.to_string(),
                                            err: format!("{e}"),
                                        })
                                    }
                                    Ok(deser_frame) => frame = Some(deser_frame),
                                };
                            }
                        }
                        break;
                    }
                }
            }
        }

        ensure!(
            has_epoch,
            MissingDataSnafu {
                which: "Epoch (UTC)"
            }
        );

        ensure!(
            frame.is_some(),
            MissingDataSnafu {
                which: "Frame in metadata"
            }
        );

        for (field, exists) in found_fields.iter().take(found_fields.len() - 1) {
            ensure!(
                exists,
                MissingDataSnafu {
                    which: format!("Missing `{}` field", field.to_field(None).name())
                }
            );
        }

        let sc_compat = found_fields.last().unwrap().1;

        let expected_type = std::any::type_name::<Spacecraft>()
            .split("::")
            .last()
            .unwrap();

        if expected_type == "Spacecraft" {
            ensure!(
                sc_compat,
                MissingDataSnafu {
                    which: format!(
                        "Missing `{}` field",
                        found_fields.last().unwrap().0.to_field(None).name()
                    )
                }
            );
        } else if sc_compat {
            // Not a spacecraft, remove the prop mass
            if let Some(last_field) = found_fields.last_mut() {
                if last_field.0 == StateParameter::PropMass() && last_field.1 {
                    last_field.1 = false;
                }
            }
        }

        // At this stage, we know that the measurement is valid and the conversion is supported.
        let mut traj = Traj::default();

        // Now convert each batch on the fly
        for maybe_batch in reader {
            let batch = maybe_batch.unwrap();

            let epochs = batch
                .column_by_name("Epoch (UTC)")
                .unwrap()
                .as_any()
                .downcast_ref::<StringArray>()
                .unwrap();

            let mut shared_data = vec![];

            for (field, _) in found_fields.iter().take(found_fields.len() - 1) {
                shared_data.push(
                    batch
                        .column_by_name(field.to_field(None).name())
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap(),
                );
            }

            if expected_type == "Spacecraft" {
                // Read the prop only if this is a spacecraft we're building
                shared_data.push(
                    batch
                        .column_by_name("prop_mass (kg)")
                        .unwrap()
                        .as_any()
                        .downcast_ref::<Float64Array>()
                        .unwrap(),
                );
            }

            // Grab the frame -- it should have been serialized with all of the data so we don't need to reload it.

            // Build the states
            for i in 0..batch.num_rows() {
                let mut state = Spacecraft::zeros();
                state.set_epoch(Epoch::from_gregorian_str(epochs.value(i)).map_err(|e| {
                    InputOutputError::Inconsistency {
                        msg: format!("{e} when parsing epoch"),
                    }
                })?);
                state.set_frame(frame.unwrap()); // We checked it was set above with an ensure! call
                state.unset_stm(); // We don't have any STM data, so let's unset this.

                for (j, (param, exists)) in found_fields.iter().enumerate() {
                    if *exists {
                        state.set_value(*param, shared_data[j].value(i)).unwrap();
                    }
                }

                traj.states.push(state);
            }
        }

        // Remove any duplicates that may exist in the imported trajectory.
        traj.finalize();

        Ok(traj)
    }
}

#[cfg(test)]
mod ut_ccsds_oem {

    use crate::md::prelude::{OrbitalDynamics, Propagator, SpacecraftDynamics};
    use crate::time::{Epoch, TimeUnits};
    use crate::Spacecraft;
    use crate::{io::ExportCfg, md::prelude::Traj, Orbit};
    use anise::almanac::Almanac;
    use anise::constants::frames::MOON_J2000;
    use pretty_env_logger;
    use std::env;
    use std::str::FromStr;
    use std::sync::Arc;
    use std::{collections::HashMap, path::PathBuf};

    #[test]
    fn test_load_oem_leo() {
        // All three samples were taken from https://github.com/bradsease/oem/blob/main/tests/samples/real/
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "data",
            "03_tests",
            "ccsds",
            "oem",
            "LEO_10s.oem",
        ]
        .iter()
        .collect();

        let _ = pretty_env_logger::try_init();

        let traj: Traj<Spacecraft> = Traj::from_oem_file(path, None).unwrap();

        // This trajectory has two duplicate epochs, which should be removed by the call to finalize()
        assert_eq!(traj.states.len(), 361);
        assert_eq!(traj.name.unwrap(), "0000-000A".to_string());
    }

    #[test]
    fn test_load_oem_meo() {
        // All three samples were taken from https://github.com/bradsease/oem/blob/main/tests/samples/real/
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "data",
            "03_tests",
            "ccsds",
            "oem",
            "MEO_60s.oem",
        ]
        .iter()
        .collect();

        let _ = pretty_env_logger::try_init();

        let traj: Traj<Spacecraft> = Traj::from_oem_file(path, None).unwrap();

        assert_eq!(traj.states.len(), 61);
        assert_eq!(traj.name.unwrap(), "0000-000A".to_string());
    }

    #[test]
    fn test_load_oem_geo() {
        use pretty_env_logger;
        use std::env;

        // All three samples were taken from https://github.com/bradsease/oem/blob/main/tests/samples/real/
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "data",
            "03_tests",
            "ccsds",
            "oem",
            "GEO_20s.oem",
        ]
        .iter()
        .collect();

        let _ = pretty_env_logger::try_init();

        let traj: Traj<Spacecraft> = Traj::from_oem_file(path, None).unwrap();

        assert_eq!(traj.states.len(), 181);
        assert_eq!(traj.name.as_ref().unwrap(), &"0000-000A".to_string());

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
            "data",
            "04_output",
            "GEO_20s_rebuilt.oem",
        ]
        .iter()
        .collect();

        traj.to_oem_file(
            &path,
            "0000-000A".to_string(),
            Some("Test Suite".to_string()),
            Some("TEST_OBJ".to_string()),
            cfg,
        )
        .unwrap();
        // And reload, make sure we have the same data.
        let traj_reloaded: Traj<Spacecraft> = Traj::from_oem_file(&path, None).unwrap();

        assert_eq!(traj_reloaded, traj);

        // Now export after trimming one state on either end
        let cfg = ExportCfg::builder()
            .timestamp(true)
            .metadata(HashMap::from([
                ("originator".to_string(), "Test suite".to_string()),
                ("object_name".to_string(), "TEST_OBJ".to_string()),
            ]))
            .step(20.seconds())
            .start_epoch(traj.first().orbit.epoch + 1.seconds())
            .end_epoch(traj.last().orbit.epoch - 1.seconds())
            .build();
        traj.to_oem_file(
            &path,
            "TEST-OBJ-ID".to_string(),
            Some("Test Suite".to_string()),
            Some("TEST_OBJ".to_string()),
            cfg,
        )
        .unwrap();
        // And reload, make sure we have the same data.
        let traj_reloaded: Traj<Spacecraft> = Traj::from_oem_file(path, None).unwrap();

        // Note that the number of states has changed because we interpolated with a step similar to the original one but
        // we started with a different time.
        assert_eq!(traj_reloaded.states.len(), traj.states.len() - 1);
        assert_eq!(
            traj_reloaded.first().orbit.epoch,
            traj.first().orbit.epoch + 1.seconds()
        );
        // Note: because we used a fixed step, the last epoch is actually an offset of step size - end offset
        // from the original trajectory
        assert_eq!(
            traj_reloaded.last().orbit.epoch,
            traj.last().orbit.epoch - 19.seconds()
        );
    }

    #[test]
    fn test_moon_frame_long_prop() {
        use std::path::PathBuf;

        let manifest_dir =
            PathBuf::from(std::env::var("CARGO_MANIFEST_DIR").unwrap_or(".".to_string()));

        let almanac = Almanac::new(
            &manifest_dir
                .clone()
                .join("data/01_planetary/pck08.pca")
                .to_string_lossy(),
        )
        .unwrap()
        .load(
            &manifest_dir
                .join("data/01_planetary/de440s.bsp")
                .to_string_lossy(),
        )
        .unwrap();

        let epoch = Epoch::from_str("2022-06-13T12:00:00").unwrap();
        let orbit = Orbit::try_keplerian_altitude(
            350.0,
            0.02,
            30.0,
            45.0,
            85.0,
            0.0,
            epoch,
            almanac.frame_info(MOON_J2000).unwrap(),
        )
        .unwrap();

        let mut traj =
            Propagator::default_dp78(SpacecraftDynamics::new(OrbitalDynamics::two_body()))
                .with(orbit.into(), Arc::new(almanac))
                .for_duration_with_traj(45.days())
                .unwrap()
                .1;
        // Set the name of this object
        traj.name = Some("TEST_MOON_OBJ".to_string());

        // Export CCSDS OEM file
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "data",
            "04_output",
            "moon_45days.oem",
        ]
        .iter()
        .collect();

        traj.to_oem_file(
            &path,
            "TEST_MOON_OBJ".to_string(),
            Some("Test Suite".to_string()),
            Some("TEST_MOON_OBJ".to_string()),
            ExportCfg::default(),
        )
        .unwrap();

        // And reload
        let traj_reloaded: Traj<Spacecraft> = Traj::from_oem_file(path, None).unwrap();

        assert_eq!(traj, traj_reloaded);
    }
}
