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
use crate::cosmic::{Cosm, Frame, Orbit, Spacecraft};
use crate::errors::NyxError;
use crate::md::prelude::StateParameter;
use crate::md::EventEvaluator;
use crate::time::{Duration, TimeUnits};
use std::collections::HashMap;
use std::error::Error;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::time::Instant;

impl Traj<Spacecraft> {
    /// Allows converting the source trajectory into the (almost) equivalent trajectory in another frame
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
            traj.states
                .push(state.with_orbit(cosm.frame_chg(&state.orbit, new_frame)));
        }
        traj.finalize();

        info!(
            "Converted trajectory from {} to {} in {} ms: {traj}",
            self.first().orbit.frame,
            new_frame,
            (Instant::now() - start_instant).as_millis()
        );
        Ok(traj)
    }

    /// A shortcut to `to_parquet_with_cfg`
    pub fn to_parquet_with_step<P: AsRef<Path>>(
        &self,
        path: P,
        step: Duration,
    ) -> Result<(), Box<dyn Error>> {
        self.to_parquet_with_cfg(
            path,
            ExportCfg {
                step: Some(step),
                ..Default::default()
            },
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

    /// Convert this spacecraft trajectory into an Orbit trajectory, loosing all references to the spacecraft
    pub fn downcast(&self) -> Traj<Orbit> {
        let mut out = Traj::new();
        for sc_state in &self.states {
            out.states.push(sc_state.orbit);
        }
        out
    }

    /// Initialize a new spacecraft trajectory from the path to a CCSDS OEM file.
    ///
    /// CCSDS OEM only contains the orbit information, so you must provide a template spacecraft since we'll upcast the orbit trajectory into a spacecraft trajectory.
    pub fn from_oem_file<P: AsRef<Path>>(path: P, template: Spacecraft) -> Result<Self, NyxError> {
        let traj = Traj::<Orbit>::from_oem_file(path)?;

        Ok(traj.upcast(template))
    }
}
