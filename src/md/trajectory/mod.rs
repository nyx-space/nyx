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

mod interpolatable;
mod traj;
mod traj_it;

pub use interpolatable::Interpolatable;
pub(crate) use interpolatable::INTERPOLATION_SAMPLES;
use serde::{Deserialize, Serialize};
pub use traj::Traj;

use super::StateParameter;
use crate::time::{Duration, Epoch};

use std::collections::HashMap;
use std::error::Error;
use std::fmt;

#[derive(Clone, PartialEq, Eq, Debug)]
pub enum TrajError {
    EventNotFound {
        start: Epoch,
        end: Epoch,
        event: String,
    },
    NoInterpolationData(Epoch),
    CreationError(String),
    OutOfSpline {
        req_epoch: Epoch,
        req_dur: Duration,
        spline_dur: Duration,
    },
}

impl fmt::Display for TrajError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            Self::EventNotFound { start, end, event } => {
                write!(f, "Event {event} not found between {start} and {end}")
            }
            Self::CreationError(reason) => write!(f, "Failed to create trajectory: {reason}"),
            Self::NoInterpolationData(e) => write!(f, "No interpolation data at {e}"),
            Self::OutOfSpline {
                req_epoch,
                req_dur,
                spline_dur,
            } => {
                write!(f, "Probable bug: Requested epoch {req_epoch}, corresponding to an offset of {req_dur} in a spline of duration {spline_dur}")
            }
        }
    }
}

impl Error for TrajError {}

/// Configuration for exporting a trajectory to parquet.
#[derive(Clone, Default, Serialize, Deserialize)]
pub struct ExportCfg {
    /// Fields to export, if unset, defaults to all possible fields.
    pub fields: Option<Vec<StateParameter>>,
    /// Start epoch to export, defaults to the start of the trajectory
    pub start_epoch: Option<Epoch>,
    /// End epoch to export, defaults to the end of the trajectory
    pub end_epoch: Option<Epoch>,
    /// An optional step, defaults to every state in the trajectory (which likely isn't equidistant)
    pub step: Option<Duration>,
    /// Additional metadata to store in the Parquet metadata
    pub metadata: Option<HashMap<String, String>>,
    /// Set to true to append the timestamp to the filename
    pub timestamp: bool,
}

impl ExportCfg {
    /// Initialize a new configuration with the given metadata entries.
    pub fn from_metadata(metadata: Vec<(String, String)>) -> Self {
        let mut me = ExportCfg {
            metadata: Some(HashMap::new()),
            ..Default::default()
        };
        for (k, v) in metadata {
            me.metadata.as_mut().unwrap().insert(k, v);
        }
        me
    }

    pub fn append_field(&mut self, field: StateParameter) {
        if let Some(fields) = self.fields.as_mut() {
            fields.push(field);
        } else {
            self.fields = Some(vec![field]);
        }
    }

    pub fn set_step(&mut self, step: Duration) {
        self.step = Some(step);
    }
}
