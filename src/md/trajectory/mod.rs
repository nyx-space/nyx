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
pub use traj::Traj;

pub use crate::io::ExportCfg;

use super::StateParameter;
use crate::time::{Duration, Epoch};

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
