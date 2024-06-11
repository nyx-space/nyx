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

use anise::math::interpolation::InterpolationError;
use snafu::prelude::*;

mod interpolatable;
// mod orbit_traj;
mod sc_traj;
mod traj;
mod traj_it;

pub use interpolatable::Interpolatable;
pub(crate) use interpolatable::INTERPOLATION_SAMPLES;
pub use traj::Traj;

pub use crate::io::ExportCfg;

use super::StateParameter;
use crate::time::{Duration, Epoch};

#[derive(Clone, PartialEq, Debug, Snafu)]
pub enum TrajError {
    #[snafu(display("Event {event} not found between {start} and {end}"))]
    EventNotFound {
        start: Epoch,
        end: Epoch,
        event: String,
    },
    #[snafu(display("No interpolation data at {epoch}"))]
    NoInterpolationData { epoch: Epoch },
    #[snafu(display("Failed to create trajectory: {msg}"))]
    CreationError { msg: String },
    #[snafu(display("Probable bug: Requested epoch {req_epoch}, corresponding to an offset of {req_dur} in a spline of duration {spline_dur}"))]
    OutOfSpline {
        req_epoch: Epoch,
        req_dur: Duration,
        spline_dur: Duration,
    },
    #[snafu(display("Interpolation failed: {source}"))]
    Interpolation { source: InterpolationError },
}
