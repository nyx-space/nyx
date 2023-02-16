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

pub use crate::dynamics::{Dynamics, NyxError};
pub use crate::{cosmic::Cosm, State, TimeTagged};
use hifitime::{Duration, Epoch};
use serde_derive::Deserialize;
use std::fmt::Debug;

pub mod arc;
pub mod trackdata;

#[derive(Debug, Deserialize)]
pub enum Schedule {
    Continuous,
    Intermittent { on: Duration, off: Duration },
}

pub enum StartMode {
    /// Start the tracking schedule with respect to when the receiver is visible
    Visible,
    /// Start the tracking schedule with respect to the provided epoch
    Epoch(Epoch),
}

/// Stores a tracking configuration, there is one per tracking data simulator (e.g. one for ground station #1 and another for #2)
pub struct TrkConfig {
    pub start_mode: StartMode,
    pub schedule: Schedule,
    pub sampling: Duration,
}
