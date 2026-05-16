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

use super::{GuidanceError, GuidanceLaw};
use crate::cosmic::{GuidanceMode, Spacecraft};
use crate::linalg::Vector3;
use crate::md::Trajectory;
use crate::State;
use anise::prelude::Almanac;
use std::fmt;
use std::sync::Arc;

/// Replays an inertial thrust-direction profile stored on a spacecraft trajectory.
///
/// The command is held piecewise-constant between stored trajectory samples.
#[derive(Clone)]
pub struct ThrustDirectionReplay {
    pub profile: Trajectory,
}

impl ThrustDirectionReplay {
    pub fn from_trajectory(profile: Trajectory) -> Arc<Self> {
        Arc::new(Self { profile })
    }

    fn command_at_or_before(&self, state: &Spacecraft) -> Option<&Spacecraft> {
        if self.profile.states.is_empty() {
            return None;
        }

        let epoch = state.epoch();
        if epoch < self.profile.first().epoch() || epoch > self.profile.last().epoch() {
            return None;
        }

        let command = match self
            .profile
            .states
            .binary_search_by(|sample| sample.epoch().cmp(&epoch))
        {
            Ok(idx) => self.profile.states.get(idx),
            Err(0) => None,
            Err(idx) => self.profile.states.get(idx - 1),
        };

        if let Some(command) = command {
            if command.thrust_direction().is_some() || command.mode() != GuidanceMode::Thrust {
                return Some(command);
            }
        }

        if self.profile.first().mode() == GuidanceMode::Thrust {
            return self
                .profile
                .states
                .iter()
                .find(|sample| sample.epoch() >= epoch && sample.thrust_direction().is_some());
        }

        command
    }
}

impl fmt::Display for ThrustDirectionReplay {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Thrust direction replay from {} to {} ({} samples)",
            self.profile.first().epoch(),
            self.profile.last().epoch(),
            self.profile.states.len()
        )
    }
}

impl GuidanceLaw for ThrustDirectionReplay {
    fn direction(&self, osc_state: &Spacecraft) -> Result<Vector3<f64>, GuidanceError> {
        Ok(self
            .command_at_or_before(osc_state)
            .and_then(|sample| sample.thrust_direction())
            .unwrap_or_else(Vector3::zeros))
    }

    fn throttle(&self, osc_state: &Spacecraft) -> Result<f64, GuidanceError> {
        Ok(
            if self
                .command_at_or_before(osc_state)
                .and_then(|sample| sample.thrust_direction())
                .is_some()
            {
                1.0
            } else {
                0.0
            },
        )
    }

    fn next(&self, next_state: &mut Spacecraft, _almanac: Arc<Almanac>) {
        let thrust_direction = self
            .command_at_or_before(next_state)
            .and_then(|sample| sample.thrust_direction());
        next_state.mut_thrust_direction(thrust_direction);
        next_state.mut_mode(if thrust_direction.is_some() {
            GuidanceMode::Thrust
        } else {
            GuidanceMode::Coast
        });
    }

    fn achieved(&self, osc_state: &Spacecraft) -> Result<bool, GuidanceError> {
        Ok(osc_state.epoch() > self.profile.last().epoch())
    }
}
