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

use anise::math::{cartesian::CartesianState, interpolation::InterpolationError};
use snafu::prelude::*;

mod interpolatable;
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

/// Smooth the RIC differences using an in-line median filter.
/// This avoids allocations and operates directly on the Cartesian components.
fn smooth_state_diff_in_place(ric_diff: &mut [CartesianState], window_size: usize) {
    assert!(
        window_size % 2 == 1,
        "Window size must be odd for proper median calculation"
    );
    let half_window = window_size / 2;

    // Temporary buffer to store sorted values for median calculation
    let mut temp_buffer = vec![0.0; window_size];

    // Iterate over each state in the array
    for i in 0..ric_diff.len() {
        let start = i.saturating_sub(half_window);
        let end = (i + half_window + 1).min(ric_diff.len());

        // Smooth each component independently
        for component in 0..6 {
            // Fill the temporary buffer with values from the current window
            for (j, idx) in (start..end).enumerate() {
                temp_buffer[j] = match component {
                    0 => ric_diff[idx].radius_km.x,
                    1 => ric_diff[idx].radius_km.y,
                    2 => ric_diff[idx].radius_km.z,
                    3 => ric_diff[idx].velocity_km_s.x,
                    4 => ric_diff[idx].velocity_km_s.y,
                    5 => ric_diff[idx].velocity_km_s.z,
                    _ => unreachable!(),
                };
            }

            // Sort the buffer to find the median
            temp_buffer[..end - start].sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());

            // Replace the current value with the median
            let median = temp_buffer[(end - start) / 2];
            match component {
                0 => ric_diff[i].radius_km.x = median,
                1 => ric_diff[i].radius_km.y = median,
                2 => ric_diff[i].radius_km.z = median,
                3 => ric_diff[i].velocity_km_s.x = median,
                4 => ric_diff[i].velocity_km_s.y = median,
                5 => ric_diff[i].velocity_km_s.z = median,
                _ => unreachable!(),
            }
        }
    }
}
