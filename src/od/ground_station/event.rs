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

use super::GroundStation;
use crate::md::EventEvaluator;
use crate::{errors::EventError, md::prelude::Interpolatable};
use anise::prelude::Almanac;
use hifitime::{Duration, Unit};
use nalgebra::{allocator::Allocator, DefaultAllocator};
use std::sync::Arc;

impl<S: Interpolatable> EventEvaluator<S> for &GroundStation
where
    DefaultAllocator: Allocator<S::Size> + Allocator<S::Size, S::Size> + Allocator<S::VecLength>,
{
    /// Compute the elevation in the SEZ frame. This call will panic if the frame of the input state does not match that of the ground station.
    fn eval(&self, rx_gs_frame: &S, almanac: Arc<Almanac>) -> Result<f64, EventError> {
        let dt = rx_gs_frame.epoch();
        // Then, compute the rotation matrix from the body fixed frame of the ground station to its topocentric frame SEZ.
        let tx_gs_frame = self.to_orbit(dt, &almanac).unwrap();

        let from = tx_gs_frame.frame.orientation_id * 1_000 + 1;
        let dcm_topo2fixed = tx_gs_frame
            .dcm_from_topocentric_to_body_fixed(from)
            .unwrap()
            .transpose();

        // Now, rotate the spacecraft in the SEZ frame to compute its elevation as seen from the ground station.
        // We transpose the DCM so that it's the fixed to topocentric rotation.
        let rx_sez = (dcm_topo2fixed * rx_gs_frame.orbit()).unwrap();
        let tx_sez = (dcm_topo2fixed * tx_gs_frame).unwrap();
        // Now, let's compute the range ρ.
        let rho_sez = (rx_sez - tx_sez).unwrap();

        // Finally, compute the elevation (math is the same as declination)
        // Source: Vallado, section 4.4.3
        // Only the sine is needed as per Vallado, and the formula is the same as the declination
        // because we're in the SEZ frame.
        Ok(rho_sez.declination_deg() - self.elevation_mask_deg)
    }

    fn eval_string(&self, state: &S, almanac: Arc<Almanac>) -> Result<String, EventError> {
        Ok(format!(
            "Elevation from {} is {:.6} deg on {}",
            self.name,
            self.eval(state, almanac)? + self.elevation_mask_deg,
            state.epoch()
        ))
    }

    /// Epoch precision of the election evaluator is 1 ms
    fn epoch_precision(&self) -> Duration {
        1 * Unit::Second
    }

    /// Angle precision of the elevation evaluator is 1 millidegree.
    fn value_precision(&self) -> f64 {
        1e-3
    }
}
