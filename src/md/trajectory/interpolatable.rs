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

use anise::math::interpolation::{hermite_eval, InterpolationError};

pub(crate) const INTERPOLATION_SAMPLES: usize = 13;

use super::StateParameter;
use crate::cosmic::Frame;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::time::Epoch;
use crate::{Orbit, Spacecraft, State};

use enum_iterator::all;

/// States that can be interpolated should implement this trait.
pub trait Interpolatable: State
where
    Self: Sized,
    DefaultAllocator: Allocator< Self::Size>
        + Allocator< Self::Size, Self::Size>
        + Allocator< Self::VecLength>,
{
    /// Interpolates a new state at the provided epochs given a slice of states.
    fn interpolate(self, epoch: Epoch, states: &[Self]) -> Result<Self, InterpolationError>;

    /// Returns the frame of this state
    fn frame(&self) -> Frame;

    /// Sets the frame of this state
    fn set_frame(&mut self, frame: Frame);

    /// List of state parameters that will be exported to a trajectory file in addition to the epoch (provided in this different formats).
    fn export_params() -> Vec<StateParameter>;

    /// Returns the orbit
    fn orbit(&self) -> &Orbit;
}

impl Interpolatable for Spacecraft {
    fn interpolate(mut self, epoch: Epoch, states: &[Self]) -> Result<Self, InterpolationError> {
        // Interpolate the Orbit first
        // Statically allocated arrays of the maximum number of samples
        let mut epochs_tdb = [0.0; INTERPOLATION_SAMPLES];
        let mut xs = [0.0; INTERPOLATION_SAMPLES];
        let mut ys = [0.0; INTERPOLATION_SAMPLES];
        let mut zs = [0.0; INTERPOLATION_SAMPLES];
        let mut vxs = [0.0; INTERPOLATION_SAMPLES];
        let mut vys = [0.0; INTERPOLATION_SAMPLES];
        let mut vzs = [0.0; INTERPOLATION_SAMPLES];

        for (cno, state) in states.iter().enumerate() {
            xs[cno] = state.orbit.radius_km.x;
            ys[cno] = state.orbit.radius_km.y;
            zs[cno] = state.orbit.radius_km.z;
            vxs[cno] = state.orbit.velocity_km_s.x;
            vys[cno] = state.orbit.velocity_km_s.y;
            vzs[cno] = state.orbit.velocity_km_s.z;
            epochs_tdb[cno] = state.epoch().to_et_seconds();
        }

        // Ensure that if we don't have enough states, we only interpolate using what we have instead of INTERPOLATION_SAMPLES
        let n = states.len();

        let (x_km, vx_km_s) =
            hermite_eval(&epochs_tdb[..n], &xs[..n], &vxs[..n], epoch.to_et_seconds())?;

        let (y_km, vy_km_s) =
            hermite_eval(&epochs_tdb[..n], &ys[..n], &vys[..n], epoch.to_et_seconds())?;

        let (z_km, vz_km_s) =
            hermite_eval(&epochs_tdb[..n], &zs[..n], &vzs[..n], epoch.to_et_seconds())?;

        self.orbit = Orbit::new(
            x_km,
            y_km,
            z_km,
            vx_km_s,
            vy_km_s,
            vz_km_s,
            epoch,
            self.orbit.frame,
        );

        // Fuel is linearly interpolated -- should really be a Lagrange interpolation here
        let first = states.first().unwrap();
        let last = states.last().unwrap();
        let fuel_kg_dt =
            (last.fuel_mass_kg - first.fuel_mass_kg) / (last.epoch() - first.epoch()).to_seconds();

        self.fuel_mass_kg += fuel_kg_dt * (epoch - first.epoch()).to_seconds();

        Ok(self)
    }

    fn frame(&self) -> Frame {
        self.orbit.frame
    }

    fn set_frame(&mut self, frame: Frame) {
        self.orbit.frame = frame;
    }

    fn export_params() -> Vec<StateParameter> {
        // Build all of the orbital parameters but keep the Cartesian state first
        let orbit_params = all::<StateParameter>()
            .filter(|p| {
                p.is_orbital()
                    && !p.is_b_plane()
                    && !matches!(
                        p,
                        StateParameter::X
                            | StateParameter::Y
                            | StateParameter::Z
                            | StateParameter::VX
                            | StateParameter::VY
                            | StateParameter::VZ
                            | StateParameter::HyperbolicAnomaly
                            | StateParameter::Height
                            | StateParameter::Latitude
                            | StateParameter::Longitude
                    )
            })
            .collect::<Vec<StateParameter>>();

        let sc_params = all::<StateParameter>()
            .filter(|p| p.is_for_spacecraft())
            .collect::<Vec<StateParameter>>();

        [
            vec![
                StateParameter::X,
                StateParameter::Y,
                StateParameter::Z,
                StateParameter::VX,
                StateParameter::VY,
                StateParameter::VZ,
            ],
            orbit_params,
            sc_params,
        ]
        .concat()
    }

    fn orbit(&self) -> &Orbit {
        &self.orbit
    }
}
