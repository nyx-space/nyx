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

pub(crate) const INTERPOLATION_SAMPLES: usize = 13;

use super::StateParameter;
use crate::cosmic::Frame;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::polyfit::hermite::hermite_eval;
use crate::time::Epoch;
use crate::{Orbit, Spacecraft, State};

use enum_iterator::all;

/// States that can be interpolated should implement this trait.
pub trait Interpolatable: State
where
    Self: Sized,
    DefaultAllocator: Allocator<f64, Self::Size>
        + Allocator<f64, Self::Size, Self::Size>
        + Allocator<f64, Self::VecLength>,
{
    /// Interpolates a new state at the provided epochs given a slice of states.
    fn interpolate(self, epoch: Epoch, states: &[Self]) -> Self;

    /// Returns the frame of this state
    fn frame(&self) -> Frame;

    /// Sets the frame of this state
    fn set_frame(&mut self, frame: Frame);

    /// List of state parameters that will be exported to a trajectory file in addition to the epoch (provided in this different formats).
    fn export_params() -> Vec<StateParameter>;

    /// Returns the orbit
    fn orbit(&self) -> &Orbit;
}

impl Interpolatable for Orbit {
    fn interpolate(self, epoch: Epoch, states: &[Self]) -> Self {
        // This is copied from ANISE

        // Statically allocated arrays of the maximum number of samples
        let mut epochs_tdb = [0.0; INTERPOLATION_SAMPLES + 1];
        let mut xs = [0.0; INTERPOLATION_SAMPLES + 1];
        let mut ys = [0.0; INTERPOLATION_SAMPLES + 1];
        let mut zs = [0.0; INTERPOLATION_SAMPLES + 1];
        let mut vxs = [0.0; INTERPOLATION_SAMPLES + 1];
        let mut vys = [0.0; INTERPOLATION_SAMPLES + 1];
        let mut vzs = [0.0; INTERPOLATION_SAMPLES + 1];

        for (cno, state) in states.iter().enumerate() {
            xs[cno] = state.x_km;
            ys[cno] = state.y_km;
            zs[cno] = state.z_km;
            vxs[cno] = state.vx_km_s;
            vys[cno] = state.vy_km_s;
            vzs[cno] = state.vz_km_s;
            epochs_tdb[cno] = state.epoch().to_tdb_seconds();
        }

        // TODO: Once I switch to using ANISE, this should use the same function as ANISE and not a clone. -- https://github.com/nyx-space/nyx/issues/86
        let (x_km, vx_km_s) = hermite_eval(
            &epochs_tdb[..states.len()],
            &xs[..states.len()],
            &vxs[..states.len()],
            epoch.to_et_seconds(),
        )
        .unwrap();

        let (y_km, vy_km_s) = hermite_eval(
            &epochs_tdb[..states.len()],
            &ys[..states.len()],
            &vys[..states.len()],
            epoch.to_et_seconds(),
        )
        .unwrap();

        let (z_km, vz_km_s) = hermite_eval(
            &epochs_tdb[..states.len()],
            &zs[..states.len()],
            &vzs[..states.len()],
            epoch.to_et_seconds(),
        )
        .unwrap();

        // And build the result
        let mut me = self;
        me.x_km = x_km;
        me.y_km = y_km;
        me.z_km = z_km;
        me.vx_km_s = vx_km_s;
        me.vy_km_s = vy_km_s;
        me.vz_km_s = vz_km_s;
        me.set_epoch(epoch);

        me
    }

    fn frame(&self) -> Frame {
        self.frame
    }

    fn set_frame(&mut self, frame: Frame) {
        self.frame = frame;
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
                            | StateParameter::GeodeticHeight
                            | StateParameter::GeodeticLatitude
                            | StateParameter::GeodeticLongitude
                    )
            })
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
        ]
        .concat()
    }

    fn orbit(&self) -> &Orbit {
        self
    }
}

impl Interpolatable for Spacecraft {
    fn interpolate(self, epoch: Epoch, states: &[Self]) -> Self {
        // Use the Orbit interpolation first.
        let orbit = Orbit::interpolate(
            self.orbit,
            epoch,
            &states.iter().map(|sc| sc.orbit).collect::<Vec<_>>(),
        );

        // Fuel is linearly interpolated -- should really be a Lagrange interpolation here
        let fuel_kg_dt = (states.last().unwrap().fuel_mass_kg
            - states.first().unwrap().fuel_mass_kg)
            / (states.last().unwrap().epoch().to_tdb_seconds()
                - states.first().unwrap().epoch().to_tdb_seconds());

        let mut me = self.with_orbit(orbit);
        me.fuel_mass_kg += fuel_kg_dt
            * (epoch.to_tdb_seconds() - states.first().unwrap().epoch().to_tdb_seconds());

        me
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
                            | StateParameter::GeodeticHeight
                            | StateParameter::GeodeticLatitude
                            | StateParameter::GeodeticLongitude
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
