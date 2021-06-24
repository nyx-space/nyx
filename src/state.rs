/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use crate::cosmic::{Frame, GuidanceMode, Orbit, Spacecraft, StmKind};
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{Const, DefaultAllocator, DimName, Matrix6, OMatrix, OVector, Vector1};
use crate::errors::NyxError;
use crate::time::{Duration, Epoch};
use std::fmt;
use std::ops::Add;

/// A trait allowing for something to have an epoch
pub trait TimeTagged {
    /// Retrieve the Epoch
    fn epoch(&self) -> Epoch;
    /// Set the Epoch
    fn set_epoch(&mut self, epoch: Epoch);

    /// Shift this epoch by a duration (can be negative)
    fn shift_by(&mut self, duration: Duration) {
        self.set_epoch(self.epoch() + duration);
    }
}

/// A trait for generate propagation and estimation state.
/// The first parameter is the size of the state, the second is the size of the propagated state including STM and extra items.
pub trait State:
    TimeTagged + Copy + Clone + PartialEq + fmt::Display + fmt::LowerExp + Send + Sync
where
    Self: Sized,
    DefaultAllocator: Allocator<f64, Self::Size>
        + Allocator<f64, Self::Size, Self::Size>
        + Allocator<f64, Self::VecLength>,
{
    /// Size of the state and its STM
    type Size: DimName;
    type VecLength: DimName;
    /// Initialize an empty state
    fn zeros() -> Self;

    /// Return this state as a vector for the propagation/estimation
    fn as_vector(&self) -> Result<OVector<f64, Self::VecLength>, NyxError>;

    /// Return this state as a vector for the propagation/estimation
    fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, NyxError>;

    /// Set this state
    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Self::VecLength>)
        -> Result<(), NyxError>;

    /// Reconstruct a new State from the provided delta time in seconds compared to the current state
    /// and with the provided vector.
    fn ctor_from(self, delta_t_s: f64, vector: &OVector<f64, Self::VecLength>) -> Self
    where
        DefaultAllocator: Allocator<f64, Self::VecLength>,
    {
        let mut me = self;
        me.set(me.epoch() + delta_t_s, vector).unwrap();
        me
    }

    fn add(self, other: OVector<f64, Self::Size>) -> Self;
}

/// Implementation of Orbit as a State for orbital dynamics with STM
impl State for Orbit {
    type Size = Const<6>;
    type VecLength = Const<42>;

    /// Returns a state whose position, velocity and frame are zero, and STM is I_{6x6}.
    fn zeros() -> Self {
        let frame = Frame::Celestial {
            gm: 1.0,
            ephem_path: [None, None, None],
            frame_path: [None, None, None],
        };

        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
            vx: 0.0,
            vy: 0.0,
            vz: 0.0,
            dt: Epoch::from_tai_seconds(0.0),
            frame,
            stm: Some(Matrix6::identity()),
            stm_kind: StmKind::Step,
        }
    }

    fn as_vector(&self) -> Result<OVector<f64, Const<42>>, NyxError> {
        let mut as_vec = OVector::<f64, Const<42>>::zeros();
        as_vec[0] = self.x;
        as_vec[1] = self.y;
        as_vec[2] = self.z;
        as_vec[3] = self.vx;
        as_vec[4] = self.vy;
        as_vec[5] = self.vz;
        let mut stm_idx = 6;
        if let Some(stm) = self.stm {
            for i in 0..6 {
                for j in 0..6 {
                    as_vec[stm_idx] = stm[(i, j)];
                    stm_idx += 1;
                }
            }
        }
        Ok(as_vec)
    }

    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Const<42>>) -> Result<(), NyxError> {
        self.set_epoch(epoch);
        self.x = vector[0];
        self.y = vector[1];
        self.z = vector[2];
        self.vx = vector[3];
        self.vy = vector[4];
        self.vz = vector[5];
        // And update the STM if applicable
        match self.stm_kind {
            StmKind::Step => {
                let mut stm_prev = self.stm.unwrap();
                let stm_k_to_0 = Matrix6::from_row_slice(&vector.as_slice()[6..]);

                if !stm_prev.try_inverse_mut() {
                    error!("STM not invertible: {}", stm_prev);
                    return Err(NyxError::SingularStateTransitionMatrix);
                }
                self.stm = Some(stm_k_to_0 * stm_prev);
            }
            StmKind::Traj => {
                let stm_k_to_0 = Matrix6::from_row_slice(&vector.as_slice()[6..]);
                self.stm = Some(stm_k_to_0)
            }
            StmKind::Unset => {}
        };
        Ok(())
    }

    fn stm(&self) -> Result<Matrix6<f64>, NyxError> {
        match self.stm {
            Some(stm) => Ok(stm),
            None => Err(NyxError::StateTransitionMatrixUnset),
        }
    }

    fn add(self, other: OVector<f64, Self::Size>) -> Self {
        self + other
    }
}

impl Add<OVector<f64, Const<6>>> for Orbit {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(self, other: OVector<f64, Const<6>>) -> Self {
        let mut me = self;
        me.x += other[0];
        me.y += other[1];
        me.z += other[2];
        me.vx += other[3];
        me.vy += other[4];
        me.vz += other[5];

        me
    }
}

impl TimeTagged for Spacecraft {
    fn epoch(&self) -> Epoch {
        self.orbit.dt
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.orbit.dt = epoch
    }
}

impl State for Spacecraft {
    type Size = Const<7>;
    type VecLength = Const<43>;

    fn zeros() -> Self {
        Self {
            orbit: Orbit::zeros(),
            dry_mass_kg: 0.0,
            fuel_mass_kg: 0.0,
            srp_area_m2: 0.0,
            drag_area_m2: 0.0,
            cr: 0.0,
            cd: 0.0,
            thruster: None,
            mode: GuidanceMode::Coast,
        }
    }

    fn as_vector(&self) -> Result<OVector<f64, Const<43>>, NyxError> {
        let orb_vec: OVector<f64, Const<42>> = self.orbit.as_vector()?;
        Ok(OVector::<f64, Const<43>>::from_iterator(
            orb_vec
                .iter()
                .chain(Vector1::new(self.fuel_mass_kg).iter())
                .cloned(),
        ))
    }

    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Const<43>>) -> Result<(), NyxError> {
        self.set_epoch(epoch);
        let orbit_vec = vector.fixed_rows::<42>(0).into_owned();
        self.orbit.set(epoch, &orbit_vec)?;
        self.fuel_mass_kg = vector[43 - 1];
        Ok(())
    }

    /// WARNING: Currently the STM assumes that the fuel mass is constant at ALL TIMES!
    fn stm(&self) -> Result<OMatrix<f64, Const<7>, Const<7>>, NyxError> {
        match self.orbit.stm {
            Some(stm) => {
                let mut rtn = OMatrix::<f64, Const<7>, Const<7>>::zeros();
                for i in 0..6 {
                    for j in 0..6 {
                        rtn[(i, j)] = stm[(i, j)];
                    }
                }
                rtn[(6, 6)] = 0.0;
                Ok(rtn)
            }
            None => Err(NyxError::StateTransitionMatrixUnset),
        }
    }

    fn add(self, other: OVector<f64, Self::Size>) -> Self {
        self + other
    }
}

impl Add<OVector<f64, Const<7>>> for Spacecraft {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(self, other: OVector<f64, Const<7>>) -> Self {
        let mut me = self;
        me.orbit.x += other[0];
        me.orbit.y += other[1];
        me.orbit.z += other[2];
        me.orbit.vx += other[3];
        me.orbit.vy += other[4];
        me.orbit.vz += other[5];
        me.fuel_mass_kg += other[6];

        me
    }
}

impl Add<OVector<f64, Const<6>>> for Spacecraft {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(self, other: OVector<f64, Const<6>>) -> Self {
        let mut me = self;
        me.orbit.x += other[0];
        me.orbit.y += other[1];
        me.orbit.z += other[2];
        me.orbit.vx += other[3];
        me.orbit.vy += other[4];
        me.orbit.vz += other[5];

        me
    }
}

#[test]
fn test_set_state() {
    let delta_t_s: f64 = 0.0;
    let mut new_vec = vec![0.0; 42];
    new_vec[0] = 1.0;
    new_vec[1] = 2.0;
    new_vec[2] = 3.0;
    new_vec[3] = 4.0;
    new_vec[4] = 5.0;
    new_vec[5] = 6.0;

    let state = OVector::<f64, Const<42>>::from_iterator(new_vec);
    let dummy_frame = Frame::Celestial {
        gm: 398600.4415,
        ephem_path: [None, None, None],
        frame_path: [None, None, None],
    };
    let ctx = Orbit::cartesian(
        6678.1363,
        2.0,
        3.0,
        4.0,
        7.725760634075587,
        6.0,
        Epoch::from_tai_days(0.0),
        dummy_frame,
    );
    let osc = ctx.ctor_from(delta_t_s, &state);
    println!("{}", osc);
}
