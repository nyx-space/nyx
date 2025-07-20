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

use crate::cosmic::AstroError;
use crate::errors::LambertError;
use crate::linalg::Vector3;
use std::f64::consts::PI;

mod godding;
mod izzo;

use anise::errors::PhysicsError;
use anise::prelude::Orbit;
pub use godding::gooding;
use hifitime::Duration;
pub use izzo::izzo;

const TAU: f64 = 2.0 * PI;
const LAMBERT_EPSILON: f64 = 1e-4; // General epsilon
const LAMBERT_EPSILON_TIME: f64 = 1e-4; // Time epsilon
const LAMBERT_EPSILON_RAD: f64 = (5e-5 / 180.0) * PI; // 0.00005 degrees
/// Maximum number of iterations allowed in the Lambert problem solver.
/// This is a safety measure to prevent infinite loops in case a solution cannot be found.
const MAX_ITERATIONS: usize = 1000;

/// Define the transfer kind for a Lambert
pub enum TransferKind {
    Auto,
    ShortWay,
    LongWay,
    NRevs(u8),
}

impl TransferKind {
    /// Calculate the direction multiplier based on the transfer kind.
    ///
    /// # Arguments
    ///
    /// * `r_final` - The final radius vector.
    /// * `r_init` - The initial radius vector.
    ///
    /// # Returns
    ///
    /// * `Result<f64, NyxError>` - The direction multiplier or an error if the transfer kind is not supported.
    fn direction_of_motion(
        self,
        r_final: &Vector3<f64>,
        r_init: &Vector3<f64>,
    ) -> Result<f64, LambertError> {
        match self {
            TransferKind::Auto => {
                let mut dnu = r_final[1].atan2(r_final[0]) - r_init[1].atan2(r_final[1]);
                if dnu > TAU {
                    dnu -= TAU;
                } else if dnu < 0.0 {
                    dnu += TAU;
                }

                if dnu > std::f64::consts::PI {
                    Ok(-1.0)
                } else {
                    Ok(1.0)
                }
            }
            TransferKind::ShortWay => Ok(1.0),
            TransferKind::LongWay => Ok(-1.0),
            _ => Err(LambertError::MultiRevNotSupported),
        }
    }
}

#[derive(Debug)]
pub struct LambertInput {
    pub r_init: Vector3<f64>,
    pub r_final: Vector3<f64>,
    pub tof: Duration,
    pub mu_km2_s3: f64,
}

impl LambertInput {
    /// Declination of the difference between the final and initial vector
    pub fn declination_deg(&self) -> f64 {
        let r_delta = self.r_final - self.r_init;
        (r_delta.z / r_delta.norm()).asin().to_degrees()
    }

    /// Right ascension of the difference between the final and initial vector
    pub fn right_ascension_deg(&self) -> f64 {
        let r_delta = self.r_final - self.r_init;
        r_delta.y.atan2(r_delta.x).to_degrees()
    }

    pub fn from_planetary_states(
        initial_state: Orbit,
        final_state: Orbit,
    ) -> Result<Self, AstroError> {
        if final_state.frame != initial_state.frame {
            return Err(AstroError::AstroPhysics {
                source: PhysicsError::FrameMismatch {
                    action: "Lambert solver requires both states to be in the same frame",
                    frame1: final_state.frame.into(),
                    frame2: initial_state.frame.into(),
                },
            });
        }
        Ok(Self {
            r_init: initial_state.radius_km,
            r_final: final_state.radius_km,
            tof: final_state.epoch - initial_state.epoch,
            mu_km2_s3: initial_state
                .frame
                .mu_km3_s2()
                .map_err(|e| AstroError::AstroPhysics { source: e })?,
        })
    }
}

#[derive(Debug)]
pub struct LambertSolution {
    pub v_init: Vector3<f64>,
    pub v_final: Vector3<f64>,
    pub phi: f64,
}

#[test]
fn test_lambert_vallado_lonway() {
    let ri = Vector3::new(15945.34, 0.0, 0.0);
    let rf = Vector3::new(12214.83899, 10249.46731, 0.0);
    let tof_s = 76.0 * 60.0;
    let mu_km3_s2 = 3.98600433e5;

    let exp_vi = Vector3::new(-3.811158, -2.003854, 0.0);
    let exp_vf = Vector3::new(4.207569, 0.914724, 0.0);

    let sol = gooding(ri, rf, tof_s, mu_km3_s2, TransferKind::LongWay).unwrap();

    assert!((sol.v_init - exp_vi).norm() < 1e-6);
    assert!((sol.v_final - exp_vf).norm() < 1e-6);
}
