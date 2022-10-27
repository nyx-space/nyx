/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

pub use super::Objective;
use super::{
    unit_vector_from_plane_angles, Frame, GuidanceLaw, GuidanceMode, NyxError, Orbit, Spacecraft,
    Vector3,
};
pub use crate::md::StateParameter;
use crate::utils::between_pm_180;
use crate::State;
use std::f64::consts::FRAC_PI_2 as half_pi;
use std::fmt;
use std::sync::Arc;

/// QLaw defines the Petropoulos refined laws from AAS/AIAA Space Flight Mechanics Conference, 2005
#[derive(Copy, Clone, Debug, Default)]
pub struct QLaw {
    /// Stores the objectives
    pub objectives: [Option<Objective>; 5],
    /// Stores the minimum efficiency to correct a given orbital element, defaults to zero (i.e. always correct)
    pub ηthresholds: [f64; 5],
    /// The "best quadratic time-to-go" weight, nominally 1 (for now, all W_oe = 1)
    pub w_p: f64,
    // init_state: Orbit,
}

impl QLaw {
    /// Creates a new QLaw locally optimal control as an Arc
    /// Note: this returns an Arc so it can be plugged into the Spacecraft dynamics directly.
    pub fn new(objectives: &[Objective], initial: Orbit) -> Result<Arc<Self>, NyxError> {
        Self::with_ηthresholds(objectives, &[0.0; 5], initial)
    }

    /// Creates a new QLaw locally optimal control as an Arc
    /// Note: this returns an Arc so it can be plugged into the Spacecraft dynamics directly.
    pub fn with_ηthresholds(
        objectives: &[Objective],
        ηthresholds: &[f64],
        _initial: Orbit,
    ) -> Result<Arc<Self>, NyxError> {
        let mut objs: [Option<Objective>; 5] = [None, None, None, None, None];
        let mut eff: [f64; 5] = [0.0; 5];
        if objectives.len() > 5 || objectives.is_empty() {
            return Err(NyxError::GuidanceConfigError(format!(
                "Must provide between 1 and 5 objectives (included), provided {}",
                objectives.len()
            )));
        } else if objectives.len() > ηthresholds.len() {
            return Err(NyxError::GuidanceConfigError(format!(
                "Must provide at least {} efficiency threshold values, provided {}",
                objectives.len(),
                ηthresholds.len()
            )));
        }

        for (i, obj) in objectives.iter().enumerate() {
            if [
                StateParameter::SMA,
                StateParameter::Eccentricity,
                StateParameter::Inclination,
                StateParameter::RAAN,
                StateParameter::AoP,
            ]
            .contains(&obj.parameter)
            {
                objs[i] = Some(*obj);
            } else {
                return Err(NyxError::GuidanceConfigError(format!(
                    "Objective {} not supported in Ruggerio",
                    obj.parameter
                )));
            }
        }
        for i in 0..objectives.len() {
            objs[i] = Some(objectives[i]);
            eff[i] = ηthresholds[i];
        }
        Ok(Arc::new(Self {
            objectives: objs,
            ηthresholds: eff,
            w_p: 1.0,
        }))
    }

    /// Penalty function
    pub fn penalty(osc: &Orbit) -> f64 {
        1.0 + (1.0 - osc.periapsis() / (150.0 + osc.frame.equatorial_radius())).exp()
    }

    /// Returns the distance of a given orbital element
    fn distance(&self, obj: &Objective, osc: &Orbit) -> f64 {
        let tgt_val = obj.desired_value;
        let tgt_tol = obj.tolerance;

        // NOTE: This function will modulo the angle errors +/- 180 deg, but paper recommends 0-180
        let dist = match obj.parameter {
            StateParameter::SMA => osc.sma() - tgt_val,
            StateParameter::Inclination => osc.inc() - tgt_val,
            StateParameter::Eccentricity => osc.ecc() - tgt_val,
            StateParameter::AoP => (between_pm_180(osc.aop() - tgt_val).to_radians())
                .cos()
                .acos(),
            StateParameter::RAAN => (between_pm_180(osc.raan() - tgt_val).to_radians())
                .cos()
                .acos(),
            _ => unreachable!(),
        };

        if (dist - tgt_tol).abs() < f64::EPSILON {
            0.0
        } else {
            dist
        }
    }

    fn weighting(&self, obj: &Objective, osc: &Orbit, _η_threshold: f64) -> f64 {
        let s = if obj.parameter == StateParameter::SMA {
            (1.0 + ((osc.sma() - obj.desired_value) / (3.0 * obj.desired_value)).powi(4)).powf(0.5)
        } else {
            1.0
        };
        let oe_xx = match obj.parameter {
            StateParameter::SMA => {
                2.0 * ((osc.sma().powi(3) * (1.0 + osc.ecc()))
                    / (osc.frame.gm() * (1.0 - osc.ecc())))
                .sqrt()
            }
            StateParameter::Eccentricity => 2.0 * osc.semi_parameter() / osc.hmag(),
            _ => unreachable!(),
        };
        s * (self.distance(obj, osc) / oe_xx).powi(2)
    }
}

impl fmt::Display for QLaw {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "QLaw with {} objectives", self.objectives.len())
    }
}

impl GuidanceLaw<GuidanceMode> for QLaw {
    /// Returns whether the guidance law has achieved all goals
    fn achieved(&self, state: &Spacecraft) -> Result<bool, NyxError> {
        for obj in self.objectives.iter().flatten() {
            if !obj.assess_raw(state.orbit.value(&obj.parameter)?).0 {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn direction(&self, sc: &Spacecraft) -> Vector3<f64> {
        if sc.mode() == GuidanceMode::Thrust {
            let osc = sc.orbit;
            let mut steering = Vector3::zeros();
            for (i, obj) in self.objectives.iter().flatten().enumerate() {
                let weight = self.weighting(obj, &osc, self.ηthresholds[i]);

                match obj.parameter {
                    StateParameter::SMA => {
                        let num = osc.ecc() * osc.ta().to_radians().sin();
                        let denom = 1.0 + osc.ecc() * osc.ta().to_radians().cos();
                        let alpha = num.atan2(denom);
                        // For SMA, we must multiply the weight by the thrust acceleration magnitude
                        steering += unit_vector_from_plane_angles(alpha, 0.0)
                            * weight
                            * sc.thruster.unwrap().thrust_N;
                    }
                    StateParameter::Eccentricity => {
                        let num = osc.ta().to_radians().sin();
                        let denom = osc.ta().to_radians().cos() + osc.ea().to_radians().cos();
                        let alpha = num.atan2(denom);
                        steering += unit_vector_from_plane_angles(alpha, 0.0) * weight;
                    }
                    StateParameter::Inclination => {
                        let beta = half_pi.copysign(((osc.ta() + osc.aop()).to_radians()).cos());
                        steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                    }
                    StateParameter::RAAN => {
                        let beta = half_pi.copysign(((osc.ta() + osc.aop()).to_radians()).sin());
                        steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                    }
                    StateParameter::AoP => {
                        let oe2 = 1.0 - osc.ecc().powi(2);
                        let e3 = osc.ecc().powi(3);
                        // Compute the optimal true anomaly for in-plane thrusting
                        let sqrt_val = (0.25 * (oe2 / e3).powi(2) + 1.0 / 27.0).sqrt();
                        let opti_ta_alpha = ((oe2 / (2.0 * e3) + sqrt_val).powf(1.0 / 3.0)
                            - (-oe2 / (2.0 * e3) + sqrt_val).powf(1.0 / 3.0)
                            - 1.0 / osc.ecc())
                        .acos();
                        // Compute the optimal true anomaly for out of plane thrusting
                        let opti_ta_beta = (-osc.ecc() * osc.aop().to_radians().cos()).acos()
                            - osc.aop().to_radians();
                        // And choose whether to do an in-plane or out of plane thrust
                        if (osc.ta().to_radians() - opti_ta_alpha).abs()
                            < (osc.ta().to_radians() - opti_ta_beta).abs()
                        {
                            // In plane
                            let p = osc.semi_parameter();
                            let (sin_ta, cos_ta) = osc.ta().to_radians().sin_cos();
                            let alpha = (-p * cos_ta).atan2((p + osc.rmag()) * sin_ta);
                            steering += unit_vector_from_plane_angles(alpha, 0.0) * weight;
                        } else {
                            // Out of plane
                            let beta = half_pi
                                .copysign(-(osc.ta().to_radians() + osc.aop().to_radians()).sin())
                                * osc.inc().to_radians().cos();
                            steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                        };
                    }
                    _ => unreachable!(),
                }
            }

            // Add the penalty factor
            steering *= Self::penalty(&osc);

            // Return a normalized vector
            steering = if steering.norm() > 0.0 {
                steering / steering.norm()
            } else {
                steering
            };
            // Convert to inertial -- this whole guidance law is computed in the RCN frame
            osc.dcm_from_traj_frame(Frame::RCN).unwrap() * steering
        } else {
            Vector3::zeros()
        }
    }

    // Either thrust full power or not at all
    fn throttle(&self, sc: &Spacecraft) -> f64 {
        if sc.mode() == GuidanceMode::Thrust {
            let osc = sc.orbit;
            for (i, obj) in self.objectives.iter().flatten().enumerate() {
                let weight = self.weighting(obj, &osc, self.ηthresholds[i]);
                if weight.abs() > 0.0 {
                    return 1.0;
                }
            }
            0.0
        } else {
            0.0
        }
    }

    /// Update the state for the next iteration
    fn next(&self, sc: &mut Spacecraft) {
        if sc.mode() != GuidanceMode::Inhibit {
            if !self.achieved(sc).unwrap() {
                if sc.mode() == GuidanceMode::Coast {
                    info!("enabling steering: {:x}", sc.orbit);
                }
                sc.mut_mode(GuidanceMode::Thrust);
            } else {
                if sc.mode() == GuidanceMode::Thrust {
                    info!("disabling steering: {:x}", sc.orbit);
                }
                sc.mut_mode(GuidanceMode::Coast);
            }
        }
    }
}
