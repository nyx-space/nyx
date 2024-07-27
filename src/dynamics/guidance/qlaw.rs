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

use anise::astro::PhysicsResult;
use snafu::ResultExt;

pub use super::Objective;
use super::{
    unit_vector_from_plane_angles, GuidStateSnafu, GuidanceError, GuidanceLaw, GuidanceMode,
    NyxError, Orbit, Spacecraft, Vector3,
};
use crate::dynamics::guidance::GuidancePhysicsSnafu;
pub use crate::md::StateParameter;
use crate::utils::between_pm_180;
use crate::State;
// use std::f64::consts::FRAC_PI_2 as half_pi;
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
}

impl QLaw {
    /// Creates a new QLaw locally optimal control as an Arc
    /// Note: this returns an Arc so it can be plugged into the Spacecraft dynamics directly.
    pub fn new(objectives: &[Objective], initial: Spacecraft) -> Result<Arc<Self>, NyxError> {
        Self::with_ηthresholds(objectives, &[0.0; 5], initial)
    }

    /// Creates a new QLaw locally optimal control as an Arc
    /// Note: this returns an Arc so it can be plugged into the Spacecraft dynamics directly.
    pub fn with_ηthresholds(
        objectives: &[Objective],
        ηthresholds: &[f64],
        _initial: Spacecraft,
    ) -> Result<Arc<Self>, NyxError> {
        let mut objs: [Option<Objective>; 5] = [None, None, None, None, None];
        let mut eff: [f64; 5] = [0.0; 5];
        if objectives.len() > 5 || objectives.is_empty() {
            return Err(NyxError::GuidanceConfigError {
                msg: format!(
                    "Must provide between 1 and 5 objectives (included), provided {}",
                    objectives.len()
                ),
            });
        } else if objectives.len() > ηthresholds.len() {
            return Err(NyxError::GuidanceConfigError {
                msg: format!(
                    "Must provide at least {} efficiency threshold values, provided {}",
                    objectives.len(),
                    ηthresholds.len()
                ),
            });
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
                return Err(NyxError::GuidanceConfigError {
                    msg: format!("Objective {} not supported in Ruggerio", obj.parameter),
                });
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
    pub fn penalty(osc: &Orbit) -> PhysicsResult<f64> {
        Ok(1.0
            + (1.0 - osc.periapsis_km()? / (150.0 + osc.frame.mean_equatorial_radius_km()?)).exp())
    }

    /// Returns the distance of a given orbital element
    fn distance(&self, obj: &Objective, osc: &Orbit) -> PhysicsResult<f64> {
        let tgt_val = obj.desired_value;
        let tgt_tol = obj.tolerance;

        // NOTE: This function will modulo the angle errors +/- 180 deg, but paper recommends 0-180
        let dist = match obj.parameter {
            StateParameter::SMA => osc.sma_km()? - tgt_val,
            StateParameter::Inclination => osc.inc_deg()? - tgt_val,
            StateParameter::Eccentricity => osc.ecc()? - tgt_val,
            StateParameter::AoP => (between_pm_180(osc.aop_deg()? - tgt_val).to_radians())
                .cos()
                .acos(),
            StateParameter::RAAN => (between_pm_180(osc.raan_deg()? - tgt_val).to_radians())
                .cos()
                .acos(),
            _ => unreachable!(),
        };

        if (dist - tgt_tol).abs() < f64::EPSILON {
            Ok(0.0)
        } else {
            Ok(dist)
        }
    }

    fn weighting(&self, obj: &Objective, osc: &Orbit, _η_threshold: f64) -> PhysicsResult<f64> {
        let s = if obj.parameter == StateParameter::SMA {
            (1.0 + ((osc.sma_km()? - obj.desired_value) / (3.0 * obj.desired_value)).powi(4))
                .powf(0.5)
        } else {
            1.0
        };
        let oe_xx = Self::oe_xx(osc, obj.parameter)?;
        Ok(s * (self.distance(obj, osc)? / oe_xx).powi(2))
    }

    fn oe_xx(osc: &Orbit, param: StateParameter) -> PhysicsResult<f64> {
        match param {
            StateParameter::SMA => Ok(2.0
                * ((osc.sma_km()?.powi(3) * (1.0 + osc.ecc()?))
                    / (osc.frame.mu_km3_s2()? * (1.0 - osc.ecc()?)))
                .sqrt()),
            StateParameter::Eccentricity => Ok(2.0 * osc.semi_parameter_km()? / osc.hmag()?),
            _ => unreachable!(),
        }
    }

    /// Computes the partial derivative of Q with respect to the provided orbital element
    /// This was computed using Sympy, pasted here for information.
    ///
    /// # SMA
    /// ```python
    /// >>> ecc = Symbol('e'); sma = Symbol('a'); mu = Symbol('mu'); h = Symbol('h'); ta = Symbol('ta'); fr = Symbol('fr'); ft = Symbol('ft'); fh = Symbol('fh')
    /// >>> f = (fr**2 + ft **2 + fh**2)**0.5
    /// >>> sma_xx = 2 * f * sqrt( (sma**3 * (1+ecc)) / (mu * (1-ecc)))
    /// >>> (sma/sma_xx**2).diff(sma)
    /// -mu*(1 - e)/(2*a**3*(e + 1)*(fh**2 + fr**2 + ft**2)**1.0)
    /// ```
    ///
    /// # ECC
    /// ```python
    /// >>> p = sma * (1.0 - ecc**2)
    /// >>> ecc_xx = 2 * f * p / h
    /// >>> (ecc/ecc_xx).diff(ecc)
    /// e**2*h/(a*(1.0 - e**2)**2*(fh**2 + fr**2 + ft**2)**0.5) + h/(2*a*(1.0 - e**2)*(fh**2 + fr**2 + ft**2)**0.5
    /// ```
    fn dq_doe(osc: &Orbit, param: StateParameter, thrust: f64) -> PhysicsResult<f64> {
        match param {
            StateParameter::SMA => Ok(-osc.frame.mu_km3_s2()? * (1.0 - osc.ecc()?)
                / (2.0 * osc.sma_km()?.powi(3) * (1.0 + osc.ecc()?) * thrust)),
            StateParameter::Eccentricity => Ok(osc.ecc()?.powi(2) * osc.hmag()?
                / (osc.sma_km()? * (1.0 - osc.ecc()?.powi(2)).powi(2) * thrust)
                + osc.hmag()? / (2.0 * osc.sma_km()? * (1.0 - osc.ecc()?.powi(2))) * thrust),
            _ => unreachable!(),
        }
    }

    /// Compute the partial derivatives over the thrust vector for each orbital element.
    /// Returns doe/dfr, doe/dftheta, deo/dfh
    ///
    /// # SMA
    /// ```python
    /// In [14]: ta = Symbol('ta')
    /// In [15]: fr = Symbol('fr')
    /// In [16]: ft = Symbol('ft')
    /// In [17]: fh = Symbol('fh')
    /// In [19]: r = Symbol('r')
    /// In [20]: gvop_sma = (2*sma**2 / h) * (ecc * sin(ta) * fr + p/r * ft)
    /// In [21]: gvop_sma.diff(fr)
    /// Out[21]: 2*a**2*e*sin(ta)/h
    /// In [22]: gvop_sma.diff(ft)
    /// Out[22]: 2*a**3*(1.0 - e**2)/(h*r)
    /// In [23]: gvop_sma.diff(fh)
    /// Out[23]: 0
    /// ```
    ///
    /// # ECC
    /// ```python
    /// In [24]: gvop_ecc = 1/h * ( p*sin(ta) *fr + ( (p+r) * cos(ta) + r*ecc )* ft )
    /// In [25]: gvop_ecc.diff(fr)
    /// Out[25]: a*(1.0 - e**2)*sin(ta)/h
    /// In [26]: gvop_ecc.diff(ft)
    /// Out[26]: (e*r + (a*(1.0 - e**2) + r)*cos(ta))/h
    /// In [27]: gvop_ecc.diff(fh)
    /// Out[27]: 0
    /// ```
    fn gaussian_vop_doe(osc: &Orbit, param: StateParameter) -> PhysicsResult<Vector3<f64>> {
        match param {
            StateParameter::SMA => Ok(Vector3::new(
                2.0 * osc.sma_km()?.powi(2) * osc.ecc()? * osc.ta_deg()?.to_radians().sin()
                    / osc.hmag()?,
                2.0 * osc.sma_km()?.powi(3) * (1.0 - osc.ecc()?.powi(2))
                    / (osc.hmag()? * osc.rmag_km()),
                0.0,
            )),
            StateParameter::Eccentricity => Ok(Vector3::new(
                osc.sma_km()? * (1.0 - osc.ecc()?.powi(2)) * osc.ta_deg()?.to_radians().sin()
                    / osc.hmag()?,
                osc.ecc()? * osc.rmag_km()
                    + (osc.sma_km()?
                        * (1.0 - osc.ecc()?.powi(2)
                            + osc.rmag_km() * osc.ta_deg()?.to_radians().cos()))
                        / osc.hmag()?,
                0.0,
            )),
            _ => unreachable!(),
        }
    }

    /// Returns whether the guidance law has achieved all goals
    pub fn status(&self, state: &Spacecraft) -> Vec<String> {
        self.objectives
            .iter()
            .flatten()
            .map(|obj| {
                let (ok, err) = obj.assess_raw(state.value(obj.parameter).unwrap());
                format!(
                    "{} achieved: {}\t error = {:.5} {}",
                    obj,
                    ok,
                    err,
                    obj.parameter.unit()
                )
            })
            .collect::<Vec<String>>()
    }
}

impl fmt::Display for QLaw {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "QLaw with {} objectives", self.objectives.len())
    }
}

impl GuidanceLaw for QLaw {
    /// Returns whether the guidance law has achieved all goals
    fn achieved(&self, state: &Spacecraft) -> Result<bool, GuidanceError> {
        for obj in self.objectives.iter().flatten() {
            if !obj
                .assess_raw(state.value(obj.parameter).context(GuidStateSnafu)?)
                .0
            {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn direction(&self, sc: &Spacecraft) -> Result<Vector3<f64>, GuidanceError> {
        if sc.mode() == GuidanceMode::Thrust {
            let osc = sc.orbit;
            // let mut steering = Vector3::zeros();
            let mut d = Vector3::zeros();
            let mut q = 0.0;
            for (i, obj) in self.objectives.iter().flatten().enumerate() {
                q += self.weighting(obj, &osc, self.ηthresholds[i]).context(
                    GuidancePhysicsSnafu {
                        action: "computing Q-Law weights",
                    },
                )?;

                d += (Self::dq_doe(&osc, obj.parameter, sc.thruster.unwrap().thrust_N).context(
                    GuidancePhysicsSnafu {
                        action: "computing Q-Law partial",
                    },
                )? - sc.value(obj.parameter).unwrap())
                    * Self::gaussian_vop_doe(&osc, obj.parameter).context(
                        GuidancePhysicsSnafu {
                            action: "computing Q-Law partial",
                        },
                    )?;

                // match obj.parameter {
                //     StateParameter::SMA => {

                //         let num = osc.ecc()? * osc.ta_deg()?.to_radians().sin();
                //         let denom = 1.0 + osc.ecc()? * osc.ta_deg()?.to_radians().cos();
                //         let alpha = num.atan2(denom);
                //         // For SMA, we must multiply the weight by the thrust acceleration magnitude
                //         steering += unit_vector_from_plane_angles(alpha, 0.0)
                //             * weight
                //             * sc.thruster.unwrap().thrust_N;
                //     }
                //     StateParameter::Eccentricity => {
                //         let num = osc.ta_deg()?.to_radians().sin();
                //         let denom = osc.ta_deg()?.to_radians().cos() + osc.ea().to_radians().cos();
                //         let alpha = num.atan2(denom);
                //         steering += unit_vector_from_plane_angles(alpha, 0.0) * weight;
                //     }
                //     StateParameter::Inclination => {
                //         let beta = half_pi.copysign(((osc.ta() + osc.aop_deg()).to_radians()).cos());
                //         steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                //     }
                //     StateParameter::RAAN => {
                //         let beta = half_pi.copysign(((osc.ta() + osc.aop_deg()).to_radians()).sin());
                //         steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                //     }
                //     StateParameter::AoP => {
                //         let oe2 = 1.0 - osc.ecc()?.powi(2);
                //         let e3 = osc.ecc()?.powi(3);
                //         // Compute the optimal true anomaly for in-plane thrusting
                //         let sqrt_val = (0.25 * (oe2 / e3).powi(2) + 1.0 / 27.0).sqrt();
                //         let opti_ta_alpha = ((oe2 / (2.0 * e3) + sqrt_val).powf(1.0 / 3.0)
                //             - (-oe2 / (2.0 * e3) + sqrt_val).powf(1.0 / 3.0)
                //             - 1.0 / osc.ecc())
                //         .acos();
                //         // Compute the optimal true anomaly for out of plane thrusting
                //         let opti_ta_beta = (-osc.ecc()? * osc.aop_deg().to_radians().cos()).acos()
                //             - osc.aop_deg().to_radians();
                //         // And choose whether to do an in-plane or out of plane thrust
                //         if (osc.ta_deg()?.to_radians() - opti_ta_alpha).abs()
                //             < (osc.ta_deg()?.to_radians() - opti_ta_beta).abs()
                //         {
                //             // In plane
                //             let p = osc.semi_parameter();
                //             let (sin_ta, cos_ta) = osc.ta_deg()?.to_radians().sin_cos();
                //             let alpha = (-p * cos_ta).atan2((p + osc.rmag_km()) * sin_ta);
                //             steering += unit_vector_from_plane_angles(alpha, 0.0) * weight;
                //         } else {
                //             // Out of plane
                //             let beta = half_pi
                //                 .copysign(-(osc.ta_deg()?.to_radians() + osc.aop_deg().to_radians()).sin())
                //                 * osc.inc_deg().to_radians().cos();
                //             steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                //         };
                //     }
                //     _ => unreachable!(),
                // }
            }

            // Add the penalty factor
            d *= Self::penalty(&osc).context(GuidancePhysicsSnafu {
                action: "computing Q-Law penalty",
            })?;

            // Solve for the optimal angles (we swap the D1 and D2 because our function return the thrust in r, theta, h, not theta, r, h)
            let alpha = -d[0].atan2(-d[1]);
            let beta = (-d[2] / (d[1].powi(2) + d[0].powi(2)).sqrt()).atan();

            // dbg!(alpha, beta);
            // dbg!(q);
            // if q.abs() < 5.0e-3 {
            //     println!("{:x}", osc);
            //     panic!();
            // }

            let steering = unit_vector_from_plane_angles(alpha, beta);

            // Convert to inertial -- this whole guidance law is computed in the RCN frame
            Ok(osc
                .dcm_from_rcn_to_inertial()
                .context(GuidancePhysicsSnafu {
                    action: "computing RCN frame",
                })?
                * steering)
        } else {
            Ok(Vector3::zeros())
        }
    }

    // Either thrust full power or not at all
    fn throttle(&self, sc: &Spacecraft) -> Result<f64, GuidanceError> {
        if sc.mode() == GuidanceMode::Thrust {
            let osc = sc.orbit;
            for (i, obj) in self.objectives.iter().flatten().enumerate() {
                let weight = self.weighting(obj, &osc, self.ηthresholds[i]).context(
                    GuidancePhysicsSnafu {
                        action: "computing Q-Law weights",
                    },
                )?;
                if weight.abs() > 0.0 {
                    return Ok(1.0);
                }
            }
            Ok(0.0)
        } else {
            Ok(0.0)
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
