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
        let oe_xx = Self::oe_xx(osc, obj.parameter);
        s * (self.distance(obj, osc) / oe_xx).powi(2)
    }

    fn oe_xx(osc: &Orbit, param: StateParameter) -> f64 {
        match param {
            StateParameter::SMA => {
                2.0 * ((osc.sma().powi(3) * (1.0 + osc.ecc()))
                    / (osc.frame.gm() * (1.0 - osc.ecc())))
                .sqrt()
            }
            StateParameter::Eccentricity => 2.0 * osc.semi_parameter() / osc.hmag(),
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
    fn dq_doe(osc: &Orbit, param: StateParameter, thrust: f64) -> f64 {
        match param {
            StateParameter::SMA => {
                -osc.frame.gm() * (1.0 - osc.ecc())
                    / (2.0 * osc.sma().powi(3) * (1.0 + osc.ecc()) * thrust)
            }
            StateParameter::Eccentricity => {
                osc.ecc().powi(2) * osc.hmag()
                    / (osc.sma() * (1.0 - osc.ecc().powi(2)).powi(2) * thrust)
                    + osc.hmag() / (2.0 * osc.sma() * (1.0 - osc.ecc().powi(2))) * thrust
            }
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
    fn gaussian_vop_doe(osc: &Orbit, param: StateParameter) -> Vector3<f64> {
        match param {
            StateParameter::SMA => Vector3::new(
                2.0 * osc.sma().powi(2) * osc.ecc() * osc.ta().to_radians().sin() / osc.hmag(),
                2.0 * osc.sma().powi(3) * (1.0 - osc.ecc().powi(2)) / (osc.hmag() * osc.rmag()),
                0.0,
            ),
            StateParameter::Eccentricity => Vector3::new(
                osc.sma() * (1.0 - osc.ecc().powi(2)) * osc.ta().to_radians().sin() / osc.hmag(),
                osc.ecc() * osc.rmag()
                    + (osc.sma()
                        * (1.0 - osc.ecc().powi(2) + osc.rmag() * osc.ta().to_radians().cos()))
                        / osc.hmag(),
                0.0,
            ),
            _ => unreachable!(),
        }
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
            // let mut steering = Vector3::zeros();
            let mut d = Vector3::zeros();
            let mut q = 0.0;
            for (i, obj) in self.objectives.iter().flatten().enumerate() {
                q += self.weighting(obj, &osc, self.ηthresholds[i]);

                d += (Self::dq_doe(&osc, obj.parameter, sc.thruster.unwrap().thrust_N)
                    - osc.value(&obj.parameter).unwrap())
                    * Self::gaussian_vop_doe(&osc, obj.parameter);

                // match obj.parameter {
                //     StateParameter::SMA => {

                //         let num = osc.ecc() * osc.ta().to_radians().sin();
                //         let denom = 1.0 + osc.ecc() * osc.ta().to_radians().cos();
                //         let alpha = num.atan2(denom);
                //         // For SMA, we must multiply the weight by the thrust acceleration magnitude
                //         steering += unit_vector_from_plane_angles(alpha, 0.0)
                //             * weight
                //             * sc.thruster.unwrap().thrust_N;
                //     }
                //     StateParameter::Eccentricity => {
                //         let num = osc.ta().to_radians().sin();
                //         let denom = osc.ta().to_radians().cos() + osc.ea().to_radians().cos();
                //         let alpha = num.atan2(denom);
                //         steering += unit_vector_from_plane_angles(alpha, 0.0) * weight;
                //     }
                //     StateParameter::Inclination => {
                //         let beta = half_pi.copysign(((osc.ta() + osc.aop()).to_radians()).cos());
                //         steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                //     }
                //     StateParameter::RAAN => {
                //         let beta = half_pi.copysign(((osc.ta() + osc.aop()).to_radians()).sin());
                //         steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                //     }
                //     StateParameter::AoP => {
                //         let oe2 = 1.0 - osc.ecc().powi(2);
                //         let e3 = osc.ecc().powi(3);
                //         // Compute the optimal true anomaly for in-plane thrusting
                //         let sqrt_val = (0.25 * (oe2 / e3).powi(2) + 1.0 / 27.0).sqrt();
                //         let opti_ta_alpha = ((oe2 / (2.0 * e3) + sqrt_val).powf(1.0 / 3.0)
                //             - (-oe2 / (2.0 * e3) + sqrt_val).powf(1.0 / 3.0)
                //             - 1.0 / osc.ecc())
                //         .acos();
                //         // Compute the optimal true anomaly for out of plane thrusting
                //         let opti_ta_beta = (-osc.ecc() * osc.aop().to_radians().cos()).acos()
                //             - osc.aop().to_radians();
                //         // And choose whether to do an in-plane or out of plane thrust
                //         if (osc.ta().to_radians() - opti_ta_alpha).abs()
                //             < (osc.ta().to_radians() - opti_ta_beta).abs()
                //         {
                //             // In plane
                //             let p = osc.semi_parameter();
                //             let (sin_ta, cos_ta) = osc.ta().to_radians().sin_cos();
                //             let alpha = (-p * cos_ta).atan2((p + osc.rmag()) * sin_ta);
                //             steering += unit_vector_from_plane_angles(alpha, 0.0) * weight;
                //         } else {
                //             // Out of plane
                //             let beta = half_pi
                //                 .copysign(-(osc.ta().to_radians() + osc.aop().to_radians()).sin())
                //                 * osc.inc().to_radians().cos();
                //             steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                //         };
                //     }
                //     _ => unreachable!(),
                // }
            }

            // Add the penalty factor
            d *= Self::penalty(&osc);

            // Solve for the optimal angles (we swap the D1 and D2 because our function return the thrust in r, theta, h, not theta, r, h)
            let alpha = -d[0].atan2(-d[1]);
            let beta = (-d[2] / (d[1].powi(2) + d[0].powi(2)).sqrt()).atan();

            // dbg!(alpha, beta);
            dbg!(q);
            if q.abs() < 5.0e-3 {
                println!("{:x}", osc);
                panic!();
            }

            // // Return a normalized vector
            // steering = if steering.norm() > 0.0 {
            //     steering / steering.norm()
            // } else {
            //     steering
            // };
            let steering = unit_vector_from_plane_angles(alpha, beta);
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
