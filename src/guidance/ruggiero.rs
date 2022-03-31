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

use super::{
    unit_vector_from_plane_angles, Frame, GuidanceLaw, GuidanceMode, NyxError, Orbit, Spacecraft,
    Vector3,
};
pub use crate::md::objective::Objective;
pub use crate::md::StateParameter;
use crate::State;
use std::f64::consts::FRAC_PI_2 as half_pi;
use std::fmt;
use std::sync::Arc;

/// Ruggiero defines the closed loop guidance law from IEPC 2011-102
#[derive(Copy, Clone, Default, Debug)]
pub struct Ruggiero {
    /// Stores the objectives
    pub objectives: [Option<Objective>; 5],
    /// Stores the minimum efficiency to correct a given orbital element, defaults to zero (i.e. always correct)
    pub ηthresholds: [f64; 5],
    init_state: Orbit,
}

/// The Ruggiero is a locally optimal guidance law of a state for specific osculating elements.
/// NOTE: The efficency parameters for AoP is NOT implemented: the paper's formulation is broken.
/// WARNING: Objectives must be in degrees!
impl Ruggiero {
    /// Creates a new Ruggiero locally optimal control as an Arc
    /// Note: this returns an Arc so it can be plugged into the Spacecraft dynamics directly.
    pub fn new(objectives: &[Objective], initial: Orbit) -> Result<Arc<Self>, NyxError> {
        Self::with_ηthresholds(objectives, &[0.0; 5], initial)
    }

    /// Creates a new Ruggiero locally optimal control as an Arc
    /// Note: this returns an Arc so it can be plugged into the Spacecraft dynamics directly.
    pub fn with_ηthresholds(
        objectives: &[Objective],
        ηthresholds: &[f64],
        initial: Orbit,
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
            init_state: initial,
            ηthresholds: eff,
        }))
    }

    /// Returns the efficency η ∈ [0; 1] of correcting a specific orbital element at the provided osculating orbit
    pub fn efficency(parameter: &StateParameter, osc_orbit: &Orbit) -> Result<f64, NyxError> {
        let e = osc_orbit.ecc();
        match parameter {
            StateParameter::SMA => {
                let a = osc_orbit.sma();
                let μ = osc_orbit.frame.gm();
                Ok(osc_orbit.vmag() * ((a * (1.0 - e)) / (μ * (1.0 + e))).sqrt())
            }
            StateParameter::Eccentricity => {
                let ν_ta = osc_orbit.ta().to_radians();
                let num = 1.0 + 2.0 * e * ν_ta.cos() + ν_ta.cos().powi(2);
                let denom = 1.0 + e * ν_ta.cos();
                // NOTE: There is a typo in IEPC 2011 102: the max of this efficiency function is at ν=0
                // where it is equal to 2*(2+2e) / (1+e). Therefore, I think the correct formulation should be
                // _divided_ by two, not multiplied by two.
                Ok(num / (2.0 * denom))
            }
            StateParameter::Inclination => {
                let ν_ta = osc_orbit.ta().to_radians();
                let ω = osc_orbit.aop().to_radians();
                let num = (ω + ν_ta).cos().abs()
                    * ((1.0 - e.powi(2) * ω.sin().powi(2)).sqrt() - e * ω.cos().abs());
                let denom = 1.0 + e * ν_ta.cos();
                Ok(num / denom)
            }
            StateParameter::RAAN => {
                let ν_ta = osc_orbit.ta().to_radians();
                let ω = osc_orbit.aop().to_radians();
                let num = (ω + ν_ta).sin().abs()
                    * ((1.0 - e.powi(2) * ω.cos().powi(2)).sqrt() - e * ω.sin().abs());
                let denom = 1.0 + e * ν_ta.cos();
                Ok(num / denom)
            }
            StateParameter::AoP => Ok(1.0),
            _ => Err(NyxError::StateParameterUnavailable),
        }
    }

    /// Computes the weight at which to correct this orbital element, will be zero if the current efficency is below the threshold
    fn weighting(&self, obj: &Objective, osc_orbit: &Orbit, η_threshold: f64) -> f64 {
        let init = self.init_state.value(&obj.parameter).unwrap();
        let osc = osc_orbit.value(&obj.parameter).unwrap();
        let target = obj.desired_value;
        let tol = obj.tolerance;

        // Calculate the efficiency of correcting this specific orbital element
        let η = Self::efficency(&obj.parameter, osc_orbit).unwrap();

        if (osc - target).abs() < tol || η < η_threshold {
            0.0
        } else {
            // Let's add the tolerance to the initial value if we want to keep a parameter fixed (i.e. target and initial are equal)
            (target - osc)
                / (target
                    - if (init - target).abs() < tol {
                        init + tol
                    } else {
                        init
                    })
                .abs()
        }
    }
}

impl fmt::Display for Ruggiero {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Ruggiero with {} objectives", self.objectives.len())
    }
}

impl GuidanceLaw<GuidanceMode> for Ruggiero {
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
                if weight.abs() <= 0.0 {
                    continue;
                }

                match obj.parameter {
                    StateParameter::SMA => {
                        let num = osc.ecc() * osc.ta().to_radians().sin();
                        let denom = 1.0 + osc.ecc() * osc.ta().to_radians().cos();
                        let alpha = num.atan2(denom);
                        steering += unit_vector_from_plane_angles(alpha, 0.0) * weight;
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

#[test]
fn ruggiero_weight() {
    use crate::cosmic::Cosm;
    use crate::time::Epoch;
    let mut cosm = Cosm::de438_raw();
    cosm.frame_mut_gm("EME2000", 398_600.433);
    let eme2k = cosm.frame("EME2000");
    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let orbit = Orbit::keplerian(7378.1363, 0.01, 0.05, 0.0, 0.0, 1.0, start_time, eme2k);

    // Define the objectives
    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42164.0, 1.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.01, 5e-5),
    ];

    let ruggiero = Ruggiero::new(objectives, orbit).unwrap();
    // 7301.597157 201.699933 0.176016 -0.202974 7.421233 0.006476 298.999726
    let osc = Orbit::cartesian(
        7_303.253_461_441_64f64,
        127.478_714_816_381_75,
        0.111_246_193_227_445_4,
        -0.128_284_025_765_195_6,
        7.422_889_151_816_439,
        0.006_477_694_429_837_2,
        start_time,
        eme2k,
    );

    let mut osc_sc = Spacecraft::new(osc, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    // Must set the guidance mode to thrusting otherwise the direction will be set to zero.
    osc_sc.mut_mode(GuidanceMode::Thrust);

    let expected = Vector3::new(
        -0.017_279_636_133_108_3,
        0.999_850_315_226_803,
        0.000_872_534_222_883_2,
    );

    let got = ruggiero.direction(&osc_sc);

    assert!(
        dbg!(expected - got).norm() < 1e-12,
        "incorrect direction computed"
    );
}
