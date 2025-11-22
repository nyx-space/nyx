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

use anise::prelude::Almanac;
use log::debug;
use serde::{Deserialize, Serialize};
use snafu::ResultExt;

use super::{
    unit_vector_from_plane_angles, GuidStateSnafu, GuidanceError, GuidanceLaw, GuidanceMode,
    GuidancePhysicsSnafu, NyxError, Orbit, Spacecraft, Vector3,
};
use crate::cosmic::eclipse::EclipseLocator;
pub use crate::md::objective::Objective;
pub use crate::md::StateParameter;
use crate::State;
use std::f64::consts::FRAC_PI_2 as half_pi;
use std::fmt;
use std::sync::Arc;

/// Ruggiero defines the closed loop guidance law from IEPC 2011-102
#[derive(Copy, Clone, Default, Serialize, Deserialize)]
pub struct Ruggiero {
    /// Stores the objectives
    pub objectives: [Option<Objective>; 5],
    /// Stores the minimum efficiency to correct a given orbital element, defaults to zero (i.e. always correct)
    pub ηthresholds: [f64; 5],
    /// If defined, coast until vehicle is out of the provided eclipse state.
    pub max_eclipse_prct: Option<f64>,
    init_state: Spacecraft,
}

/// The Ruggiero is a locally optimal guidance law of a state for specific osculating elements.
/// NOTE: The efficiency parameters for AoP is NOT implemented: the paper's formulation is broken.
/// WARNING: Objectives must be in degrees!
impl Ruggiero {
    /// Creates a new Ruggiero locally optimal control as an Arc
    /// Note: this returns an Arc so it can be plugged into the Spacecraft dynamics directly.
    pub fn simple(objectives: &[Objective], initial: Spacecraft) -> Result<Arc<Self>, NyxError> {
        Self::from_ηthresholds(objectives, &[0.0; 5], initial)
    }

    /// Creates a new Ruggiero locally optimal control with the provided efficiency threshold.
    /// If the efficiency to correct the mapped orbital element is greater than the threshold, then the control law will be applied to this orbital element.
    /// Note: this returns an Arc so it can be plugged into the Spacecraft dynamics directly.
    pub fn from_ηthresholds(
        objectives: &[Objective],
        ηthresholds: &[f64],
        initial: Spacecraft,
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
            init_state: initial,
            ηthresholds: eff,
            max_eclipse_prct: None,
        }))
    }

    /// Creates a new Ruggiero locally optimal control as an Arc
    /// Note: this returns an Arc so it can be plugged into the Spacecraft dynamics directly.
    pub fn from_max_eclipse(
        objectives: &[Objective],
        initial: Spacecraft,
        max_eclipse: f64,
    ) -> Result<Arc<Self>, NyxError> {
        let mut objs: [Option<Objective>; 5] = [None, None, None, None, None];
        let eff: [f64; 5] = [0.0; 5];
        if objectives.len() > 5 || objectives.is_empty() {
            return Err(NyxError::GuidanceConfigError {
                msg: format!(
                    "Must provide between 1 and 5 objectives (included), provided {}",
                    objectives.len()
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
        }
        Ok(Arc::new(Self {
            objectives: objs,
            init_state: initial,
            ηthresholds: eff,
            max_eclipse_prct: Some(max_eclipse),
        }))
    }

    /// Sets the maximum eclipse during which we can thrust.
    pub fn set_max_eclipse(&mut self, max_eclipse: f64) {
        self.max_eclipse_prct = Some(max_eclipse);
    }

    /// Returns the efficiency η ∈ [0; 1] of correcting a specific orbital element at the provided osculating orbit
    pub fn efficiency(parameter: &StateParameter, osc_orbit: &Orbit) -> Result<f64, GuidanceError> {
        let e = osc_orbit.ecc().context(GuidancePhysicsSnafu {
            action: "computing Ruggiero efficiency",
        })?;

        let ν_ta = osc_orbit
            .ta_deg()
            .context(GuidancePhysicsSnafu {
                action: "computing Ruggiero efficiency",
            })?
            .to_radians();

        let ω = osc_orbit
            .aop_deg()
            .context(GuidancePhysicsSnafu {
                action: "computing Ruggiero efficiency",
            })?
            .to_radians();

        match parameter {
            StateParameter::SMA => {
                let a = osc_orbit.sma_km().context(GuidancePhysicsSnafu {
                    action: "computing Ruggiero efficiency",
                })?;

                let μ = osc_orbit.frame.mu_km3_s2().context(GuidancePhysicsSnafu {
                    action: "computing Ruggiero efficiency",
                })?;
                Ok(osc_orbit.vmag_km_s() * ((a * (1.0 - e)) / (μ * (1.0 + e))).sqrt())
            }
            StateParameter::Eccentricity => {
                let num = 1.0 + 2.0 * e * ν_ta.cos() + ν_ta.cos().powi(2);
                let denom = 1.0 + e * ν_ta.cos();
                // NOTE: There is a typo in IEPC 2011 102: the max of this efficiency function is at ν=0
                // where it is equal to 2*(2+2e) / (1+e). Therefore, I think the correct formulation should be
                // _divided_ by two, not multiplied by two.
                Ok(num / (2.0 * denom))
            }
            StateParameter::Inclination => {
                let num = (ω + ν_ta).cos().abs()
                    * ((1.0 - e.powi(2) * ω.sin().powi(2)).sqrt() - e * ω.cos().abs());
                let denom = 1.0 + e * ν_ta.cos();
                Ok(num / denom)
            }
            StateParameter::RAAN => {
                let num = (ω + ν_ta).sin().abs()
                    * ((1.0 - e.powi(2) * ω.cos().powi(2)).sqrt() - e * ω.sin().abs());
                let denom = 1.0 + e * ν_ta.cos();
                Ok(num / denom)
            }
            StateParameter::AoP => Ok(1.0),
            _ => Err(GuidanceError::InvalidControl { param: *parameter }),
        }
    }

    /// Computes the weight at which to correct this orbital element, will be zero if the current efficiency is below the threshold
    fn weighting(&self, obj: &Objective, osc_sc: &Spacecraft, η_threshold: f64) -> f64 {
        let init = self.init_state.value(obj.parameter).unwrap();
        let osc = osc_sc.value(obj.parameter).unwrap();
        let target = obj.desired_value;
        let tol = obj.tolerance;

        // Calculate the efficiency of correcting this specific orbital element
        let η = Self::efficiency(&obj.parameter, &osc_sc.orbit).unwrap();

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

    /// Returns whether the guidance law has achieved all goals
    pub fn status(&self, state: &Spacecraft) -> Vec<String> {
        self.objectives
            .iter()
            .flatten()
            .map(|obj| {
                let (ok, err) = obj.assess(state).unwrap();
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

impl fmt::Display for Ruggiero {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let obj_msg = self
            .objectives
            .iter()
            .flatten()
            .map(|obj| format!("{obj}"))
            .collect::<Vec<String>>();
        write!(
            f,
            "Ruggiero Controller (max eclipse: {}): \n {}",
            match self.max_eclipse_prct {
                Some(eclp) => format!("{eclp}"),
                None => "None".to_string(),
            },
            obj_msg.join("\n")
        )
    }
}

impl GuidanceLaw for Ruggiero {
    /// Returns whether the guidance law has achieved all goals
    fn achieved(&self, state: &Spacecraft) -> Result<bool, GuidanceError> {
        for obj in self.objectives.iter().flatten() {
            if !obj
                .assess_value(state.value(obj.parameter).context(GuidStateSnafu)?)
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
            let mut steering = Vector3::zeros();
            for (i, obj) in self.objectives.iter().flatten().enumerate() {
                let weight = self.weighting(obj, sc, self.ηthresholds[i]);
                if weight.abs() <= 0.0 {
                    continue;
                }

                // Compute all of the orbital elements here to unclutter the algorithm
                let ecc = osc.ecc().context(GuidancePhysicsSnafu {
                    action: "computing Ruggiero guidance",
                })?;

                let ta_rad = osc
                    .ta_deg()
                    .context(GuidancePhysicsSnafu {
                        action: "computing Ruggiero guidance",
                    })?
                    .to_radians();

                let inc_rad = osc
                    .inc_deg()
                    .context(GuidancePhysicsSnafu {
                        action: "computing Ruggiero guidance",
                    })?
                    .to_radians();

                let aop_rad = osc
                    .aop_deg()
                    .context(GuidancePhysicsSnafu {
                        action: "computing Ruggiero guidance",
                    })?
                    .to_radians();

                let ea_rad = osc
                    .ea_deg()
                    .context(GuidancePhysicsSnafu {
                        action: "computing Ruggiero guidance",
                    })?
                    .to_radians();

                match obj.parameter {
                    StateParameter::SMA => {
                        let num = ecc * ta_rad.sin();
                        let denom = 1.0 + ecc * ta_rad.cos();
                        let alpha = num.atan2(denom);
                        steering += unit_vector_from_plane_angles(alpha, 0.0) * weight;
                    }
                    StateParameter::Eccentricity => {
                        let num = ta_rad.sin();
                        let denom = ta_rad.cos() + ea_rad.cos();
                        let alpha = num.atan2(denom);
                        steering += unit_vector_from_plane_angles(alpha, 0.0) * weight;
                    }
                    StateParameter::Inclination => {
                        let beta = half_pi.copysign((ta_rad + aop_rad).cos());
                        steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                    }
                    StateParameter::RAAN => {
                        let beta = half_pi.copysign((ta_rad + aop_rad).sin());
                        steering += unit_vector_from_plane_angles(0.0, beta) * weight;
                    }
                    StateParameter::AoP => {
                        let oe2 = 1.0 - ecc.powi(2);
                        let e3 = ecc.powi(3);
                        // Compute the optimal true anomaly for in-plane thrusting
                        let sqrt_val = (0.25 * (oe2 / e3).powi(2) + 1.0 / 27.0).sqrt();
                        let opti_ta_alpha = ((oe2 / (2.0 * e3) + sqrt_val).powf(1.0 / 3.0)
                            - (-oe2 / (2.0 * e3) + sqrt_val).powf(1.0 / 3.0)
                            - 1.0 / ecc)
                            .acos();
                        // Compute the optimal true anomaly for out of plane thrusting
                        let opti_ta_beta = (-ecc * aop_rad.cos()).acos() - aop_rad;
                        // And choose whether to do an in-plane or out of plane thrust
                        if (ta_rad - opti_ta_alpha).abs() < (ta_rad - opti_ta_beta).abs() {
                            // In plane
                            let p = osc.semi_parameter_km().context(GuidancePhysicsSnafu {
                                action: "computing Ruggiero guidance",
                            })?;
                            let (sin_ta, cos_ta) = ta_rad.sin_cos();
                            let alpha = (-p * cos_ta).atan2((p + osc.rmag_km()) * sin_ta);
                            steering += unit_vector_from_plane_angles(alpha, 0.0) * weight;
                        } else {
                            // Out of plane
                            let beta = half_pi.copysign(-(ta_rad + aop_rad).sin()) * inc_rad.cos();
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
            if self.direction(sc)?.norm() > 0.0 {
                Ok(1.0)
            } else {
                Ok(0.0)
            }
        } else {
            Ok(0.0)
        }
    }

    /// Update the state for the next iteration
    fn next(&self, sc: &mut Spacecraft, almanac: Arc<Almanac>) {
        if sc.mode() != GuidanceMode::Inhibit {
            if !self.achieved(sc).unwrap() {
                // Check eclipse state if applicable.
                if let Some(max_eclipse) = self.max_eclipse_prct {
                    let locator = EclipseLocator::cislunar(almanac.clone());
                    if locator
                        .compute(sc.orbit, almanac)
                        .expect("cannot compute eclipse")
                        .percentage
                        > max_eclipse
                    {
                        // Coast in eclipse
                        sc.mode = GuidanceMode::Coast;
                    } else {
                        sc.mode = GuidanceMode::Thrust;
                    }
                } else if sc.mode() == GuidanceMode::Coast {
                    debug!("enabling steering: {:x}", sc.orbit);
                }
                sc.mut_mode(GuidanceMode::Thrust);
            } else {
                if sc.mode() == GuidanceMode::Thrust {
                    debug!("disabling steering: {:x}", sc.orbit);
                }
                sc.mut_mode(GuidanceMode::Coast);
            }
        }
    }
}

#[test]
fn ruggiero_weight() {
    use crate::time::Epoch;
    use anise::constants::frames::EARTH_J2000;

    let eme2k = EARTH_J2000.with_mu_km3_s2(398_600.433);
    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
    let orbit = Orbit::keplerian(7378.1363, 0.01, 0.05, 0.0, 0.0, 1.0, start_time, eme2k);
    let sc = Spacecraft::new(orbit, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);

    // Define the objectives
    let objectives = &[
        Objective::within_tolerance(StateParameter::SMA, 42164.0, 1.0),
        Objective::within_tolerance(StateParameter::Eccentricity, 0.01, 5e-5),
    ];

    let ruggiero = Ruggiero::simple(objectives, sc).unwrap();
    // 7301.597157 201.699933 0.176016 -0.202974 7.421233 0.006476 298.999726
    let osc = Orbit::new(
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

    let got = ruggiero.direction(&osc_sc).unwrap();

    assert!(
        dbg!(expected - got).norm() < 1e-12,
        "incorrect direction computed"
    );
}
