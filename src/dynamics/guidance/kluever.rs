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

use anise::analysis::prelude::OrbitalElement;
use anise::prelude::Almanac;
use log::debug;
use serde::{Deserialize, Serialize};
use snafu::ResultExt;

use super::{
    GuidStateSnafu, GuidanceError, GuidanceLaw, GuidanceMode, GuidancePhysicsSnafu, Spacecraft,
    Vector3,
};
use crate::cosmic::eclipse::ShadowModel;
use crate::dynamics::guidance::ObjectiveWeight;
pub use crate::md::objective::Objective;
pub use crate::md::StateParameter;
use crate::State;
use std::fmt;
use std::sync::Arc;

/// Kluever defines a blended closed loop guidance law for low-thrust orbit transfers.
#[derive(Clone, Serialize, Deserialize)]
pub struct Kluever {
    /// Stores the objectives and their associated weights
    pub objectives: Vec<ObjectiveWeight>,
    /// If defined, coast until vehicle is out of the provided eclipse state.
    pub max_eclipse_prct: Option<f64>,
}

impl Kluever {
    /// Creates a new Kluever blended control law.
    pub fn new(objectives: &[Objective], weights: &[f64]) -> Arc<Self> {
        Arc::new(Self {
            objectives: objectives
                .iter()
                .copied()
                .zip(weights.iter().copied())
                .map(|(obj, w)| ObjectiveWeight {
                    objective: obj,
                    weight: w,
                })
                .collect(),
            max_eclipse_prct: None,
        })
    }

    /// Creates a new Kluever blended control law, specifying a maximum allowable eclipse for thrusting.
    pub fn from_max_eclipse(
        objectives: &[Objective],
        weights: &[f64],
        max_eclipse: f64,
    ) -> Arc<Self> {
        Arc::new(Self {
            objectives: objectives
                .iter()
                .copied()
                .zip(weights.iter().copied())
                .map(|(obj, w)| ObjectiveWeight {
                    objective: obj,
                    weight: w,
                })
                .collect(),
            max_eclipse_prct: Some(max_eclipse),
        })
    }

    /// Returns whether the guidance law has achieved all goals
    pub fn status(&self, state: &Spacecraft) -> Vec<String> {
        self.objectives
            .iter()
            .map(|obj| {
                let (ok, err) = obj.objective.assess(state).unwrap();
                format!(
                    "{} achieved: {}\t error = {:.5} {}",
                    obj.objective,
                    ok,
                    err,
                    obj.objective.parameter.unit()
                )
            })
            .collect::<Vec<String>>()
    }
}

impl fmt::Display for Kluever {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let obj_msg = self
            .objectives
            .iter()
            .map(|objective| {
                let obj = objective.objective;
                let weight = objective.weight;
                format!("{obj} (weight: {weight:.3})")
            })
            .collect::<Vec<String>>();
        write!(
            f,
            "Kluever Guidance (max eclipse: {}): \n {}",
            match self.max_eclipse_prct {
                Some(eclp) => format!("{eclp}"),
                None => "None".to_string(),
            },
            obj_msg.join("\n")
        )
    }
}

impl GuidanceLaw for Kluever {
    /// Returns whether the guidance law has achieved all goals
    fn achieved(&self, state: &Spacecraft) -> Result<bool, GuidanceError> {
        for obj in &self.objectives {
            if !obj
                .objective
                .assess_value(
                    state
                        .value(obj.objective.parameter)
                        .context(GuidStateSnafu)?,
                )
                .0
            {
                return Ok(false);
            }
        }
        Ok(true)
    }

    fn direction(&self, sc: &Spacecraft) -> Result<Vector3<f64>, GuidanceError> {
        if sc.mode() != GuidanceMode::Thrust {
            return Ok(Vector3::zeros());
        }

        let osc = sc.orbit;
        let mut sum_weighted_num_alpha = 0.0;
        let mut sum_weighted_den_alpha = 0.0;
        let mut sum_weighted_num_beta = 0.0;

        // Elements for calculations
        let ecc = osc.ecc().context(GuidancePhysicsSnafu {
            action: "Kluever guidance",
        })?;
        let ta_rad = osc
            .ta_deg()
            .context(GuidancePhysicsSnafu {
                action: "Kluever guidance",
            })?
            .to_radians();
        let aop_rad = osc
            .aop_deg()
            .context(GuidancePhysicsSnafu {
                action: "Kluever guidance",
            })?
            .to_radians();
        let raan_rad = osc
            .raan_deg()
            .context(GuidancePhysicsSnafu {
                action: "Kluever guidance",
            })?
            .to_radians();

        let u_rad = ta_rad + aop_rad;
        let l_rad = u_rad + raan_rad;
        let (sin_l, cos_l) = l_rad.sin_cos();

        for objective in &self.objectives {
            let obj = objective.objective;
            let weight = objective.weight;
            if weight == 0.0 {
                continue;
            }

            let osc_val = sc.value(obj.parameter).context(GuidStateSnafu)?;
            let error = obj.desired_value - osc_val;
            if error.abs() < obj.tolerance {
                continue;
            }
            let weight = weight * error.signum();

            match obj.parameter {
                StateParameter::Element(OrbitalElement::SemiMajorAxis) => {
                    // Maximize rate of change of energy
                    sum_weighted_num_alpha += weight * (ecc * ta_rad.sin());
                    sum_weighted_den_alpha += weight * (1.0 + ecc * ta_rad.cos());
                }
                StateParameter::Element(OrbitalElement::Eccentricity) => {
                    // Optimal alpha for eccentricity
                    let (sin_ta, cos_ta) = ta_rad.sin_cos();
                    sum_weighted_num_alpha += weight * sin_ta;
                    sum_weighted_den_alpha +=
                        weight * (cos_ta + (ecc + cos_ta) / (1.0 + ecc * cos_ta));
                }
                StateParameter::Element(OrbitalElement::Inclination) => {
                    // Purely out-of-plane requirement
                    let beta_opt = if u_rad.cos() >= 0.0 { 1.0 } else { -1.0 };
                    sum_weighted_num_beta += weight * beta_opt;
                }
                StateParameter::Element(OrbitalElement::RAAN) => {
                    let beta_opt = if u_rad.sin() >= 0.0 { 1.0 } else { -1.0 };
                    sum_weighted_num_beta += weight * beta_opt;
                }

                StateParameter::Element(OrbitalElement::EquinoctialH) => {
                    // H = e * sin(omega + RAAN)
                    let h = sc
                        .value(StateParameter::Element(OrbitalElement::EquinoctialH))
                        .context(GuidStateSnafu)?;
                    let k = sc
                        .value(StateParameter::Element(OrbitalElement::EquinoctialK))
                        .context(GuidStateSnafu)?;
                    sum_weighted_num_alpha += weight * cos_l;
                    sum_weighted_den_alpha +=
                        weight * (sin_l + (h + sin_l) / (1.0 + h * sin_l + k * cos_l));
                }

                StateParameter::Element(OrbitalElement::EquinoctialK) => {
                    // K = e * cos(omega + RAAN)
                    let h = sc
                        .value(StateParameter::Element(OrbitalElement::EquinoctialH))
                        .context(GuidStateSnafu)?;
                    let k = sc
                        .value(StateParameter::Element(OrbitalElement::EquinoctialK))
                        .context(GuidStateSnafu)?;
                    sum_weighted_num_alpha += weight * (-sin_l);
                    sum_weighted_den_alpha +=
                        weight * (cos_l + (k + cos_l) / (1.0 + h * sin_l + k * cos_l));
                }

                StateParameter::Element(OrbitalElement::EquinoctialP) => {
                    // P = tan(i/2) * sin(RAAN)
                    let beta_opt = if sin_l >= 0.0 { 1.0 } else { -1.0 };
                    sum_weighted_num_beta += weight * beta_opt;
                }

                StateParameter::Element(OrbitalElement::EquinoctialQ) => {
                    // Q = tan(i/2) * cos(RAAN)
                    let beta_opt = if cos_l >= 0.0 { 1.0 } else { -1.0 };
                    sum_weighted_num_beta += weight * beta_opt;
                }

                StateParameter::Element(OrbitalElement::EquinoctialLambda) => {
                    // Phasing / True Longitude
                    sum_weighted_den_alpha += weight * 1.0;
                }

                _ => {
                    return Err(GuidanceError::InvalidControl {
                        param: obj.parameter,
                    })
                }
            }
        }

        // Calculate blended steering angles
        let alpha = sum_weighted_num_alpha.atan2(sum_weighted_den_alpha);

        // Beta blending (simplified Kluever approach)
        let denom_beta = (sum_weighted_num_alpha.powi(2) + sum_weighted_den_alpha.powi(2)).sqrt();
        let beta = sum_weighted_num_beta.atan2(denom_beta);

        // Construct unit vector in RCN (Radial, Circum-normal, Normal) frame
        let (s_a, c_a) = alpha.sin_cos();
        let (s_b, c_b) = beta.sin_cos();
        let steering_rcn = Vector3::new(s_a * c_b, c_a * c_b, s_b);

        // Convert RCN to Inertial frame
        Ok(osc
            .dcm_from_rcn_to_inertial()
            .context(GuidancePhysicsSnafu {
                action: "RCN to Inertial",
            })?
            * steering_rcn)
    }

    /// Either thrust full power or not at all
    fn throttle(&self, sc: &Spacecraft) -> Result<f64, GuidanceError> {
        if sc.mode() == GuidanceMode::Thrust {
            Ok(1.0)
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
                    let locator = ShadowModel::cislunar(almanac.clone());
                    if locator
                        .compute(sc.orbit, almanac)
                        .expect("cannot compute eclipse")
                        .percentage
                        > max_eclipse
                    {
                        // Coast in eclipse
                        sc.mut_mode(GuidanceMode::Coast);
                    } else {
                        sc.mut_mode(GuidanceMode::Thrust);
                    }
                } else {
                    if sc.mode() == GuidanceMode::Coast {
                        debug!("enabling steering: {:x}", sc.orbit);
                    }
                    sc.mut_mode(GuidanceMode::Thrust);
                }
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
fn kluever_direction() {
    use crate::time::Epoch;
    use anise::analysis::prelude::OrbitalElement;
    use anise::constants::frames::EARTH_J2000;

    let eme2k = EARTH_J2000.with_mu_km3_s2(398_600.433);
    let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

    // Define the objectives
    let objectives = &[
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::SemiMajorAxis),
            42164.0,
            1.0,
        ),
        Objective::within_tolerance(
            StateParameter::Element(OrbitalElement::Eccentricity),
            0.01,
            5e-5,
        ),
    ];
    let weights = &[1.0, 1.0];

    let kluever = Kluever::new(objectives, weights);
    // 7301.597157 201.699933 0.176016 -0.202974 7.421233 0.006476 298.999726
    let osc = crate::Orbit::new(
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

    let got = kluever.direction(&osc_sc).unwrap();

    // Verification of the exact value might be tricky without a reference,
    // but we can check it's normalized and non-zero.
    assert!(got.norm() > 0.0);
    assert!((got.norm() - 1.0).abs() < 1e-12);

    println!("Kluever direction: {got}");
}
