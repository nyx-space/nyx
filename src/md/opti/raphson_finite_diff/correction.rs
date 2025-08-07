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

use crate::md::opti::targeter::Targeter;
use crate::dynamics::guidance::{Maneuver, MnvrRepr};
use crate::errors::TargetingError;
use crate::linalg::{SVector, Vector6};
use crate::md::prelude::*;
use crate::cosmic::AstroPhysicsSnafu;
use crate::md::{AstroSnafu, GuidanceSnafu};
pub use crate::md::Vary;
use hifitime::TimeUnits;
use snafu::ResultExt;

#[allow(clippy::too_many_arguments)]
pub(crate) fn apply_correction<const V: usize, const O: usize>(
    targeter: &Targeter<V, O>,
    delta: &SVector<f64, V>,
    xi: &mut Spacecraft,
    mnvr: &mut Maneuver,
    correction_epoch: Epoch,
    achievement_epoch: Epoch,
    is_initial_guess: bool,
) -> Result<(), TargetingError> {
    let mut state_correction = Vector6::<f64>::zeros();
    for (i, var) in targeter.variables.iter().enumerate() {
        let corr = if is_initial_guess {
            var.init_guess
        } else {
            delta[i]
        };

        if var.component.is_finite_burn() {
            // Modify the maneuver
            match var.component {
                Vary::Duration => {
                    if corr.abs() > 1e-3 {
                        let init_duration_s =
                            (correction_epoch - achievement_epoch).to_seconds();
                        let acceptable_corr = var.apply_bounds(init_duration_s).seconds();
                        mnvr.end = mnvr.start + acceptable_corr;
                    }
                }
                Vary::EndEpoch => {
                    if corr.abs() > 1e-3 {
                        let total_end_corr =
                            (mnvr.end + corr.seconds() - achievement_epoch).to_seconds();
                        let acceptable_corr = var.apply_bounds(total_end_corr).seconds();
                        mnvr.end += acceptable_corr;
                    }
                }
                Vary::StartEpoch => {
                    if corr.abs() > 1e-3 {
                        let total_start_corr =
                            (mnvr.start + corr.seconds() - correction_epoch).to_seconds();
                        let acceptable_corr = var.apply_bounds(total_start_corr).seconds();
                        mnvr.end += acceptable_corr;

                        mnvr.start += corr.seconds()
                    }
                }
                Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => {
                    match mnvr.representation {
                        MnvrRepr::Angles { azimuth, elevation } => {
                            let azimuth = azimuth
                                .add_val_in_order(corr, var.component.vec_index())
                                .unwrap();
                            mnvr.representation = MnvrRepr::Angles { azimuth, elevation };
                        }
                        _ => unreachable!(),
                    };
                }
                Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                    match mnvr.representation {
                        MnvrRepr::Angles { azimuth, elevation } => {
                            let elevation = elevation
                                .add_val_in_order(corr, var.component.vec_index())
                                .unwrap();
                            mnvr.representation = MnvrRepr::Angles { azimuth, elevation };
                        }
                        _ => unreachable!(),
                    };
                }
                Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => {
                    let mut vector = mnvr.direction();
                    vector[var.component.vec_index()] += corr;
                    var.ensure_bounds(&mut vector[var.component.vec_index()]);
                    mnvr.set_direction(vector).context(GuidanceSnafu)?;
                }
                Vary::ThrustRateX | Vary::ThrustRateY | Vary::ThrustRateZ => {
                    let mut vector = mnvr.rate();
                    let idx = (var.component.vec_index() - 1) % 3;
                    vector[idx] += corr;
                    var.ensure_bounds(&mut vector[idx]);
                    mnvr.set_rate(vector).context(GuidanceSnafu)?;
                }
                Vary::ThrustAccelX | Vary::ThrustAccelY | Vary::ThrustAccelZ => {
                    let mut vector = mnvr.accel();
                    let idx = (var.component.vec_index() - 1) % 3;
                    vector[idx] += corr;
                    var.ensure_bounds(&mut vector[idx]);
                    mnvr.set_accel(vector).context(GuidanceSnafu)?;
                }
                Vary::ThrustLevel => {
                    mnvr.thrust_prct += corr;
                    var.ensure_bounds(&mut mnvr.thrust_prct);
                }
                _ => unreachable!(),
            }
        } else {
            // State correction
            let mut current_corr = corr;
            if !is_initial_guess {
                // Choose the minimum step between the provided max step and the correction.
                if current_corr.abs() > var.max_step.abs() {
                    current_corr = var.max_step.abs() * current_corr.signum();
                } else if current_corr > var.max_value {
                    current_corr = var.max_value;
                } else if current_corr < var.min_value {
                    current_corr = var.min_value;
                }
            }
            state_correction[var.component.vec_index()] += current_corr;
        }
    }

    // Now, let's apply the correction to the initial state
    if let Some(frame) = targeter.correction_frame {
        let dcm_vnc2inertial = frame
            .dcm_to_inertial(xi.orbit)
            .context(AstroPhysicsSnafu)
            .context(AstroSnafu)?
            .rot_mat;

        let velocity_correction = dcm_vnc2inertial * state_correction.fixed_rows::<3>(3);
        xi.orbit.apply_dv_km_s(velocity_correction);
    } else {
        *xi = *xi + state_correction;
    }

    Ok(())
}
