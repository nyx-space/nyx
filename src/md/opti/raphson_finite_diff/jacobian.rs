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
use crate::cosmic::{AstroAlmanacSnafu, AstroPhysicsSnafu};
use crate::dynamics::guidance::{Maneuver, MnvrRepr};
use crate::errors::TargetingError;
use crate::linalg::{SMatrix, Vector6};
use crate::md::prelude::*;
use crate::md::{AstroSnafu, StateParameter};
pub use crate::md::Vary;
use hifitime::TimeUnits;
use rayon::prelude::*;
use snafu::ResultExt;

use crate::md::objective::Objective;

#[allow(clippy::too_many_arguments)]
pub(crate) fn compute_jacobian<const V: usize, const O: usize>(
    targeter: &Targeter<V, O>,
    xi: &Spacecraft,
    mnvr: &Maneuver,
    finite_burn_target: bool,
    almanac: &Arc<Almanac>,
    achievement_epoch: Epoch,
    jac: &mut SMatrix<f64, O, V>,
    i: usize,
    obj: &Objective,
    achieved: f64,
    is_bplane_tgt: bool,
) -> Result<(), TargetingError> {
    let mut pert_calc: Vec<_> = targeter
        .variables
        .iter()
        .enumerate()
        .map(|(j, var)| (j, var, 0.0_f64))
        .collect();

    pert_calc.par_iter_mut().for_each(|(_, var, jac_val)| {
        let mut this_xi = *xi;

        let mut this_prop = targeter.prop.clone();
        let mut this_mnvr = *mnvr;

        let mut opposed_pert = false;

        if var.component.is_finite_burn() {
            // Modify the burn itself
            let pert = var.perturbation;
            // Modify the maneuver, but do not change the epochs of the maneuver unless the change is greater than one millisecond
            match var.component {
                Vary::Duration => {
                    if pert.abs() > 1e-3 {
                        this_mnvr.end = mnvr.start + pert.seconds()
                    }
                }
                Vary::EndEpoch => {
                    if pert.abs() > 1e-3 {
                        this_mnvr.end = mnvr.end + pert.seconds()
                    }
                }
                Vary::StartEpoch => {
                    if pert.abs() > 1e-3 {
                        this_mnvr.start = mnvr.start + pert.seconds()
                    }
                }
                Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => {
                    match mnvr.representation {
                        MnvrRepr::Angles { azimuth, elevation } => {
                            let azimuth = azimuth
                                .add_val_in_order(pert, var.component.vec_index())
                                .unwrap();
                            this_mnvr.representation =
                                MnvrRepr::Angles { azimuth, elevation };
                        }
                        _ => unreachable!(),
                    };
                }
                Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                    match mnvr.representation {
                        MnvrRepr::Angles { azimuth, elevation } => {
                            let elevation = elevation
                                .add_val_in_order(pert, var.component.vec_index())
                                .unwrap();
                            this_mnvr.representation =
                                MnvrRepr::Angles { azimuth, elevation };
                        }
                        _ => unreachable!(),
                    };
                }
                Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => {
                    let mut vector = this_mnvr.direction();
                    vector[var.component.vec_index()] += var.perturbation;
                    if !var.check_bounds(vector[var.component.vec_index()]).1 {
                        // Oops, bound was hit, go the other way
                        vector[var.component.vec_index()] -= 2.0 * var.perturbation;
                        opposed_pert = true;
                    }
                    this_mnvr.set_direction(vector).unwrap();
                }
                Vary::ThrustRateX | Vary::ThrustRateY | Vary::ThrustRateZ => {
                    let mut vector = this_mnvr.rate();
                    vector[(var.component.vec_index() - 1) % 3] += var.perturbation;
                    if !var
                        .check_bounds(vector[(var.component.vec_index() - 1) % 3])
                        .1
                    {
                        // Oops, bound was hit, go the other way
                        vector[(var.component.vec_index() - 1) % 3] -=
                            2.0 * var.perturbation;
                        opposed_pert = true;
                    }
                    this_mnvr.set_rate(vector).unwrap();
                }
                Vary::ThrustAccelX | Vary::ThrustAccelY | Vary::ThrustAccelZ => {
                    let mut vector = this_mnvr.accel();
                    vector[(var.component.vec_index() - 1) % 3] += var.perturbation;
                    if !var
                        .check_bounds(vector[(var.component.vec_index() - 1) % 3])
                        .1
                    {
                        // Oops, bound was hit, go the other way
                        vector[(var.component.vec_index() - 1) % 3] -=
                            2.0 * var.perturbation;
                        opposed_pert = true;
                    }
                    this_mnvr.set_accel(vector).unwrap();
                }
                Vary::ThrustLevel => {
                    this_mnvr.thrust_prct += var.perturbation;
                    this_mnvr.thrust_prct = this_mnvr.thrust_prct.clamp(0.0, 1.0);
                }
                _ => unreachable!(),
            }
        } else {
            let mut state_correction = Vector6::<f64>::zeros();
            state_correction[var.component.vec_index()] += var.perturbation;
            // Now, let's apply the correction to the initial state
            if let Some(frame) = targeter.correction_frame {
                // The following will error if the frame is not local
                let dcm_vnc2inertial = frame
                    .dcm_to_inertial(this_xi.orbit)
                    .context(AstroPhysicsSnafu)
                    .context(AstroSnafu)
                    .unwrap()
                    .rot_mat;

                let velocity_correction =
                    dcm_vnc2inertial * state_correction.fixed_rows::<3>(3);
                this_xi.orbit.apply_dv_km_s(velocity_correction);
            } else {
                this_xi = *xi + state_correction;
            }
        }

        let this_xf = if finite_burn_target {
            // Propagate normally until start of maneuver
            let pre_mnvr = this_prop
                .with(*xi, almanac.clone())
                .until_epoch(this_mnvr.start)
                .unwrap();
            // Add this maneuver to the dynamics, make sure that we don't over-step this maneuver
            let prop_opts = this_prop.opts;
            this_prop.set_max_step(this_mnvr.duration());
            this_prop.dynamics =
                this_prop.dynamics.with_guidance_law(Arc::new(this_mnvr));
            let post_mnvr = this_prop
                .with(
                    pre_mnvr.with_guidance_mode(GuidanceMode::Thrust),
                    almanac.clone(),
                )
                .until_epoch(this_mnvr.end)
                .unwrap();
            // Reset the propagator options to their previous configuration
            this_prop.opts = prop_opts;
            // And propagate until the achievement epoch
            this_prop
                .with(post_mnvr, almanac.clone())
                .until_epoch(achievement_epoch)
                .unwrap()
                .orbit
        } else {
            this_prop
                .with(this_xi, almanac.clone())
                .until_epoch(achievement_epoch)
                .unwrap()
                .orbit
        };

        let xf_dual_obj_frame = match &targeter.objective_frame {
            Some(frame) => {
                let orbit_obj_frame = almanac
                    .transform_to(this_xf, *frame, None)
                    .context(AstroAlmanacSnafu)
                    .context(AstroSnafu)
                    .unwrap();

                OrbitDual::from(orbit_obj_frame)
            }
            None => OrbitDual::from(this_xf),
        };

        let b_plane = if is_bplane_tgt {
            Some(BPlane::from_dual(xf_dual_obj_frame).unwrap())
        } else {
            None
        };

        let partial = if obj.parameter.is_b_plane() {
            match obj.parameter {
                StateParameter::BdotR => b_plane.unwrap().b_r,
                StateParameter::BdotT => b_plane.unwrap().b_t,
                StateParameter::BLTOF => b_plane.unwrap().ltof_s,
                _ => unreachable!(),
            }
        } else {
            xf_dual_obj_frame.partial_for(obj.parameter).unwrap()
        };

        let this_achieved = partial.real();
        *jac_val = (this_achieved - achieved) / var.perturbation;
        if opposed_pert {
            // We opposed the perturbation to ensure we don't over step a min/max bound
            *jac_val = -*jac_val;
        }
    });

    for (j, var, jac_val) in &pert_calc {
        // If this is a thrust level, we oppose the value so that the correction can still be positive.
        jac[(i, *j)] = if var.component == Vary::ThrustLevel {
            -*jac_val
        } else {
            *jac_val
        };
    }
    Ok(())
}
