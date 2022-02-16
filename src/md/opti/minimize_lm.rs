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

extern crate levenberg_marquardt;

use super::optimizer::Optimizer;
// use super::solution::TargeterSolution;
use crate::dynamics::guidance::Mnvr;
use crate::errors::TargetingError;
use crate::linalg::{storage::Owned, Const, SMatrix, SVector, Vector6};
use crate::linalg::{DimMax, DimMin, ToTypenum};
use crate::md::rayon::prelude::*;
use crate::md::ui::*;
use crate::md::StateParameter;
pub use crate::md::{Variable, Vary};
use crate::polyfit::CommonPolynomial;
use crate::propagators::error_ctrl::ErrorCtrl;
use hifitime::TimeUnits;
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};
// use std::time::Instant;

/// N: number of variables; M: number of objectives
pub struct OptimizerInstance<'a, E: ErrorCtrl, const N: usize, const M: usize>
where
    Const<N>: ToTypenum,
    Const<M>: ToTypenum,
    Const<M>: DimMin<Const<N>, Output = Const<N>> + DimMax<Const<N>, Output = Const<N>>,
{
    /// The propagator setup (kind, stages, etc.)
    pub prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
    /// The list of objectives of this targeter
    pub objectives: [Objective; M],
    /// An optional frame (and Cosm) to compute the objectives in.
    /// Needed if the propagation frame is separate from objectives frame (e.g. for B Plane targeting).
    pub objective_frame: Option<(Frame, Arc<Cosm>)>,
    /// The kind of correction to apply to achieve the objectives
    pub variables: [Variable; N],
    /// The frame in which the correction should be applied, must be either a local frame or inertial
    pub correction_frame: Option<Frame>,
    /// The starting state of the optimizer
    pub spacecraft: Spacecraft,
    /// Epoch at which the objectives must be achieved
    /// TODO: Convert this to an Either<Epoch, Event> (and create an `until` function on the propagator that takes an Either as well).
    pub achievement_epoch: Epoch,
    /// Epoch at which the trajectory is corrected
    pub correction_epoch: Epoch,
    /// The control solution to this problem.
    pub control: SVector<f64, N>,
    pub residuals: SVector<f64, M>,
    // pub jacobian: SMatrix<f64, O, V>,
    pub jacobian: SMatrix<f64, M, N>,
}

// We implement a trait for every problem we want to solve
impl<'a, E: ErrorCtrl, const V: usize, const O: usize> LeastSquaresProblem<f64, Const<O>, Const<V>>
    for OptimizerInstance<'a, E, V, O>
where
    Const<V>: ToTypenum,
    Const<O>: ToTypenum,
    Const<O>: DimMin<Const<V>, Output = Const<V>> + DimMax<Const<V>, Output = Const<V>>,
{
    type ResidualStorage = Owned<f64, Const<O>>;
    type ParameterStorage = Owned<f64, Const<V>>;
    type JacobianStorage = Owned<f64, Const<O>, Const<V>>;

    fn set_params(&mut self, attempted_control: &SVector<f64, V>) {
        // TODO: Switch methods based on whether the finite differencing is requested
        // do common calculations for residuals and the Jacobian here

        println!("Ctrl: {}", attempted_control);
        self.control = *attempted_control;

        let mut is_bplane_tgt = false;
        for obj in &self.objectives {
            if obj.parameter.is_b_plane() {
                is_bplane_tgt = true;
                break;
            }
        }

        // Now we know that the problem is correctly defined, so let's propagate as is to the epoch
        // where the correction should be applied.
        let xi_start = self
            .prop
            .with(self.spacecraft)
            .until_epoch(self.correction_epoch)
            .unwrap();

        debug!("xi_start = {}", xi_start);

        let mut xi = xi_start;
        // We'll store the initial state correction here.
        let mut state_correction = Vector6::<f64>::zeros();

        // Store the total correction in Vector3
        let mut total_correction = SVector::<f64, V>::zeros();

        // Create a default maneuver that will only be used if a finite burn is being targeted
        let mut mnvr = Mnvr {
            start: self.correction_epoch,
            end: self.achievement_epoch,
            thrust_lvl: 1.0,
            alpha_inplane_radians: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
            delta_outofplane_radians: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
            frame: Frame::RCN,
        };

        let mut finite_burn_target = false;

        // Apply the initial guess
        for (i, var) in self.variables.iter().enumerate() {
            // Check the validity (this function will report to log and raise an error)
            var.valid().unwrap();
            // Check that there is no attempt to target a position in a local frame
            if self.correction_frame.is_some() && var.component.vec_index() < 3 {
                // Then this is a position correction, which is not allowed if a frame is provided!
                let msg = format!(
                    "Variable is in frame {} but that frame cannot be used for a {:?} correction",
                    self.correction_frame.unwrap(),
                    var.component
                );
                error!("{}", msg);
                panic!();
            }

            // Check that a thruster is provided since we'll be changing that and the burn duration
            if var.component.is_finite_burn() {
                if xi_start.thruster.is_none() {
                    // Can't do any conversion to finite burns without a thruster
                    // return Err(NyxError::CtrlExistsButNoThrusterAvail);
                    panic!();
                }
                finite_burn_target = true;
                // Modify the default maneuver
                match var.component {
                    Vary::Duration => mnvr.end = mnvr.start + attempted_control[i].seconds(),
                    Vary::EndEpoch => mnvr.end += attempted_control[i].seconds(),
                    Vary::StartEpoch => mnvr.start += attempted_control[i].seconds(),
                    Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => {
                        mnvr.alpha_inplane_radians = mnvr
                            .alpha_inplane_radians
                            .add_val_in_order(attempted_control[i], var.component.vec_index())
                            .unwrap();
                    }
                    Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                        mnvr.delta_outofplane_radians = mnvr
                            .delta_outofplane_radians
                            .add_val_in_order(attempted_control[i], var.component.vec_index())
                            .unwrap();
                    }
                    Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => {
                        let mut vector = mnvr.vector(mnvr.start);
                        vector[var.component.vec_index()] = attempted_control[i];
                        mnvr.set_direction(vector).unwrap();
                    }
                    Vary::ThrustLevel => {
                        mnvr.thrust_lvl += attempted_control[i];
                    }
                    _ => unreachable!(),
                }
                info!("Initial maneuver guess: {}", mnvr);
            } else {
                state_correction[var.component.vec_index()] -= attempted_control[i];
                // Now, let's apply the correction to the initial state
                if let Some(frame) = self.correction_frame {
                    // The following will error if the frame is not local
                    let dcm_vnc2inertial = xi.orbit.dcm_from_traj_frame(frame).unwrap();
                    let velocity_correction =
                        dcm_vnc2inertial * state_correction.fixed_rows::<3>(3);
                    xi.orbit.apply_dv(velocity_correction);
                } else {
                    xi.orbit.x += state_correction[0];
                    xi.orbit.y += state_correction[1];
                    xi.orbit.z += state_correction[2];
                    xi.orbit.vx += state_correction[3];
                    xi.orbit.vy += state_correction[4];
                    xi.orbit.vz += state_correction[5];
                }
            }

            total_correction[i] += attempted_control[i];
        }

        // Determine padding in debugging info
        // For the width, we find the largest desired values and multiply it by the order of magnitude of its tolerance
        let max_obj_val = self
            .objectives
            .iter()
            .map(|obj| {
                (obj.desired_value.abs().ceil() as i32
                    * 10_i32.pow(obj.tolerance.abs().log10().ceil() as u32)) as i32
            })
            .max()
            .unwrap();

        let max_obj_tol = self
            .objectives
            .iter()
            .map(|obj| obj.tolerance.log10().abs().ceil() as usize)
            .max()
            .unwrap();

        let width = f64::from(max_obj_val).log10() as usize + 2 + max_obj_tol;

        // Modify each variable by the desired perturbatino, propagate, compute the final parameter, and store how modifying that variable affects the final parameter
        let cur_xi = xi;

        // If we are targeting a finite burn, let's set propagate in several steps to make sure we don't miss the burn
        let xf = if finite_burn_target {
            info!("{}", mnvr);
            let mut prop = self.prop.clone();
            let prop_opts = prop.opts;
            let pre_mnvr = prop.with(cur_xi).until_epoch(mnvr.start).unwrap();
            prop.dynamics = prop.dynamics.with_guidance_law_no_decr(Arc::new(mnvr));
            prop.set_max_step(mnvr.end - mnvr.start);
            let post_mnvr = prop
                .with(pre_mnvr.with_guidance_mode(GuidanceMode::Thrust))
                .until_epoch(mnvr.end)
                .unwrap();
            // Reset the propagator options to their previous configuration
            prop.opts = prop_opts;
            // And propagate until the achievement epoch
            prop.with(post_mnvr)
                .until_epoch(self.achievement_epoch)
                .unwrap()
                .orbit
        } else {
            self.prop
                .with(cur_xi)
                .until_epoch(self.achievement_epoch)
                .unwrap()
                .orbit
        };

        let xf_dual_obj_frame = match &self.objective_frame {
            Some((frame, cosm)) => {
                let orbit_obj_frame = cosm.frame_chg(&xf, *frame);
                OrbitDual::from(orbit_obj_frame)
            }
            None => OrbitDual::from(xf),
        };

        // Build the B-Plane once, if needed, and always in the objective frame
        let b_plane = if is_bplane_tgt {
            Some(BPlane::from_dual(xf_dual_obj_frame).unwrap())
        } else {
            None
        };

        // Build debugging information
        let mut objmsg = Vec::with_capacity(self.objectives.len());

        // The Jacobian includes the sensitivity of each objective with respect to each variable for the whole trajectory.
        // As such, it includes the STM of that variable for the whole propagation arc.
        // let mut jac = DMatrix::from_element(self.objectives.len(), self.variables.len(), 0.0);

        for (i, obj) in self.objectives.iter().enumerate() {
            let partial = if obj.parameter.is_b_plane() {
                match obj.parameter {
                    StateParameter::BdotR => b_plane.unwrap().b_r,
                    StateParameter::BdotT => b_plane.unwrap().b_t,
                    StateParameter::BLTOF => b_plane.unwrap().ltof_s,
                    _ => unreachable!(),
                }
            } else {
                xf_dual_obj_frame.partial_for(&obj.parameter).unwrap()
            };

            let achieved = partial.real();

            self.residuals[i] = obj.assess_raw(achieved).1;

            objmsg.push(format!(
                "\t{:?}: achieved = {:>width$.prec$}\t desired = {:>width$.prec$}\t scaled error = {:>width$.prec$}",
                obj.parameter,
                achieved,
                obj.desired_value,
                self.residuals[i], width=width, prec=max_obj_tol
            ));

            let mut pert_calc: Vec<_> = self
                .variables
                .iter()
                .enumerate()
                .map(|(j, var)| (j, var, 0.0_f64))
                .collect();

            pert_calc.par_iter_mut().for_each(|(_, var, jac_val)| {
                let mut this_xi = xi;

                let mut this_prop = self.prop.clone();
                let mut this_mnvr = mnvr;

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
                            this_mnvr.alpha_inplane_radians = mnvr
                                .alpha_inplane_radians
                                .add_val_in_order(pert, var.component.vec_index())
                                .unwrap();
                        }
                        Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                            this_mnvr.delta_outofplane_radians = mnvr
                                .delta_outofplane_radians
                                .add_val_in_order(pert, var.component.vec_index())
                                .unwrap();
                        }
                        Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => {
                            let mut vector = this_mnvr.vector(self.correction_epoch);
                            vector[var.component.vec_index()] += pert;
                            this_mnvr.set_direction(vector).unwrap();
                        }
                        Vary::ThrustLevel => {
                            this_mnvr.thrust_lvl += pert;
                        }
                        _ => unreachable!(),
                    }
                } else {
                    let mut state_correction = Vector6::<f64>::zeros();
                    state_correction[var.component.vec_index()] += var.perturbation;
                    // Now, let's apply the correction to the initial state
                    if let Some(frame) = self.correction_frame {
                        // The following will error if the frame is not local
                        let dcm_vnc2inertial = this_xi.orbit.dcm_from_traj_frame(frame).unwrap();
                        let velocity_correction =
                            dcm_vnc2inertial * state_correction.fixed_rows::<3>(3);
                        this_xi.orbit.apply_dv(velocity_correction);
                    } else {
                        this_xi = xi + state_correction;
                    }
                }

                let this_xf = if finite_burn_target {
                    // Propagate normally until start of maneuver
                    let pre_mnvr = this_prop.with(cur_xi).until_epoch(this_mnvr.start).unwrap();
                    // Add this maneuver to the dynamics, make sure that we don't over-step this maneuver
                    let prop_opts = this_prop.opts;
                    this_prop.set_max_step(this_mnvr.duration());
                    this_prop.dynamics = this_prop.dynamics.with_guidance_law(Arc::new(this_mnvr));
                    let post_mnvr = this_prop
                        .with(pre_mnvr.with_guidance_mode(GuidanceMode::Thrust))
                        .until_epoch(this_mnvr.end)
                        .unwrap();
                    // Reset the propagator options to their previous configuration
                    this_prop.opts = prop_opts;
                    // And propagate until the achievement epoch
                    this_prop
                        .with(post_mnvr)
                        .until_epoch(self.achievement_epoch)
                        .unwrap()
                        .orbit
                } else {
                    this_prop
                        .with(this_xi)
                        .until_epoch(self.achievement_epoch)
                        .unwrap()
                        .orbit
                };

                let xf_dual_obj_frame = match &self.objective_frame {
                    Some((frame, cosm)) => {
                        let orbit_obj_frame = cosm.frame_chg(&this_xf, *frame);
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
                    xf_dual_obj_frame.partial_for(&obj.parameter).unwrap()
                };

                let this_achieved = partial.real();
                *jac_val = (this_achieved - achieved) / var.perturbation;
            });

            for (j, var, jac_val) in &pert_calc {
                println!(
                    "jac[({}, {})] = {} for {:?} and {:?}",
                    i, *j, jac_val, var, obj
                );
                self.jacobian[(i, *j)] = *jac_val;
            }
        }

        println!("resid: {}", self.residuals);
    }

    fn params(&self) -> SVector<f64, V> {
        self.control
    }

    fn residuals(&self) -> Option<SVector<f64, O>> {
        Some(self.residuals)
    }

    fn jacobian(&self) -> Option<SMatrix<f64, O, V>> {
        Some(self.jacobian)
        // Some(pseudo_inverse!(self.jacobian).unwrap())
    }
}

impl<'a, E: ErrorCtrl, const V: usize, const O: usize> Optimizer<'a, E, V, O>
where
    Const<V>: ToTypenum,
    Const<O>: ToTypenum,
    Const<O>: DimMin<Const<V>, Output = Const<V>> + DimMax<Const<V>, Output = Const<V>>,
{
    /// Differential correction using finite differencing
    #[allow(clippy::comparison_chain)]
    pub fn minimize(
        &self,
        initial_state: Spacecraft,
        correction_epoch: Epoch,
        achievement_epoch: Epoch,
    ) -> Result<(), NyxError> {
        // Check the variables and builds the initial guess
        let mut initial_control = SVector::<f64, V>::zeros();
        // Apply the initial guess
        for (i, var) in self.variables.iter().enumerate() {
            // Check the validity (this function will report to log and raise an error)
            var.valid()?;
            // Check that there is no attempt to target a position in a local frame
            if self.correction_frame.is_some() && var.component.vec_index() < 3 {
                // Then this is a position correction, which is not allowed if a frame is provided!
                let msg = format!(
                    "Variable is in frame {} but that frame cannot be used for a {:?} correction",
                    self.correction_frame.unwrap(),
                    var.component
                );
                error!("{}", msg);
                return Err(NyxError::Targeter(TargetingError::FrameError(msg)));
            }

            initial_control[i] = var.init_guess;
        }
        let mut instance = OptimizerInstance {
            prop: &self.prop.clone(),
            objectives: self.objectives,
            objective_frame: self.objective_frame.clone(),
            variables: self.variables,
            correction_frame: self.correction_frame,
            spacecraft: initial_state,
            achievement_epoch,
            correction_epoch,
            control: initial_control,
            // residuals: self.residuals.clone(),
            // TODO: Need a `step` function to compute the residuals without any correction.
            residuals: SVector::zeros(),
            jacobian: SMatrix::zeros(),
        };

        instance.set_params(&initial_control);
        println!("Init resid: {}", instance.residuals);
        let (result, report) = LevenbergMarquardt::new()
            .with_patience(10)
            .minimize(instance);

        println!("{:?}", report);

        println!(
            "Result correction: {}\t\t{} km/s",
            result.control,
            result.control.norm()
        );

        Ok(())
    }
}
