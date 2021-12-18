/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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
use hifitime::TimeUnitHelper;
use levenberg_marquardt::{LeastSquaresProblem, LevenbergMarquardt};

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
    pub jacobian: SMatrix<f64, M, N>,
}

impl<'a, E: ErrorCtrl, const V: usize, const O: usize> OptimizerInstance<'a, E, V, O>
where
    Const<V>: ToTypenum,
    Const<O>: ToTypenum,
    Const<O>: DimMin<Const<V>, Output = Const<V>> + DimMax<Const<V>, Output = Const<V>>,
{
    /// Compute the residuals for the given control
    pub fn residuals_for_ctrl(
        &self,
        attempted_control: &SVector<f64, V>,
    ) -> Result<SVector<f64, O>, NyxError> {
        // Start by applying the control to the initial spacecraft

        // Create a maneuver in case we need it (finite burn targeting)
        let mut mnvr = Mnvr {
            start: self.correction_epoch,
            end: self.achievement_epoch,
            thrust_lvl: 1.0,
            alpha_inplane_radians: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
            delta_outofplane_radians: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
            frame: Frame::RCN,
        };

        // Create a full state correction in case we need that too
        let mut state_correction = Vector6::<f64>::zeros();
        let mut finite_burn_target = false;

        let mut delta = *attempted_control;
        println!("ctrl = {}", attempted_control);
        for (i, var) in self.variables.iter().enumerate() {
            info!(
                "Correction {:?}{} (element {}): {}",
                var.component,
                match self.correction_frame {
                    Some(f) => format!(" in {:?}", f),
                    None => format!(""),
                },
                i,
                delta[i]
            );

            let corr = delta[i];

            if var.component.is_finite_burn() {
                finite_burn_target = true;
                // Modify the maneuver, but do not change the epochs of the maneuver unless the change is greater than one millisecond
                match var.component {
                    Vary::Duration => {
                        if corr.abs() > 1e-3 {
                            // Check that we are within the bounds
                            let init_duration_s =
                                (self.correction_epoch - self.achievement_epoch).in_seconds();
                            let acceptable_corr = var.apply_bounds(init_duration_s).seconds();
                            mnvr.end = mnvr.start + acceptable_corr;
                        }
                    }
                    Vary::EndEpoch => {
                        if corr.abs() > 1e-3 {
                            // Check that we are within the bounds
                            let total_end_corr =
                                (mnvr.end + corr.seconds() - self.achievement_epoch).in_seconds();
                            let acceptable_corr = var.apply_bounds(total_end_corr).seconds();
                            mnvr.end = mnvr.end + acceptable_corr;
                        }
                    }
                    Vary::StartEpoch => {
                        if corr.abs() > 1e-3 {
                            // Check that we are within the bounds
                            let total_start_corr =
                                (mnvr.start + corr.seconds() - self.correction_epoch).in_seconds();
                            let acceptable_corr = var.apply_bounds(total_start_corr).seconds();
                            mnvr.end = mnvr.end + acceptable_corr;

                            mnvr.start = mnvr.start + corr.seconds()
                        }
                    }
                    Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => {
                        mnvr.alpha_inplane_radians = mnvr
                            .alpha_inplane_radians
                            .add_val_in_order(corr, var.component.vec_index())
                            .unwrap();
                    }
                    Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                        mnvr.delta_outofplane_radians = mnvr
                            .delta_outofplane_radians
                            .add_val_in_order(corr, var.component.vec_index())
                            .unwrap();
                    }
                    Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => {
                        let mut vector = mnvr.direction();
                        vector[var.component.vec_index()] += corr;
                        var.ensure_bounds(&mut vector[var.component.vec_index()]);
                        mnvr.set_direction(vector);
                    }
                    Vary::ThrustRateX | Vary::ThrustRateY | Vary::ThrustRateZ => {
                        let mut vector = mnvr.rate();
                        vector[(var.component.vec_index() - 1) % 3] += corr;
                        mnvr.set_rate(vector);
                    }
                    Vary::ThrustAccelX | Vary::ThrustAccelY | Vary::ThrustAccelZ => {
                        let mut vector = mnvr.accel();
                        vector[(var.component.vec_index() - 1) % 3] += corr;
                        mnvr.set_accel(vector);
                    }
                    Vary::ThrustLevel => {
                        mnvr.thrust_lvl -= corr;
                        var.ensure_bounds(&mut mnvr.thrust_lvl);
                    }
                    _ => unreachable!(),
                }
            } else {
                state_correction[var.component.vec_index()] += delta[i];
            }
        }

        let mut xi = self.spacecraft;

        if !finite_burn_target {
            xi.orbit = xi.orbit + state_correction;
            println!("Corrected state: {}", xi);
        }

        // Now that we have updated the state or the finite burn, let's propagate until the achievement epoch.
        let xf = if finite_burn_target {
            info!("{}", mnvr);
            let mut prop = self.prop.clone();
            let prop_opts = prop.opts;
            let pre_mnvr = prop.with(xi).until_epoch(mnvr.start)?;
            prop.dynamics = prop.dynamics.with_ctrl(Arc::new(mnvr));
            prop.set_max_step(mnvr.duration());
            let post_mnvr = prop
                .with(pre_mnvr.with_guidance_mode(GuidanceMode::Thrust))
                .until_epoch(mnvr.end)?;
            // Reset the propagator options to their previous configuration
            prop.opts = prop_opts;
            // And propagate until the achievement epoch
            prop.with(post_mnvr)
                .until_epoch(self.achievement_epoch)?
                .orbit
        } else {
            self.prop
                .with(xi)
                .until_epoch(self.achievement_epoch)?
                .orbit
        };

        let xf_dual_obj_frame = match &self.objective_frame {
            Some((frame, cosm)) => {
                let orbit_obj_frame = cosm.frame_chg(&xf, *frame);
                OrbitDual::from(orbit_obj_frame)
            }
            None => OrbitDual::from(xf),
        };

        // Build the logging information
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

        // Build the error vector
        let mut residuals = SVector::<f64, O>::zeros();
        let mut converged = true;

        // Build debugging information
        let mut objmsg = Vec::with_capacity(self.objectives.len());

        for (i, obj) in self.objectives.iter().enumerate() {
            let partial = if obj.parameter.is_b_plane() {
                let b_plane = BPlane::from_dual(xf_dual_obj_frame)?;
                match obj.parameter {
                    StateParameter::BdotR => b_plane.b_r,
                    StateParameter::BdotT => b_plane.b_t,
                    StateParameter::BLTOF => b_plane.ltof_s,
                    _ => unreachable!(),
                }
            } else {
                xf_dual_obj_frame.partial_for(&obj.parameter)?
            };

            let achieved = partial.real();

            let (ok, param_err) = obj.assess_raw(achieved);
            if !ok {
                converged = false;
            }
            residuals[i] = param_err;

            objmsg.push(format!(
                "\t{:?}: achieved = {:>width$.prec$}\t desired = {:>width$.prec$}\t scaled error = {:>width$.prec$}",
                obj.parameter,
                achieved,
                obj.desired_value,
                param_err, width=width, prec=max_obj_tol
            ));
        }

        dbg!(converged);
        Ok(residuals)
    }
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
        // Start by applying the control to the initial spacecraft

        // Create a maneuver in case we need it (finite burn targeting)
        let mut mnvr = Mnvr {
            start: self.correction_epoch,
            end: self.achievement_epoch,
            thrust_lvl: 1.0,
            alpha_inplane_radians: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
            delta_outofplane_radians: CommonPolynomial::Quadratic(0.0, 0.0, 0.0),
            frame: Frame::RCN,
        };

        // Create a full state correction in case we need that too
        let mut state_correction = Vector6::<f64>::zeros();
        let mut finite_burn_target = false;

        let mut delta = *attempted_control;
        println!("ctrl = {}", attempted_control);
        for (i, var) in self.variables.iter().enumerate() {
            info!(
                "Correction {:?}{} (element {}): {}",
                var.component,
                match self.correction_frame {
                    Some(f) => format!(" in {:?}", f),
                    None => format!(""),
                },
                i,
                delta[i]
            );

            let corr = delta[i];

            if var.component.is_finite_burn() {
                finite_burn_target = true;
                // Modify the maneuver, but do not change the epochs of the maneuver unless the change is greater than one millisecond
                match var.component {
                    Vary::Duration => {
                        if corr.abs() > 1e-3 {
                            // Check that we are within the bounds
                            let init_duration_s =
                                (self.correction_epoch - self.achievement_epoch).in_seconds();
                            let acceptable_corr = var.apply_bounds(init_duration_s).seconds();
                            mnvr.end = mnvr.start + acceptable_corr;
                        }
                    }
                    Vary::EndEpoch => {
                        if corr.abs() > 1e-3 {
                            // Check that we are within the bounds
                            let total_end_corr =
                                (mnvr.end + corr.seconds() - self.achievement_epoch).in_seconds();
                            let acceptable_corr = var.apply_bounds(total_end_corr).seconds();
                            mnvr.end = mnvr.end + acceptable_corr;
                        }
                    }
                    Vary::StartEpoch => {
                        if corr.abs() > 1e-3 {
                            // Check that we are within the bounds
                            let total_start_corr =
                                (mnvr.start + corr.seconds() - self.correction_epoch).in_seconds();
                            let acceptable_corr = var.apply_bounds(total_start_corr).seconds();
                            mnvr.end = mnvr.end + acceptable_corr;

                            mnvr.start = mnvr.start + corr.seconds()
                        }
                    }
                    Vary::MnvrAlpha | Vary::MnvrAlphaDot | Vary::MnvrAlphaDDot => {
                        mnvr.alpha_inplane_radians = mnvr
                            .alpha_inplane_radians
                            .add_val_in_order(corr, var.component.vec_index())
                            .unwrap();
                    }
                    Vary::MnvrDelta | Vary::MnvrDeltaDot | Vary::MnvrDeltaDDot => {
                        mnvr.delta_outofplane_radians = mnvr
                            .delta_outofplane_radians
                            .add_val_in_order(corr, var.component.vec_index())
                            .unwrap();
                    }
                    Vary::ThrustX | Vary::ThrustY | Vary::ThrustZ => {
                        let mut vector = mnvr.direction();
                        vector[var.component.vec_index()] += corr;
                        var.ensure_bounds(&mut vector[var.component.vec_index()]);
                        mnvr.set_direction(vector);
                    }
                    Vary::ThrustRateX | Vary::ThrustRateY | Vary::ThrustRateZ => {
                        let mut vector = mnvr.rate();
                        vector[(var.component.vec_index() - 1) % 3] += corr;
                        mnvr.set_rate(vector);
                    }
                    Vary::ThrustAccelX | Vary::ThrustAccelY | Vary::ThrustAccelZ => {
                        let mut vector = mnvr.accel();
                        vector[(var.component.vec_index() - 1) % 3] += corr;
                        mnvr.set_accel(vector);
                    }
                    Vary::ThrustLevel => {
                        mnvr.thrust_lvl -= corr;
                        var.ensure_bounds(&mut mnvr.thrust_lvl);
                    }
                    _ => unreachable!(),
                }
            } else {
                state_correction[var.component.vec_index()] += delta[i];
            }
        }

        let mut xi = self.spacecraft;

        if !finite_burn_target {
            xi.orbit = xi.orbit + state_correction;
            println!("Corrected state: {}", xi);
        }

        // Now that we have updated the state or the finite burn, let's propagate until the achievement epoch.
        let xf = if finite_burn_target {
            info!("{}", mnvr);
            let mut prop = self.prop.clone();
            let prop_opts = prop.opts;
            let pre_mnvr = prop.with(xi).until_epoch(mnvr.start).unwrap();
            prop.dynamics = prop.dynamics.with_ctrl(Arc::new(mnvr));
            prop.set_max_step(mnvr.duration());
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
                .with(xi)
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

        // Build the logging information
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

        // Build the error vector
        let mut residuals = SVector::<f64, O>::zeros();
        let mut converged = true;

        // Build debugging information
        let mut objmsg = Vec::with_capacity(self.objectives.len());

        // The Jacobian includes the sensitivity of each objective with respect to each variable for the whole trajectory.
        // As such, it includes the STM of that variable for the whole propagation arc.
        // let mut jac = DMatrix::from_element(self.objectives.len(), self.variables.len(), 0.0);

        for (i, obj) in self.objectives.iter().enumerate() {
            let partial = if obj.parameter.is_b_plane() {
                let b_plane = BPlane::from_dual(xf_dual_obj_frame).unwrap();
                match obj.parameter {
                    StateParameter::BdotR => b_plane.b_r,
                    StateParameter::BdotT => b_plane.b_t,
                    StateParameter::BLTOF => b_plane.ltof_s,
                    _ => unreachable!(),
                }
            } else {
                xf_dual_obj_frame.partial_for(&obj.parameter).unwrap()
            };

            let achieved = partial.real();

            let (ok, param_err) = obj.assess_raw(achieved);
            if !ok {
                converged = false;
            }
            self.residuals[i] = param_err;

            objmsg.push(format!(
                "\t{:?}: achieved = {:>width$.prec$}\t desired = {:>width$.prec$}\t scaled error = {:>width$.prec$}",
                obj.parameter,
                achieved,
                obj.desired_value,
                param_err, width=width, prec=max_obj_tol
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
                            this_mnvr.set_direction(vector);
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
                    let pre_mnvr = this_prop
                        .with(this_xi)
                        .until_epoch(this_mnvr.start)
                        .unwrap();
                    // Add this maneuver to the dynamics, make sure that we don't over-step this maneuver
                    let prop_opts = this_prop.opts;
                    this_prop.set_max_step(this_mnvr.duration());
                    this_prop.dynamics = this_prop.dynamics.with_ctrl(Arc::new(this_mnvr));
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

                let partial = if obj.parameter.is_b_plane() {
                    let b_plane = BPlane::from_dual(xf_dual_obj_frame).unwrap();
                    match obj.parameter {
                        StateParameter::BdotR => b_plane.b_r,
                        StateParameter::BdotT => b_plane.b_t,
                        StateParameter::BLTOF => b_plane.ltof_s,
                        _ => unreachable!(),
                    }
                } else {
                    xf_dual_obj_frame.partial_for(&obj.parameter).unwrap()
                };

                let this_achieved = partial.real();
                *jac_val = (this_achieved - achieved) / var.perturbation;
            });

            for (j, _, jac_val) in &pert_calc {
                self.jacobian[(i, *j)] = *jac_val;
            }
        }

        // println!("resid: {}", self.residuals);
        for obj in &objmsg {
            info!("{}", obj);
        }
        dbg!(converged);
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
        initial_control: SVector<f64, V>,
    ) -> Result<(), NyxError> {
        let mut instance = OptimizerInstance {
            prop: &self.prop.clone(),
            objectives: self.objectives.clone(),
            objective_frame: self.objective_frame.clone(),
            variables: self.variables.clone(),
            correction_frame: self.correction_frame.clone(),
            spacecraft: initial_state,
            achievement_epoch: achievement_epoch,
            correction_epoch: correction_epoch,
            control: initial_control,
            // residuals: self.residuals.clone(),
            // TODO: Need a `step` function to compute the residuals without any correction.
            residuals: SVector::zeros(),
            jacobian: SMatrix::zeros(),
        };

        instance.set_params(&initial_control);
        println!("Init ctrl : {}", instance.control);
        println!("Init resid: {}", instance.residuals);

        let (result, report) = LevenbergMarquardt::new()
            .with_patience(50)
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
