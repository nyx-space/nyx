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

pub use super::CostFunction;
// use crate::dynamics::guidance::{FiniteBurns, Mnvr};
use crate::linalg::{DMatrix, DVector, SVector};
use crate::md::opti::solution::TargeterSolution;
use crate::md::optimizer::Optimizer;
use crate::md::ui::*;
use crate::propagators::error_ctrl::ErrorCtrl;
use crate::pseudo_inverse;
use crate::{Orbit, Spacecraft};

use std::fmt;

pub trait MultishootNode<const O: usize>: Copy + Into<[Objective; O]> {
    fn epoch(&self) -> Epoch;
    fn update_component(&mut self, component: usize, add_val: f64);
}

/// Multiple shooting is an optimization method.
/// Source of implementation: "Low Thrust Optimization in Cislunar and Translunar space", 2018 Nathan Re (Parrish)
/// OT: size of the objectives for each node (e.g. 3 if the objectives are X, Y, Z).
/// VT: size of the variables for targeter node (e.g. 4 if the objectives are thrust direction (x,y,z) and thrust level).
pub struct MultipleShooting<
    'a,
    E: ErrorCtrl,
    T: MultishootNode<OT>,
    const VT: usize,
    const OT: usize,
> {
    /// The propagator setup (kind, stages, etc.)
    pub prop: &'a Propagator<'a, SpacecraftDynamics<'a>, E>,
    /// List of nodes of the optimal trajectory
    pub targets: Vec<T>,
    /// Starting point, must be a spacecraft equipped with a thruster
    pub x0: Spacecraft,
    /// Destination (Should this be the final node?)
    pub xf: Orbit,
    pub current_iteration: usize,
    /// The maximum number of iterations allowed
    pub max_iterations: usize,
    /// Threshold after which the outer loop is considered to have converged,
    /// e.g. 0.01 means that a 1% of less improvement in case between two iterations
    /// will stop the iterations.
    pub improvement_threshold: f64,
    /// The kind of correction to apply to achieve the objectives
    pub variables: [Variable; VT],
    pub all_dvs: Vec<SVector<f64, VT>>,
}

impl<'a, E: ErrorCtrl, T: MultishootNode<OT>, const VT: usize, const OT: usize>
    MultipleShooting<'a, E, T, VT, OT>
{
    /// Solve the multiple shooting problem by finding the arrangement of nodes to minimize the cost function.
    pub fn solve(
        &mut self,
        cost: CostFunction,
    ) -> Result<MultipleShootingSolution<T, OT>, NyxError> {
        let mut prev_cost = 1e12; // We don't use infinity because we compare a ratio of cost
        for it in 0..self.max_iterations {
            let mut initial_states = Vec::with_capacity(self.targets.len());
            initial_states.push(self.x0);
            let mut outer_jacobian =
                DMatrix::from_element(3 * self.targets.len(), OT * (self.targets.len() - 1), 0.0);
            let mut cost_vec = DVector::from_element(3 * self.targets.len(), 0.0);

            // Reset the all_dvs
            self.all_dvs = Vec::with_capacity(self.all_dvs.len());

            for i in 0..self.targets.len() {
                /* ***
                 ** 1. Solve the delta-v differential corrector between each node
                 ** *** */
                let tgt = Optimizer {
                    prop: self.prop,
                    objectives: self.targets[i].into(),
                    variables: self.variables,
                    iterations: 100,
                    objective_frame: None,
                    correction_frame: None,
                };
                let sol = match tgt.try_achieve_dual(
                    initial_states[i],
                    initial_states[i].epoch(),
                    self.targets[i].epoch(),
                ) {
                    Ok(sol) => sol,
                    Err(e) => return Err(NyxError::MultipleShootingTargeter(i, Box::new(e))),
                };

                let nominal_delta_v = sol.correction;

                self.all_dvs.push(nominal_delta_v);
                // Store the Δv and the initial state for the next targeter.
                initial_states.push(sol.achieved_state);
            }
            // NOTE: We have two separate loops because we need the initial state of node i+2 for the dv computation
            // of the third entry to the outer jacobian.
            for i in 0..(self.targets.len() - 1) {
                /* ***
                 ** 2. Perturb each node and compute the partial of the Δv for the (i-1), i, and (i+1) nodes
                 ** where the partial on the i+1 -th node is just the difference between the velocity at the
                 ** achieved state and the initial state at that node.
                 ** We don't perturb the endpoint node
                 ** *** */

                for axis in 0..OT {
                    /* ***
                     ** 2.A. Perturb the i-th node
                     ** *** */
                    let mut next_node = self.targets[i].into();
                    next_node[axis].desired_value += next_node[axis].tolerance;
                    /* ***
                     ** 2.b. Rerun the targeter from the previous node to this one
                     ** Note that because the first initial_state is x0, the i-th "initial state"
                     ** is the initial state to reach the i-th node.
                     ** *** */
                    let inner_tgt_a = Optimizer::delta_v(self.prop, next_node);
                    let inner_sol_a = match inner_tgt_a.try_achieve_dual(
                        initial_states[i],
                        initial_states[i].epoch(),
                        self.targets[i].epoch(),
                    ) {
                        Ok(sol) => sol,
                        Err(e) => return Err(NyxError::MultipleShootingTargeter(i, Box::new(e))),
                    };

                    // ∂Δv_x / ∂r_x
                    outer_jacobian[(3 * i, OT * i + axis)] = (inner_sol_a.correction[0]
                        - self.all_dvs[i][0])
                        / next_node[axis].tolerance;
                    // ∂Δv_y / ∂r_x
                    outer_jacobian[(3 * i + 1, OT * i + axis)] = (inner_sol_a.correction[1]
                        - self.all_dvs[i][1])
                        / next_node[axis].tolerance;
                    // ∂Δv_z / ∂r_x
                    outer_jacobian[(3 * i + 2, OT * i + axis)] = (inner_sol_a.correction[2]
                        - self.all_dvs[i][2])
                        / next_node[axis].tolerance;

                    /* ***
                     ** 2.C. Rerun the targeter from the new state at the perturbed node to the next unpertubed node
                     ** *** */
                    let inner_tgt_b = Optimizer::delta_v(self.prop, self.targets[i + 1].into());
                    let inner_sol_b = inner_tgt_b.try_achieve_dual(
                        inner_sol_a.achieved_state,
                        inner_sol_a.achieved_state.epoch(),
                        self.targets[i + 1].epoch(),
                    )?;

                    // Compute the partials wrt the next Δv
                    // ∂Δv_x / ∂r_x
                    outer_jacobian[(3 * (i + 1), OT * i + axis)] = (inner_sol_b.correction[0]
                        - self.all_dvs[i + 1][0])
                        / next_node[axis].tolerance;
                    // ∂Δv_y / ∂r_x
                    outer_jacobian[(3 * (i + 1) + 1, OT * i + axis)] = (inner_sol_b.correction[1]
                        - self.all_dvs[i + 1][1])
                        / next_node[axis].tolerance;
                    // ∂Δv_z / ∂r_x
                    outer_jacobian[(3 * (i + 1) + 2, OT * i + axis)] = (inner_sol_b.correction[2]
                        - self.all_dvs[i + 1][2])
                        / next_node[axis].tolerance;

                    /* ***
                     ** 2.D. Compute the difference between the arrival and departure velocities and node i+1
                     ** *** */
                    if i < self.targets.len() - 3 {
                        let dv_ip1 = inner_sol_b.achieved_state.orbit.velocity()
                            - initial_states[i + 2].orbit.velocity();
                        // ∂Δv_x / ∂r_x
                        outer_jacobian[(3 * (i + 2), OT * i + axis)] =
                            dv_ip1[0] / next_node[axis].tolerance;
                        // ∂Δv_y / ∂r_x
                        outer_jacobian[(3 * (i + 2) + 1, OT * i + axis)] =
                            dv_ip1[1] / next_node[axis].tolerance;
                        // ∂Δv_z / ∂r_x
                        outer_jacobian[(3 * (i + 2) + 2, OT * i + axis)] =
                            dv_ip1[2] / next_node[axis].tolerance;
                    }
                }
            }

            // Build the cost vector
            for i in 0..self.targets.len() {
                for j in 0..3 {
                    cost_vec[3 * i + j] = self.all_dvs[i][j];
                }
            }

            // Compute the cost -- used to stop the algorithm if it does not change much.
            let new_cost = match cost {
                CostFunction::MinimumEnergy => cost_vec.dot(&cost_vec),
                CostFunction::MinimumFuel => cost_vec.dot(&cost_vec).sqrt(),
            };

            // If the new cost is greater than the previous one, then the cost improvement is negative.
            let cost_improvmt = (prev_cost - new_cost) / new_cost.abs();
            // If the cost does not improve by more than threshold stop iteration
            match cost {
                CostFunction::MinimumEnergy => info!(
                    "Multiple shooting iteration #{}\t\tCost = {:.3} km^2/s^2\timprovement = {:.2}%",
                    it,
                    new_cost,
                    100.0 * cost_improvmt
                ),
                CostFunction::MinimumFuel => info!(
                    "Multiple shooting iteration #{}\t\tCost = {:.3} km/s\timprovement = {:.2}%",
                    it,
                    new_cost,
                    100.0 * cost_improvmt
                ),
            };
            if cost_improvmt.abs() < self.improvement_threshold {
                info!("Improvement below desired threshold. Running targeter on computed nodes.");

                /* ***
                 ** FIN -- Check the impulsive burns work and return all targeter solutions
                 ** *** */
                let mut ms_sol = MultipleShootingSolution {
                    x0: self.x0,
                    xf: self.xf,
                    nodes: self.targets.clone(),
                    solutions: Vec::with_capacity(self.targets.len()),
                };
                let mut initial_states = Vec::with_capacity(self.targets.len());
                initial_states.push(self.x0);

                for (i, node) in self.targets.iter().enumerate() {
                    // Run the unpertubed targeter
                    let tgt = Optimizer::delta_v(self.prop, (*node).into());
                    let sol = tgt.try_achieve_dual(
                        initial_states[i],
                        initial_states[i].epoch(),
                        node.epoch(),
                    )?;
                    initial_states.push(sol.achieved_state);
                    ms_sol.solutions.push(sol);
                }

                return Ok(ms_sol);
            }

            prev_cost = new_cost;
            // 2. Solve for the next position of the nodes using a pseudo inverse.
            let inv_jac = match pseudo_inverse!(&outer_jacobian) {
                Ok(inv_jac) => inv_jac,
                Err(e) => {
                    error!("Singular Jacobian {:.3}", outer_jacobian);
                    return Err(e);
                }
            };
            let delta_r = inv_jac * cost_vec;
            // 3. Apply the correction to the node positions and iterator
            let node_vector = -delta_r;
            for (i, val) in node_vector.iter().enumerate() {
                let node_no = i / 3;
                let component_no = i % OT;
                self.targets[node_no].update_component(component_no, *val);
            }
            self.current_iteration += 1;
        }
        Err(NyxError::MaxIterReached(format!("{}", self.max_iterations)))
    }
}

impl<'a, E: ErrorCtrl, T: MultishootNode<OT>, const VT: usize, const OT: usize> fmt::Display
    for MultipleShooting<'a, E, T, VT, OT>
{
    #[allow(clippy::or_fun_call, clippy::clone_on_copy)]
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut nodemsg = String::from("");
        // Add the starting point too
        nodemsg.push_str(&format!(
            "[{:.3}, {:.3}, {:.3}, {}, {}, {}, {}, {}, {}],\n",
            self.x0.orbit.x,
            self.x0.orbit.y,
            self.x0.orbit.z,
            self.current_iteration,
            0.0,
            0.0,
            0.0,
            0.0,
            0
        ));

        for (i, node) in self.targets.iter().enumerate() {
            let objectives: [Objective; OT] = (*node).into();
            let mut this_nodemsg = String::from("");
            for obj in &objectives {
                this_nodemsg.push_str(&format!("{:.3}, ", obj.desired_value));
            }
            let mut this_costmsg = String::from("");
            let dv = match self.all_dvs.get(i) {
                Some(dv) => dv.clone(),
                None => SVector::<f64, VT>::zeros(),
            };
            for val in &dv {
                this_costmsg.push_str(&format!("{}, ", val));
            }
            if VT == 3 {
                // Add the norm of the control
                this_costmsg.push_str(&format!("{}, ", dv.norm()));
            }
            nodemsg.push_str(&format!(
                "[{}{}, {}{}],\n",
                this_nodemsg,
                self.current_iteration,
                this_nodemsg,
                i + 1
            ));
        }
        write!(f, "{}", nodemsg)
    }
}

#[derive(Clone, Debug)]
pub struct MultipleShootingSolution<T: MultishootNode<O>, const O: usize> {
    pub x0: Spacecraft,
    pub xf: Orbit,
    pub nodes: Vec<T>,
    pub solutions: Vec<TargeterSolution<3, O>>,
}

impl<T: MultishootNode<O>, const O: usize> fmt::Display for MultipleShootingSolution<T, O> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for sol in &self.solutions {
            write!(f, "{}", sol)?;
        }
        Ok(())
    }
}

impl<T: MultishootNode<O>, const O: usize> MultipleShootingSolution<T, O> {
    /// Allows building the trajectories between different nodes
    /// This will rebuild the targeters and apply the solutions sequentially
    pub fn build_trajectories<'a, E: ErrorCtrl>(
        &self,
        prop: &'a Propagator<'a, SpacecraftDynamics, E>,
    ) -> Result<Vec<ScTraj>, NyxError> {
        let mut trajz = Vec::with_capacity(self.nodes.len());

        for (i, node) in self.nodes.iter().enumerate() {
            let (_, traj) =
                Optimizer::delta_v(prop, (*node).into()).apply_with_traj(&self.solutions[i])?;
            trajz.push(traj);
        }

        Ok(trajz)
    }
}
