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

pub use super::CostFunction;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DMatrix, DVector, DefaultAllocator, Vector3};
use crate::md::targeter::{Objective, Targeter};
use crate::md::ui::*;
use crate::propagators::error_ctrl::ErrorCtrl;
use crate::time::Epoch;
use crate::utils::pseudo_inverse;
use crate::{Orbit, Spacecraft};

use std::fmt;

/// Multiple shooting is an optimization method.
/// Source of implementation: "Low Thrust Optimization in Cislunar and Translunar space", 2018 Nathan Re (Parrish)
pub struct MultipleShooting<'a, D: Dynamics<StateType = Spacecraft>, E: ErrorCtrl>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
{
    /// The propagator setup (kind, stages, etc.)
    pub prop: &'a Propagator<'a, D, E>,
    /// List of position nodes of the optimal trajectory
    pub nodes: Vec<[Objective; 3]>,
    /// List of epoch at which the nodes must be met
    pub epochs: Vec<Epoch>,
    /// Starting point, must be a spacecraft equipped with a thruster
    pub x0: Spacecraft,
    /// Destination (Should this be the final node?)
    pub xf: Orbit,
    /// The maximum number of iterations allowed
    pub max_iterations: usize,
}

impl<'a, D: Dynamics<StateType = Spacecraft>, E: ErrorCtrl> MultipleShooting<'a, D, E>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>
        + Allocator<f64, <D::StateType as State>::Size>,
{
    /// Builds a multiple shooting structure assuming that the optimal trajectory is a straight line
    /// between the start and end points. The position of the nodes will be update at each iteration
    /// of the outer loop.
    pub fn equidistant_nodes(
        x0: Spacecraft,
        xf: Orbit,
        node_count: usize,
        prop: &'a Propagator<'a, D, E>,
    ) -> Result<Self, NyxError> {
        if node_count < 3 {
            error!("At least three nodes are needed for a multiple shooting optimization");
            return Err(NyxError::UnderdeterminedProblem);
        }

        // Compute the direction of the objective
        let mut direction = xf.radius() - x0.orbit.radius();
        let distance_increment = direction.norm() / (node_count as f64);
        let duration_increment = (xf.epoch() - x0.epoch()) / (node_count as f64);
        direction /= direction.norm();

        // Build each node successively (includes xf)
        let mut nodes = Vec::with_capacity(node_count + 1);
        let mut epochs = Vec::with_capacity(node_count + 1);
        let mut prev_node_radius = x0.orbit.radius();
        let mut prev_node_epoch = x0.epoch();
        for _ in 0..node_count {
            // Compute the position we want.
            let this_node = prev_node_radius + distance_increment * direction;
            nodes.push([
                Objective::new(StateParameter::X, this_node[0]),
                Objective::new(StateParameter::Y, this_node[1]),
                Objective::new(StateParameter::Z, this_node[2]),
            ]);
            let this_epoch = prev_node_epoch + duration_increment;
            epochs.push(this_epoch);
            prev_node_radius = this_node;
            prev_node_epoch = this_epoch;
        }
        Ok(Self {
            prop,
            nodes,
            epochs,
            x0,
            xf,
            max_iterations: 20,
        })
    }

    /// Solve the multiple shooting problem by finding the arrangement of nodes to minimize the cost function.
    pub fn solve(&mut self, cost: CostFunction) -> Result<(), NyxError> {
        let mut prev_cost = 1e12; // We don't use infinity because we compare a ratio of cost
        for _ in 0..self.max_iterations {
            println!("{}", self);
            /* ***
             ** 1. Solve a simple differential corrector between each node
             ** *** */
            let mut initial_states = Vec::with_capacity(self.nodes.len());
            initial_states.push(self.x0);
            let mut all_dvs = Vec::with_capacity(self.nodes.len());
            let mut outer_jacobian =
                DMatrix::from_element(3 * self.nodes.len(), 3 * self.nodes.len(), 0.0);
            // Build the vector of the intermediate nodes, excluding the starting and final positions.
            let mut node_vector = DVector::from_element(3 * self.nodes.len(), 0.0);

            for i in 0..self.nodes.len() {
                // Run the unpertubed targeter
                let tgt = Targeter::delta_v(self.prop, self.nodes[i].to_vec());
                let sol = tgt.try_achieve_fd(
                    initial_states[i],
                    initial_states[i].epoch(),
                    self.epochs[i],
                )?;

                let nominal_delta_v =
                    Vector3::new(sol.correction[0], sol.correction[1], sol.correction[2]);
                println!("\n{} => {}", i, nominal_delta_v);

                all_dvs.push(nominal_delta_v);
                // Store the Δv and the initial state for the next targeter.
                initial_states.push(sol.achieved);
            }

            /* ***
             ** 2. Perturb each node and compute the partial of the Δv for the (i-1), i, and (i+1) nodes
             ** where the partial on the i+1 -th node is just the difference between the velocity at the
             ** achieved state and the initial state at that node.
             ** *** */
            for i in 1..(self.nodes.len() - 1) {
                // Add the current node info to the node vector
                node_vector[3 * (i - 1)] = self.nodes[i][0].desired_value;
                node_vector[3 * (i - 1) + 1] = self.nodes[i][1].desired_value;
                node_vector[3 * (i - 1) + 2] = self.nodes[i][2].desired_value;

                for axis in 0..3 {
                    /* ***
                     ** 2.A. Perturb the i-th node
                     ** *** */
                    let mut next_node = self.nodes[i].to_vec();
                    next_node[axis].desired_value += next_node[axis].tolerance;
                    /* ***
                     ** 2.b. Rerun the targeter from the previous node to this one
                     ** Note that because the first initial_state is x0, the i-th "initial state"
                     ** is the initial state to reach the i-th node.
                     ** *** */

                    let inner_tgt_a = Targeter::delta_v(self.prop, next_node.to_vec());
                    let inner_sol_a = inner_tgt_a.try_achieve_fd(
                        initial_states[i],
                        initial_states[i].epoch(),
                        self.epochs[i],
                    )?;

                    let idx = i - 1;
                    // ∂Δv_x / ∂r_x
                    outer_jacobian[(3 * idx + axis, 3 * idx)] =
                        (inner_sol_a.correction[0] - all_dvs[i][0]) / next_node[axis].tolerance;
                    // ∂Δv_y / ∂r_x
                    outer_jacobian[(3 * idx + axis, 3 * idx + 1)] =
                        (inner_sol_a.correction[1] - all_dvs[i][1]) / next_node[axis].tolerance;
                    // ∂Δv_z / ∂r_x
                    outer_jacobian[(3 * idx + axis, 3 * idx + 2)] =
                        (inner_sol_a.correction[2] - all_dvs[i][2]) / next_node[axis].tolerance;

                    /* ***
                     ** 2.C. Rerun the targeter from the new state at the perturbed node to the next node
                     ** *** */
                    let inner_tgt_b = Targeter::delta_v(self.prop, self.nodes[i + 1].to_vec());
                    let inner_sol_b = inner_tgt_b.try_achieve_fd(
                        inner_sol_a.achieved,
                        inner_sol_a.achieved.epoch(),
                        self.epochs[i + 1],
                    )?;

                    // Compute the partials wrt the next Δv
                    // ∂Δv_x / ∂r_x
                    outer_jacobian[(3 * idx + axis, 3 * (idx + 1))] =
                        (inner_sol_a.correction[0] - all_dvs[i][0]) / next_node[axis].tolerance;
                    // ∂Δv_y / ∂r_x
                    outer_jacobian[(3 * idx + axis, 3 * (idx + 1) + 1)] =
                        (inner_sol_a.correction[1] - all_dvs[i][1]) / next_node[axis].tolerance;
                    // ∂Δv_z / ∂r_x
                    outer_jacobian[(3 * idx + axis, 3 * (idx + 1) + 2)] =
                        (inner_sol_a.correction[2] - all_dvs[i][2]) / next_node[axis].tolerance;

                    /* ***
                     ** 2.D. Compute the difference between the arrival and departure velocities and node i+1
                     ** *** */
                    let dv_ip1 = inner_sol_b.achieved.orbit.velocity()
                        - initial_states[i + 1].orbit.velocity();

                    // ∂Δv_x / ∂r_x
                    outer_jacobian[(3 * idx + axis, 3 * (idx + 2))] =
                        dv_ip1[0] / next_node[axis].tolerance;
                    // ∂Δv_y / ∂r_x
                    outer_jacobian[(3 * idx + axis, 3 * (idx + 2) + 1)] =
                        dv_ip1[1] / next_node[axis].tolerance;
                    // ∂Δv_z / ∂r_x
                    outer_jacobian[(3 * idx + axis, 3 * (idx + 2) + 2)] =
                        dv_ip1[2] / next_node[axis].tolerance;
                }
            }

            println!("{}", outer_jacobian);
            println!("node_vector = {}", node_vector);
            // Compute the cost -- used to stop the algorithm if it does not change much.
            let mut cost_vec = DVector::from_element(3 * self.nodes.len(), 0.0);
            for (i, dv) in all_dvs.iter().enumerate() {
                for j in 0..3 {
                    cost_vec[3 * i + j] = dv[j];
                }
            }
            // println!("cost_vec = {}", cost_vec);
            let new_cost = match cost {
                CostFunction::MinimumEnergy => cost_vec.dot(&cost_vec),
                CostFunction::MinimumFuel => cost_vec.dot(&cost_vec).sqrt(),
            };
            let cost_improvmt = (new_cost - prev_cost) / new_cost.abs();
            println!(
                "new_cost = {:.3}\timprovement = {:.3}",
                new_cost, cost_improvmt
            );
            if cost_improvmt.abs() < 0.01 {
                println!("We're done!");
                // If the cost does not improve by more than 1%, stop iteration
                return Ok(());
            }

            prev_cost = new_cost;

            // 2. Solve for the next position of the nodes using a pseudo inverse.
            let inv_jac = pseudo_inverse(outer_jacobian, NyxError::SingularJacobian)?;
            let delta_r = inv_jac * cost_vec;
            println!("delta_r = {}", delta_r);

            // 3. Apply the correction to the node positions and iterator
            node_vector = -delta_r;
            for (i, val) in node_vector.iter().enumerate() {
                let node_no = i / 3;
                let component_no = i % 3;
                // dbg!(node_no, component_no, i);
                self.nodes[node_no][component_no].desired_value += val;
            }
        }
        Err(NyxError::MaxIterReached(format!("{}", self.max_iterations)))
    }
}

impl<'a, D: Dynamics<StateType = Spacecraft>, E: ErrorCtrl> fmt::Display
    for MultipleShooting<'a, D, E>
where
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>
        + Allocator<f64, <D::StateType as State>::Size>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let mut nodemsg = String::from("");
        for (i, node) in self.nodes.iter().enumerate() {
            nodemsg.push_str(&format!(
                "\n\t\t@{} = [{:.3}, {:.3}, {:.3}]",
                self.epochs[i], node[0].desired_value, node[1].desired_value, node[2].desired_value
            ));
        }
        write!(f, "{}", nodemsg)
    }
}
