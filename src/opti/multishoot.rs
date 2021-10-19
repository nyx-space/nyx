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

        // Build each node successively
        let mut nodes = Vec::with_capacity(node_count);
        let mut epochs = Vec::with_capacity(node_count);
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
            // 1. Solve a simple differential corrector between each node
            // Note that we use the hyperdual formulation to solve the differential correction because
            // it's faster and more precise for small arcs.
            let mut initial_state = self.x0;
            let mut all_dvs = DVector::from_element(3 * self.nodes.len(), 0.0);
            let mut outer_jacobian =
                DMatrix::from_element(3 * self.nodes.len(), 3 * (self.nodes.len() - 2), 0.0);
            // Build the vector of the intermediate nodes, excluding the starting and final positions.
            let mut node_vector = DVector::from_element(3 * (self.nodes.len() - 2), 0.0);

            for i in 1..self.nodes.len() {
                // Add the current node info to the node vector
                node_vector[3 * (i - 1)] = self.nodes[i][0].desired_value;
                node_vector[3 * (i - 1) + 1] = self.nodes[i][1].desired_value;
                node_vector[3 * (i - 1) + 2] = self.nodes[i][2].desired_value;
                // Build the targeter
                let tgt = Targeter::delta_v(self.prop, self.nodes[i].to_vec());
                let sol = tgt.try_achieve_from_dual(
                    initial_state,
                    initial_state.epoch(),
                    self.epochs[i],
                )?;
                // Store the Δv and the initial state for the next targeter.
                initial_state = sol.achieved;
                let nominal_delta_v =
                    Vector3::new(sol.correction[0], sol.correction[1], sol.correction[2]);

                all_dvs[3 * (i - 1)] = nominal_delta_v[0];
                all_dvs[3 * (i - 1) + 1] = nominal_delta_v[1];
                all_dvs[3 * (i - 1) + 2] = nominal_delta_v[2];

                // Now, let's perturb the position of each node component to compute the
                // partial derivative of the delta_v with respect to each component of the node position.
                for j in 0..3 {
                    let mut next_node = self.nodes[i].to_vec();
                    next_node[j].desired_value += next_node[j].tolerance;
                    let inner_tgt = Targeter::delta_v(self.prop, self.nodes[i].to_vec());
                    let inner_sol = inner_tgt.try_achieve_from_dual(
                        initial_state,
                        initial_state.epoch(),
                        self.epochs[i],
                    )?;
                    // ∂Δv_x / ∂r_x
                    outer_jacobian[(i - 1 + j, 3 * (i - 1))] =
                        (inner_sol.correction[0] - nominal_delta_v[0]) / next_node[j].tolerance;
                    // ∂Δv_y / ∂r_x
                    outer_jacobian[(i - 1 + j, 3 * (i - 1) + 1)] =
                        (inner_sol.correction[1] - nominal_delta_v[1]) / next_node[j].tolerance;
                    // ∂Δv_z / ∂r_x
                    outer_jacobian[(i - 1 + j, 3 * (i - 1) + 2)] =
                        (inner_sol.correction[2] - nominal_delta_v[2]) / next_node[j].tolerance;
                }
            }
            // Compute the cost -- used to stop the algorithm if it does not change much.
            let cost_vec = &outer_jacobian * &node_vector;
            let new_cost = match cost {
                CostFunction::MinimumEnergy => cost_vec.dot(&cost_vec),
                CostFunction::MinimumFuel => cost_vec.dot(&cost_vec).sqrt(),
            };
            let cost_improvmt = (new_cost - prev_cost) / new_cost.abs();
            println!(
                "new_cost = {:.3}\timprovement = {:.3}",
                new_cost, cost_improvmt
            );
            if cost_improvmt < 0.01 {
                // If the cost does not improve by more than 1%, stop iteration
                return Ok(());
            }

            prev_cost = new_cost;

            // 2. Solve for the next position of the nodes using a pseudo inverse.
            let delta_r = pseudo_inverse(outer_jacobian, NyxError::SingularJacobian)? * cost_vec;

            // 3. Apply the correction to the node positions and iterator
            node_vector -= delta_r;
            for (i, val) in node_vector.iter().enumerate() {
                let node_no = i % 3;
                let component_no = i - node_no;
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
