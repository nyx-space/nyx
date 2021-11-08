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
// use crate::dynamics::guidance::{FiniteBurns, Mnvr};
use super::ctrlnodes::Node;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DMatrix, DVector, DefaultAllocator, Vector3};
use crate::md::targeter::Targeter;
use crate::md::ui::*;
use crate::propagators::error_ctrl::ErrorCtrl;
use crate::utils::pseudo_inverse;
use crate::{Orbit, Spacecraft};
use std::sync::Arc;

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
    /// List of nodes of the optimal trajectory
    pub nodes: Vec<Node>,
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
    pub all_dvs: Vec<Vector3<f64>>,
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
    /// NOTE: this may cause some nodes to be below the surface of a celestial object if in low orbit
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
        if direction.norm() < 2e-16 {
            return Err(NyxError::TargetsTooClose);
        }
        let distance_increment = direction.norm() / (node_count as f64);
        let duration_increment = (xf.epoch() - x0.epoch()) / (node_count as f64);
        direction /= direction.norm();

        // Build each node successively (includes xf)
        let mut nodes = Vec::with_capacity(node_count + 1);
        let mut prev_node_radius = x0.orbit.radius();
        let mut prev_node_epoch = x0.epoch();

        for _ in 0..node_count {
            // Compute the position we want.
            let this_node = prev_node_radius + distance_increment * direction;
            let this_epoch = prev_node_epoch + duration_increment;
            nodes.push(Node {
                x: this_node[0],
                y: this_node[1],
                z: this_node[2],
                frame: x0.orbit.frame,
                epoch: this_epoch,
            });
            prev_node_radius = this_node;
            prev_node_epoch = this_epoch;
        }
        Ok(Self {
            prop,
            nodes,
            x0,
            xf,
            current_iteration: 0,
            max_iterations: 50,
            improvement_threshold: 0.01,
            all_dvs: Vec::with_capacity(node_count),
        })
    }

    /// Builds a multiple shooting structure assuming that the optimal trajectory is near a linear
    /// heuristic in geodetic altitude and direction.
    /// For example, if x0 has an altitude of 100 km and xf has an altitude
    /// of 200 km, and 10 nodes are required over 10 minutes, then node 1 will be 110 km, node 2 220km, etc.
    /// body_frame must be a body fixed frame
    pub fn linear_altitude_heuristic(
        x0: Spacecraft,
        xf: Orbit,
        node_count: usize,
        body_frame: Frame,
        prop: &'a Propagator<'a, D, E>,
        cosm: Arc<Cosm>,
    ) -> Result<Self, NyxError> {
        if node_count < 3 {
            error!("At least three nodes are needed for a multiple shooting optimization");
            return Err(NyxError::UnderdeterminedProblem);
        }

        if !body_frame.is_body_fixed() {
            return Err(NyxError::TargetError(
                "Body frame is not body fixed".to_string(),
            ));
        }

        let delta_t = xf.epoch() - x0.epoch();
        let xf_bf = cosm.frame_chg(&xf, body_frame);

        let duration_increment = (xf.epoch() - x0.epoch()) / (node_count as f64);

        let (_, traj) = prop.with(x0).for_duration_with_traj(delta_t)?;

        // Build each node successively (includes xf)
        let mut nodes = Vec::with_capacity(node_count + 1);
        let mut prev_node_epoch = x0.epoch();

        let inertial_frame = x0.orbit.frame;
        for i in 0..node_count {
            // Compute the position we want.
            let this_epoch = prev_node_epoch + duration_increment;
            let orbit_point = traj.at(this_epoch)?.orbit;
            // Convert this orbit into the body frame
            let orbit_point_bf = cosm.frame_chg(&orbit_point, body_frame);
            // Note that the altitude here might be different, so we scale the altitude change by the current altitude
            let desired_alt_i = (xf_bf.geodetic_height() - orbit_point_bf.geodetic_height())
                / ((node_count - i) as f64).sqrt();
            // Build the node in the body frame and convert that to the original frame
            let node_bf = Orbit::from_geodesic(
                orbit_point_bf.geodetic_latitude(),
                orbit_point_bf.geodetic_longitude(),
                orbit_point_bf.geodetic_height() + desired_alt_i,
                this_epoch,
                body_frame,
            );
            // Convert that back into the inertial frame
            let this_node = cosm.frame_chg(&node_bf, inertial_frame).radius();
            nodes.push(Node {
                x: this_node[0],
                y: this_node[1],
                z: this_node[2],
                frame: inertial_frame,
                epoch: this_epoch,
            });
            prev_node_epoch = this_epoch;
        }
        Ok(Self {
            prop,
            nodes,
            x0,
            xf,
            current_iteration: 0,
            max_iterations: 100,
            improvement_threshold: 0.01,
            all_dvs: Vec::with_capacity(node_count),
        })
    }

    /// Solve the multiple shooting problem by finding the arrangement of nodes to minimize the cost function.
    pub fn solve(&mut self, cost: CostFunction) -> Result<MultipleShootingSolution, NyxError> {
        let mut prev_cost = 1e12; // We don't use infinity because we compare a ratio of cost
        for it in 0..self.max_iterations {
            let mut initial_states = Vec::with_capacity(self.nodes.len());
            initial_states.push(self.x0);
            let mut outer_jacobian =
                DMatrix::from_element(3 * self.nodes.len(), 3 * (self.nodes.len() - 1), 0.0);
            let mut cost_vec = DVector::from_element(3 * self.nodes.len(), 0.0);

            // Reset the all_dvs
            self.all_dvs = Vec::with_capacity(self.all_dvs.len());

            for i in 0..self.nodes.len() {
                /* ***
                 ** 1. Solve the delta-v differential corrector between each node
                 ** *** */
                let tgt = Targeter::delta_v(self.prop, self.nodes[i].to_targeter_objective());
                let sol = match tgt.try_achieve_dual(
                    initial_states[i],
                    initial_states[i].epoch(),
                    self.nodes[i].epoch,
                ) {
                    Ok(sol) => sol,
                    Err(e) => return Err(NyxError::MultipleShootingTargeter(i, Box::new(e))),
                };

                let nominal_delta_v =
                    Vector3::new(sol.correction[0], sol.correction[1], sol.correction[2]);

                self.all_dvs.push(nominal_delta_v);
                // Store the Δv and the initial state for the next targeter.
                initial_states.push(sol.achieved_state);
            }
            // NOTE: We have two separate loops because we need the initial state of node i+2 for the dv computation
            // of the third entry to the outer jacobian.
            for i in 0..(self.nodes.len() - 1) {
                /* ***
                 ** 2. Perturb each node and compute the partial of the Δv for the (i-1), i, and (i+1) nodes
                 ** where the partial on the i+1 -th node is just the difference between the velocity at the
                 ** achieved state and the initial state at that node.
                 ** We don't perturb the endpoint node
                 ** *** */

                for axis in 0..3 {
                    /* ***
                     ** 2.A. Perturb the i-th node
                     ** *** */
                    let mut next_node = self.nodes[i].to_targeter_objective();
                    next_node[axis].desired_value += next_node[axis].tolerance;
                    /* ***
                     ** 2.b. Rerun the targeter from the previous node to this one
                     ** Note that because the first initial_state is x0, the i-th "initial state"
                     ** is the initial state to reach the i-th node.
                     ** *** */

                    let inner_tgt_a = Targeter::delta_v(self.prop, next_node.to_vec());
                    let inner_sol_a = match inner_tgt_a.try_achieve_dual(
                        initial_states[i],
                        initial_states[i].epoch(),
                        self.nodes[i].epoch,
                    ) {
                        Ok(sol) => sol,
                        Err(e) => return Err(NyxError::MultipleShootingTargeter(i, Box::new(e))),
                    };

                    // ∂Δv_x / ∂r_x
                    outer_jacobian[(3 * i, 3 * i + axis)] = (inner_sol_a.correction[0]
                        - self.all_dvs[i][0])
                        / next_node[axis].tolerance;
                    // ∂Δv_y / ∂r_x
                    outer_jacobian[(3 * i + 1, 3 * i + axis)] = (inner_sol_a.correction[1]
                        - self.all_dvs[i][1])
                        / next_node[axis].tolerance;
                    // ∂Δv_z / ∂r_x
                    outer_jacobian[(3 * i + 2, 3 * i + axis)] = (inner_sol_a.correction[2]
                        - self.all_dvs[i][2])
                        / next_node[axis].tolerance;

                    /* ***
                     ** 2.C. Rerun the targeter from the new state at the perturbed node to the next unpertubed node
                     ** *** */
                    let inner_tgt_b =
                        Targeter::delta_v(self.prop, self.nodes[i + 1].to_targeter_objective());
                    let inner_sol_b = inner_tgt_b.try_achieve_dual(
                        inner_sol_a.achieved_state,
                        inner_sol_a.achieved_state.epoch(),
                        self.nodes[i + 1].epoch,
                    )?;

                    // Compute the partials wrt the next Δv
                    // ∂Δv_x / ∂r_x
                    outer_jacobian[(3 * (i + 1), 3 * i + axis)] = (inner_sol_b.correction[0]
                        - self.all_dvs[i + 1][0])
                        / next_node[axis].tolerance;
                    // ∂Δv_y / ∂r_x
                    outer_jacobian[(3 * (i + 1) + 1, 3 * i + axis)] = (inner_sol_b.correction[1]
                        - self.all_dvs[i + 1][1])
                        / next_node[axis].tolerance;
                    // ∂Δv_z / ∂r_x
                    outer_jacobian[(3 * (i + 1) + 2, 3 * i + axis)] = (inner_sol_b.correction[2]
                        - self.all_dvs[i + 1][2])
                        / next_node[axis].tolerance;

                    /* ***
                     ** 2.D. Compute the difference between the arrival and departure velocities and node i+1
                     ** *** */
                    if i < self.nodes.len() - 3 {
                        let dv_ip1 = inner_sol_b.achieved_state.orbit.velocity()
                            - initial_states[i + 2].orbit.velocity();
                        // ∂Δv_x / ∂r_x
                        outer_jacobian[(3 * (i + 2), 3 * i + axis)] =
                            dv_ip1[0] / next_node[axis].tolerance;
                        // ∂Δv_y / ∂r_x
                        outer_jacobian[(3 * (i + 2) + 1, 3 * i + axis)] =
                            dv_ip1[1] / next_node[axis].tolerance;
                        // ∂Δv_z / ∂r_x
                        outer_jacobian[(3 * (i + 2) + 2, 3 * i + axis)] =
                            dv_ip1[2] / next_node[axis].tolerance;
                    }
                }
            }

            // Build the cost vector
            for i in 0..self.nodes.len() {
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
                    nodes: self.nodes.clone(),
                    solutions: Vec::with_capacity(self.nodes.len()),
                };
                let mut initial_states = Vec::with_capacity(self.nodes.len());
                initial_states.push(self.x0);

                for (i, node) in self.nodes.iter().enumerate() {
                    // Run the unpertubed targeter
                    let tgt = Targeter::delta_v(self.prop, node.to_targeter_objective());
                    let sol = tgt.try_achieve_dual(
                        initial_states[i],
                        initial_states[i].epoch(),
                        node.epoch,
                    )?;
                    initial_states.push(sol.achieved_state);
                    ms_sol.solutions.push(sol);
                }

                return Ok(ms_sol);
            }

            prev_cost = new_cost;
            // 2. Solve for the next position of the nodes using a pseudo inverse.
            let inv_jac = match pseudo_inverse(&outer_jacobian, NyxError::SingularJacobian) {
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
                let component_no = i % 3;
                match component_no {
                    0 => self.nodes[node_no].x += val,
                    1 => self.nodes[node_no].y += val,
                    2 => self.nodes[node_no].z += val,
                    _ => unreachable!(),
                }
            }
            self.current_iteration += 1;
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

        for (i, node) in self.nodes.iter().enumerate() {
            let dv = match self.all_dvs.get(i) {
                Some(dv) => dv.clone(),
                None => Vector3::zeros(),
            };
            nodemsg.push_str(&format!(
                "[{:.3}, {:.3}, {:.3}, {}, {}, {}, {}, {}, {}],\n",
                node.x,
                node.y,
                node.z,
                self.current_iteration,
                dv[0],
                dv[1],
                dv[2],
                dv.norm(),
                i + 1
            ));
        }
        write!(f, "{}", nodemsg)
    }
}

#[derive(Clone, Debug)]
pub struct MultipleShootingSolution {
    pub x0: Spacecraft,
    pub xf: Orbit,
    pub nodes: Vec<Node>,
    pub solutions: Vec<TargeterSolution>,
}

impl fmt::Display for MultipleShootingSolution {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        for sol in &self.solutions {
            write!(f, "{}", sol)?;
        }
        Ok(())
    }
}

impl MultipleShootingSolution {
    /// Allows building the trajectories between different nodes
    /// This will rebuild the targeters and apply the solutions sequentially
    pub fn build_trajectories<'a, D: Dynamics<StateType = Spacecraft>, E: ErrorCtrl>(
        &self,
        prop: &'a Propagator<'a, D, E>,
    ) -> Result<Vec<ScTraj>, NyxError>
    where
        DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
            + Allocator<f64, <D::StateType as State>::VecLength>,
    {
        let mut trajz = Vec::with_capacity(self.nodes.len());

        for (i, node) in self.nodes.iter().enumerate() {
            let (_, traj) = Targeter::delta_v(prop, node.to_targeter_objective())
                .apply_with_traj(&self.solutions[i])?;
            trajz.push(traj);
        }

        Ok(trajz)
    }
}
