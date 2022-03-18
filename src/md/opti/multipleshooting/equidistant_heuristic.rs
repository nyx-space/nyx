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

use super::ctrlnodes::Node;
use super::multishoot::MultipleShooting;
pub use super::CostFunction;
use crate::errors::TargetingError;
use crate::md::ui::*;
use crate::propagators::error_ctrl::ErrorCtrl;
use crate::{Orbit, Spacecraft};

impl<'a, E: ErrorCtrl> MultipleShooting<'a, E, Node, 3, 3> {
    /// Builds a multiple shooting structure assuming that the optimal trajectory is a straight line
    /// between the start and end points. The position of the nodes will be update at each iteration
    /// of the outer loop.
    /// NOTE: this may cause some nodes to be below the surface of a celestial object if in low orbit
    pub fn equidistant_nodes(
        x0: Spacecraft,
        xf: Orbit,
        node_count: usize,
        prop: &'a Propagator<'a, SpacecraftDynamics, E>,
    ) -> Result<Self, NyxError> {
        if node_count < 3 {
            error!("At least three nodes are needed for a multiple shooting optimization");
            return Err(NyxError::Targeter(TargetingError::UnderdeterminedProblem));
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
                vmag: 0.0,
                frame: x0.orbit.frame,
                epoch: this_epoch,
            });
            prev_node_radius = this_node;
            prev_node_epoch = this_epoch;
        }
        Ok(Self {
            prop,
            targets: nodes,
            x0,
            xf,
            current_iteration: 0,
            max_iterations: 50,
            improvement_threshold: 0.01,
            variables: [
                Vary::VelocityX.try_into().unwrap(),
                Vary::VelocityY.try_into().unwrap(),
                Vary::VelocityZ.try_into().unwrap(),
            ],
            all_dvs: Vec::with_capacity(node_count),
        })
    }
}
