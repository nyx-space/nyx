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

use super::ctrlnodes::Node;
use super::multishoot::MultipleShooting;
pub use super::CostFunction;
use super::{
    MultiShootAlmanacSnafu, MultiShootPhysicsSnafu, MultiShootTrajSnafu, MultipleShootingError,
    TargetingSnafu,
};
use crate::errors::TargetingError;
use crate::md::{prelude::*, PropSnafu};
use crate::{Orbit, Spacecraft};
use log::error;
use snafu::ResultExt;

impl<'a> MultipleShooting<'a, Node, 3, 3> {
    /// Builds a multiple shooting structure assuming that the optimal trajectory is near a linear
    /// heuristic in geodetic altitude and direction.
    /// For example, if x0 has an altitude of 100 km and xf has an altitude
    /// of 200 km, and 10 nodes are required over 10 minutes, then node 1 will be 110 km, node 2 220km, etc.
    /// body_frame must be a body fixed frame
    pub fn linear_altitude_heuristic(
        x0: Spacecraft,
        xf: Orbit,
        node_count: usize,
        angular_velocity_deg_s: f64,
        body_frame: Frame,
        prop: &'a Propagator<SpacecraftDynamics>,
        almanac: Arc<Almanac>,
    ) -> Result<Self, MultipleShootingError> {
        if node_count < 3 {
            error!("At least three nodes are needed for a multiple shooting optimization");
            return Err(MultipleShootingError::TargetingError {
                segment: 0_usize,
                source: Box::new(TargetingError::UnderdeterminedProblem),
            });
        }

        let delta_t = xf.epoch - x0.epoch();
        let xf_bf = almanac
            .transform_to(xf, body_frame, None)
            .context(MultiShootAlmanacSnafu {
                action: "converting node into the body frame",
            })?;

        let duration_increment = (xf.epoch - x0.epoch()) / (node_count as f64);

        let (_, traj) = prop
            .with(x0, almanac.clone())
            .for_duration_with_traj(delta_t)
            .context(PropSnafu)
            .context(TargetingSnafu { segment: 0_usize })?;

        // Build each node successively (includes xf)
        let mut nodes = Vec::with_capacity(node_count + 1);
        let mut prev_node_epoch = x0.epoch();

        let inertial_frame = x0.orbit.frame;
        for i in 0..node_count {
            // Compute the position we want.
            let this_epoch = prev_node_epoch + duration_increment;
            let orbit_point = traj.at(this_epoch).context(MultiShootTrajSnafu)?.orbit;
            // Convert this orbit into the body frame
            let orbit_point_bf = almanac
                .clone()
                .transform_to(orbit_point, body_frame, None)
                .context(MultiShootAlmanacSnafu {
                    action: "converting node into the body frame",
                })?;

            // Note that the altitude here might be different, so we scale the altitude change by the current altitude
            let desired_alt_i = (xf_bf.height_km().context(MultiShootPhysicsSnafu)?
                - orbit_point_bf.height_km().context(MultiShootPhysicsSnafu)?)
                / ((node_count - i) as f64).sqrt();
            // Build the node in the body frame and convert that to the original frame
            let node_bf = Orbit::try_latlongalt(
                orbit_point_bf
                    .latitude_deg()
                    .context(MultiShootPhysicsSnafu)?,
                orbit_point_bf.longitude_deg(),
                orbit_point_bf.height_km().context(MultiShootPhysicsSnafu)? + desired_alt_i,
                angular_velocity_deg_s,
                this_epoch,
                body_frame,
            )
            .context(MultiShootPhysicsSnafu)?;

            // Convert that back into the inertial frame
            let this_node = almanac
                .transform_to(node_bf, inertial_frame, None)
                .context(MultiShootAlmanacSnafu {
                    action: "converting node back into the inertial frame",
                })?
                .radius_km;

            nodes.push(Node {
                x: this_node[0],
                y: this_node[1],
                z: this_node[2],
                vmag: 0.0,
                frame: inertial_frame,
                epoch: this_epoch,
            });
            prev_node_epoch = this_epoch;
        }
        Ok(Self {
            prop,
            targets: nodes,
            x0,
            xf,
            current_iteration: 0,
            max_iterations: 100,
            improvement_threshold: 0.01,
            variables: [
                Vary::VelocityX.into(),
                Vary::VelocityY.into(),
                Vary::VelocityZ.into(),
            ],
            all_dvs: Vec::with_capacity(node_count),
        })
    }
}
