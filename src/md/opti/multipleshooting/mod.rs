/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use anise::errors::{AlmanacError, PhysicsError};
use snafu::Snafu;

use crate::md::{trajectory::TrajError, TargetingError};

pub mod altitude_heuristic;
pub mod ctrlnodes;
pub mod equidistant_heuristic;
pub mod multishoot;

/// Built-in cost functions to minimize
#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum CostFunction {
    /// J = ∫ \vec{u}^T\vec{u} dt
    MinimumEnergy,
    /// J = ∫ |\vec{u}| dt -- Warning, this may lead to loads to bang-coast-bang solutions
    MinimumFuel,
}

#[derive(Debug, Snafu)]
pub enum MultipleShootingError {
    #[snafu(display("segment #{segment} encountered {source}"))]
    TargetingError {
        segment: usize,
        source: TargetingError,
    },
    #[snafu(display("during a multiple shooting, encountered {source}"))]
    MultiShootTrajError { source: TrajError },
    #[snafu(display("duration a multiple shoot, issue due to Almanac: {action} {source}"))]
    MultiShootAlmanacError {
        source: AlmanacError,
        action: &'static str,
    },
    #[snafu(display("duration a multiple shoot, physics issue:  {source}"))]
    MultiShootPhysicsError { source: PhysicsError },
}
