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

use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::prelude::{Interpolatable, Traj};
use crate::time::Epoch;
use anise::analysis::specs::StateSpecTrait;
use anise::analysis::AnalysisError;
use anise::astro::orbit::Orbit;
use anise::prelude::{Aberration, Almanac};

impl<S: Interpolatable> StateSpecTrait for Traj<S>
where
    DefaultAllocator: Allocator<S::VecLength> + Allocator<S::Size> + Allocator<S::Size, S::Size>,
{
    fn ab_corr(&self) -> Option<Aberration> {
        None
    }

    fn evaluate(&self, epoch: Epoch, _almanac: &Almanac) -> Result<Orbit, AnalysisError> {
        self.at(epoch)
            .map(|state| state.orbit())
            .map_err(|e| AnalysisError::GenericAnalysisError {
                err: format!("{e}"),
            })
    }
}
