/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName};
use crate::md::trajectory::Interpolatable;
pub use crate::od::estimate::*;
pub use crate::od::*;
use msr::sensitivity::TrackerSensitivity;
use std::fmt;
use std::ops::Add;

use super::ODSolution;

impl<StateType, EstType, MsrSize, Trk> fmt::Display for ODSolution<StateType, EstType, MsrSize, Trk>
where
    StateType: Interpolatable + Add<OVector<f64, <StateType as State>::Size>, Output = StateType>,
    EstType: Estimate<StateType>,
    MsrSize: DimName,
    Trk: TrackerSensitivity<StateType, StateType>,
    <DefaultAllocator as Allocator<<StateType as State>::VecLength>>::Buffer<f64>: Send,
    DefaultAllocator: Allocator<<StateType as State>::Size>
        + Allocator<<StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <StateType as State>::Size>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<StateType as State>::Size, <StateType as State>::Size>
        + Allocator<<StateType as State>::Size, MsrSize>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.results().count() == 0 {
            writeln!(f, "Empty OD solution")
        } else {
            let kind = if self.is_filter_run() {
                "Filter"
            } else {
                "Smoother"
            };
            writeln!(
                f,
                "{kind} OD solution from {:?}, spanning {}-{}, with {} estimates, including {} accepted residuals",
                self.unique(),
                self.estimates.first().unwrap().epoch(),
                self.estimates.last().unwrap().epoch(),
                self.estimates.len(),
                self.accepted_residuals().len()
            )
        }
    }
}
