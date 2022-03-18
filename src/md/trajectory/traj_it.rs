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

use super::{InterpState, Traj};
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::time::TimeSeries;

pub struct TrajIterator<'a, S: InterpState>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    pub time_series: TimeSeries,
    /// A shared pointer to the original trajectory.
    pub traj: &'a Traj<S>,
}

impl<S: InterpState> Iterator for TrajIterator<'_, S>
where
    DefaultAllocator:
        Allocator<f64, S::VecLength> + Allocator<f64, S::Size> + Allocator<f64, S::Size, S::Size>,
{
    type Item = S;

    fn next(&mut self) -> Option<Self::Item> {
        match self.time_series.next() {
            Some(next_epoch) => match self.traj.at(next_epoch) {
                Ok(item) => Some(item),
                Err(e) => {
                    if next_epoch >= self.traj.first().epoch()
                        && next_epoch <= self.traj.last().epoch()
                    {
                        if log_enabled!(log::Level::Error) {
                            error!("[!!!] BUG [!!!]");
                            error!("[!!!]\t{}\t[!!!]", e);
                            error!("[!!!]\t{}\t[!!!]", self.traj);
                            error!("[!!!]     [!!!]");
                        } else {
                            println!("[!!!] BUG [!!!]");
                            println!("[!!!]\t{}\t[!!!]", e);
                            println!("[!!!]\t{}\t[!!!]", self.traj);
                            println!("[!!!]     [!!!]");
                        };
                    }
                    None
                }
            },
            None => None,
        }
    }
}
