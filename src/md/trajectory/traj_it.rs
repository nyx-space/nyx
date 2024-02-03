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

use super::Traj;
use crate::time::TimeSeries;
use crate::Spacecraft;

pub struct TrajIterator<'a> {
    pub time_series: TimeSeries,
    /// A shared pointer to the original trajectory.
    pub traj: &'a Traj,
}

impl Iterator for TrajIterator<'_> {
    type Item = Spacecraft;

    fn next(&mut self) -> Option<Self::Item> {
        match self.time_series.next() {
            Some(next_epoch) => match self.traj.at(next_epoch) {
                Ok(item) => Some(item),
                Err(e) => {
                    if next_epoch >= self.traj.first().epoch()
                        && next_epoch <= self.traj.last().epoch()
                    {
                        let msg = format!(
                            "{e} out of bounds in {}! Please submit bug report with exported traj",
                            self.traj
                        );
                        if log_enabled!(log::Level::Error) {
                            error!("{msg}");
                        } else {
                            eprintln!("{msg}");
                        };
                    }
                    None
                }
            },
            None => None,
        }
    }
}
