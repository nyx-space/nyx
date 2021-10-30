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

use crate::cosmic::Frame;
use crate::md::targeter::Objective;
use crate::md::ui::StateParameter;
use crate::time::Epoch;

#[derive(Copy, Clone, Debug)]
pub struct Node {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub frame: Frame,
    pub epoch: Epoch,
}

impl Node {
    pub fn to_targeter_objective(&self) -> Vec<Objective> {
        return vec![
            Objective::new(StateParameter::X, self.x),
            Objective::new(StateParameter::Y, self.y),
            Objective::new(StateParameter::Z, self.z),
        ];
    }

    pub fn rmag(&self) -> f64 {
        (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt()
    }
}
