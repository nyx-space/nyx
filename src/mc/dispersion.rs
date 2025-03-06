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

use crate::md::StateParameter;
use serde::{Deserialize, Serialize};
use typed_builder::TypedBuilder;

/// A dispersions configuration, allows specifying min/max bounds (by default, they are not set)
#[derive(Copy, Clone, Debug, TypedBuilder, Serialize, Deserialize)]
pub struct StateDispersion {
    pub param: StateParameter,
    #[builder(default, setter(strip_option))]
    pub mean: Option<f64>,
    #[builder(default, setter(strip_option))]
    pub std_dev: Option<f64>,
}

impl StateDispersion {
    pub fn zero_mean(param: StateParameter, std_dev: f64) -> Self {
        Self {
            param,
            std_dev: Some(std_dev),
            mean: Some(0.0),
        }
    }
}
