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
use crate::md::StateParameter;
use crate::State;
use rand_distr::{Distribution, Normal};

/// A dispersions configuration, allows specifying min/max bounds (by default, they are not set)
#[derive(Copy, Clone)]
pub struct Dispersion<Distr: Distribution<f64> + Copy> {
    pub param: StateParameter,
    pub distr: Distr,
    pub bound_min: Option<f64>,
    pub bound_max: Option<f64>,
}

impl<Distr: Distribution<f64> + Copy> Dispersion<Distr> {
    pub fn new(param: StateParameter, distr: Distr) -> Self {
        Self {
            param,
            distr,
            bound_min: None,
            bound_max: None,
        }
    }
}

impl Dispersion<Normal<f64>> {
    /// Initializes a new normal dispersion of zero mean from the 1σ
    pub fn from_std_dev(param: StateParameter, std_dev: f64) -> Self {
        Self::new(param, Normal::new(0.0, std_dev).unwrap())
    }

    /// Initializes a new normal dispersion of zero mean from the 3σ
    pub fn from_3std_dev(param: StateParameter, std_dev: f64) -> Self {
        Self::new(param, Normal::new(0.0, std_dev / 3.0).unwrap())
    }
}

/// A dispersed state
#[derive(Clone)]
pub struct DispersedState<S: State>
where
    DefaultAllocator: Allocator<S::Size> + Allocator<S::Size, S::Size> + Allocator<S::VecLength>,
{
    /// The dispersed state
    pub state: S,
    /// The dispersions applied to the template state (template state + self.actual_dispersions = self.state)
    pub actual_dispersions: Vec<(StateParameter, f64)>,
}
