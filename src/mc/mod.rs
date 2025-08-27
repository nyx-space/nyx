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

pub use rand::prelude::*;
pub use rand_distr::{Distribution, Normal, Uniform};
pub use rand_pcg::Pcg64Mcg;

pub mod helpers;
mod montecarlo;

pub use montecarlo::MonteCarlo;

mod dispersion;
pub use dispersion::StateDispersion;

mod generator;
pub use generator::{DispersedState, Dispersion};

mod multivariate;
pub use multivariate::MvnSpacecraft;

mod results;
pub use results::{Results, Stats};
