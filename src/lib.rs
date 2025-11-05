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

/*! # nyx-space

[Nyx](https://en.wikipedia.org/wiki/Nyx): High-fidelity, fast astrodynamics toolkit for spacecraft navigation, mission design, and Monte Carlo analyses.

Refer to the [user guide](https://nyxspace.com) for more information.
*/

// Allow confusable identifiers, as the code tries to use the literature's notation where possible.
#![allow(confusable_idents)]

extern crate log;

/// Numerical propagators and integrators.
pub mod propagators;

/// Orbital and attitude dynamics which can be elegantly combined.
pub mod dynamics;

/// Solar system bodies, state, and ephemerides management.
pub mod cosmic;

/// Shared utilities.
pub mod utils;

pub mod errors;
/// Nyx's error handling.
pub use self::errors::NyxError;

/// I/O for SPICE kernels, gravity models, and more.
pub mod io;

/// Orbit determination and spacecraft navigation.
pub mod od;

/// Mission design and analysis.
pub mod md;

/// Standalone astrodynamics tools.
pub mod tools;

/// Monte Carlo simulations.
pub mod mc;

/// Polynomial fitting and interpolation.
pub mod polyfit;

/// Re-export of hifitime
pub mod time {
    pub use hifitime::prelude::*;
}

/// Re-export nalgebra
pub mod linalg {
    pub use nalgebra::base::*;
}

/// Re-export some useful things
pub use self::cosmic::{Orbit, Spacecraft, State, TimeTagged};

/// The GMAT Earth gravitation parameter, used only for testing.
#[cfg(test)]
pub(crate) const GMAT_EARTH_GM: f64 = 398_600.441_5;
