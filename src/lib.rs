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

[Nyx](https://en.wikipedia.org/wiki/Nyx): Blazing fast high-fidelity astrodynamics for Monte Carlo analyzes of constellations, interplanetary missions, and deep space flight navigation.

Refer to [nyxspace.com](https://nyxspace.com) for a user guide, a show case, the MathSpec, and the validation data.
*/

// Allow confusable identifiers, as the code tries to use the literature's notation where possible.
#![allow(confusable_idents)]

extern crate log;

/// Provides all the propagators / integrators available in `nyx`.
pub mod propagators;

/// Provides several dynamics used for orbital mechanics and attitude dynamics, which can be elegantly combined.
pub mod dynamics;

/// Provides the solar system planets, and state and ephemerides management.
pub mod cosmic;

/// Utility functions shared by different modules, and which may be useful to engineers.
pub mod utils;

pub mod errors;
/// Nyx will (almost) never panic and functions which may fail will return an error.
pub use self::errors::NyxError;

/// All the input/output needs for this library, including loading of SPICE kernels, and gravity potential files.
pub mod io;

/// All the orbital determination and spacecraft navigation tools and functions.
pub mod od;

/// All of the mission design and mission analysis tools and functions
pub mod md;

/// Simple tools (e.g. Lambert solver)
pub mod tools;

/// Monte Carlo module
pub mod mc;

/// Polynomial and fitting module
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
