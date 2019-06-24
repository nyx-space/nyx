extern crate nalgebra as na;
/// `CelestialBody` represents a celestial body.
///
/// Note that all planets are defined as types. This leverages higher speed of execution via monomorphism.
/// The `CelestialBody`s provided in nyx use the same values as those in [GMAT 2016a](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/GmatDefaults.hpp).
/// NOTE: There is no Pluto defined in nyx because it isn't a planet: it's a collection of three (four?) small rocks orbiting each other.
pub trait CelestialBody {
    /// Returns the gravitional parameter of the given body. **Unit**: km<sup>3</sup>/s<sup>2</sup>
    fn gm() -> f64;
    /// Returns the equatorial radius of this celestial object.
    fn eq_radius() -> f64;
    /// Returns the flattening of this celestial object.
    fn flattening() -> f64;
}

// Re-Export state
mod state;
pub use self::state::*;

// Re-Export frames
mod frames;
pub use self::frames::*;

mod axb;
mod cosm;
pub use self::cosm::*;
mod exb;
mod fxb;
mod hermite;
