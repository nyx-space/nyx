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

/// Known orientation IDs defined for ease of access. All Cosm objects may be accessed via Cosm directly.
pub mod orientations {
    /// J2000 orientation frame
    pub const J2000: i32 = 1;
}

/// Known planets IDs defined for ease of access. All Cosm objects may be accessed via Cosm directly.
pub mod bodies {
    /// Sun center ID
    pub const SUN: i32 = 0;
    /// Mercury center ID
    pub const MERCURY: i32 = 1;
    /// Venus center ID
    pub const VENUS: i32 = 2;
    /// Earth center ID
    pub const EARTH: i32 = 3;
    /// Mars center ID
    pub const MARS: i32 = 4;
    /// Jupiter center ID
    pub const JUPITER: i32 = 5;
    /// Saturn center ID
    pub const SATURN: i32 = 6;
    /// Uranus center ID
    pub const URANUS: i32 = 7;
    /// Neptune center ID
    pub const NEPTUNE: i32 = 8;
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
