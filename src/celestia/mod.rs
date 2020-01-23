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
    /// Solar System Barycenter
    pub const SSB: i32 = 0;
    /// Sun center ID
    pub const SUN: i32 = 10;
    /// Mercury barycenter ID
    pub const MERCURY_BARYCENTER: i32 = 1;
    /// Mercury center ID
    pub const MERCURY: i32 = 1;
    /// Venus barycenter ID
    pub const VENUS_BARYCENTER: i32 = 2;
    /// Venus center ID
    pub const VENUS: i32 = 2;
    /// Earth barycenter ID
    pub const EARTH_BARYCENTER: i32 = 3;
    /// Earth planet ID
    pub const EARTH: i32 = 399;
    /// Earth's MOon planet ID
    pub const EARTH_MOON: i32 = 301;
    /// Mars barycenter ID
    pub const MARS_BARYCENTER: i32 = 4;
    /// Jupiter barycenter ID
    pub const JUPITER_BARYCENTER: i32 = 5;
    /// Saturn barycenter ID
    pub const SATURN_BARYCENTER: i32 = 6;
    /// Uranus barycenter ID
    pub const URANUS_BARYCENTER: i32 = 7;
    /// Neptune barycenter ID
    pub const NEPTUNE_BARYCENTER: i32 = 8;
}

// Re-Export state
mod state;
pub use self::state::*;

// Re-Export frames
mod frames;
pub use self::frames::*;

mod cosm;
mod xb;
pub use self::cosm::*;
mod hermite;

/// The eclipse module allows finding eclipses and (conversely) visibility between a state and another one (e.g. a planet or the Sun).
pub mod eclipse;

/// Speed of light in meters per second
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;
/// Speed of light in kilometers per second
pub const SPEED_OF_LIGHT_KMS: f64 = 299_792.458;

/// Astronomical unit, in kilometers, according to the [IAU](https://www.iau.org/public/themes/measuring/).
pub const AU: f64 = 149_597_870.700;
