extern crate nalgebra as na;

use super::time::Epoch;

/// A trait allowing for something to have an epoch
pub trait TimeTagged {
    /// Retrieve the Epoch
    fn epoch(&self) -> Epoch;
    /// Set the Epoch
    fn set_epoch(&mut self, dt: Epoch);
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

mod rotations;
pub use self::rotations::*;

mod cosm;
mod xb;
pub use self::cosm::*;

/// The eclipse module allows finding eclipses and (conversely) visibility between a state and another one (e.g. a planet or the Sun).
pub mod eclipse;

/// Speed of light in meters per second
pub const SPEED_OF_LIGHT: f64 = 299_792_458.0;
/// Speed of light in kilometers per second
pub const SPEED_OF_LIGHT_KMS: f64 = 299_792.458;

/// Astronomical unit, in kilometers, according to the [IAU](https://www.iau.org/public/themes/measuring/).
pub const AU: f64 = 149_597_870.700;
