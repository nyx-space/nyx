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
    /// Returns the flatenning of this celestial object.
    fn flatenning() -> f64;
}

/// `NAIF` represents an object which has a NAIF ID and can be loaded from an SPK file.
pub trait NAIF {
    /// Returns the NAIF ID of this object.
    fn id() -> i32;

    /// Returns the NAIF ID of the barycenter
    fn barycenter() -> i32 {
        Self::id() / 100
    }

    /// Returns the NAIF ID of the satellite of this planet given its position (the first satellite is number **one** (not zero)).
    ///
    /// **Warning:** the logic here is based entirely on the NAIF ID of the main planet. No guarantee is made that the requested satellite number exists,
    /// or is present in the loaded SPK files.
    fn satellite(pos: i32) -> i32 {
        Self::id() - 99 + pos
    }
}

// Re-Export the planets
mod planets;
pub use self::planets::*;

// Re-Export state
mod state;
pub use self::state::State;

// Re-Export frames
mod frames;
pub use self::frames::*;
