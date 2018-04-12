/// `Body` represents a celestial body.
pub trait Body {
    /// Returns the gravitional parameter of the given body. **Unit**: km<sup>3</sup>/s<sup>2</sup>
    fn gm(&self) -> f64;
}

/// `CustomBody` can be used to define a custom celestial object, such as an asteroid or a moon
/// which isn't provided by nyx or one of the ephemerides in use.
pub struct CustomBody {
    mu: f64,
}

impl CustomBody {
    /// Initialize a custom celestial body from the provided gravitional parameter in km<sup>3</sup>/s<sup>2</sup>.
    pub fn from_gm(mu: f64) -> CustomBody {
        CustomBody { mu: mu }
    }
}

impl Body for CustomBody {
    fn gm(&self) -> f64 {
        self.mu
    }
}

/// The `Planet`s provided in nyx use the same values as those in [GMAT 2016a](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/GmatDefaults.hpp).
/// NOTE: There is no Pluto defined in nyx because it isn't a planet: it's a collection of three (four?) small rocks orbiting each other.
pub struct Planet {
    pub mu: f64,
    pub eq_radius: f64,
    pub flatenning: f64,
    pub naif_id: u16,
}

impl Body for Planet {
    fn gm(&self) -> f64 {
        self.mu
    }
}

lazy_static! {
    pub static ref MERCURY: Planet = Planet{mu: 22032.080486418, flatenning: 0.0,
                                            eq_radius: 2.4397e3, naif_id: 199};
    pub static ref VENUS: Planet = Planet{mu: 324858.59882646, flatenning: 0.0,
                                            eq_radius: 6.0519e3, naif_id: 299};
    pub static ref EARTH: Planet = Planet{mu: 398_600.4415, flatenning: 0.00335270,
                                            eq_radius: 6.3781363E3, naif_id: 399};
    pub static ref MARS: Planet = Planet{mu: 42828.314258067, flatenning: 0.00647630,
                                            eq_radius: 3.397e3, naif_id: 499};
    pub static ref JUPITER: Planet = Planet{mu: 126712767.85780, flatenning: 0.06487439,
                                            eq_radius: 7.1492e4, naif_id: 599};
    pub static ref SATURN: Planet = Planet{mu: 37940626.061137, flatenning: 0.09796243,
                                            eq_radius: 6.0268e4, naif_id: 699};
    pub static ref URANUS: Planet = Planet{mu: 5794549.0070719, flatenning: 0.02292734,
                                            eq_radius: 2.5559e4, naif_id: 799};
    pub static ref NEPTUNE: Planet = Planet{mu: 6836534.0638793, flatenning: 0.01856029,
                                            eq_radius: 2.5269e4, naif_id: 899};
}
