use super::{CelestialBody, NAIF};

/// Planet Mercury as defined in [GMAT 2016a](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/GmatDefaults.hpp).
/// **Warning**: Keplerian dynamics are _not_ a correct representation of the orbit of Mercury
/// (cf. [this discussion](https://physics.stackexchange.com/questions/26408/what-did-general-relativity-clarify-about-mercury))
/// so one should take into account that general relativity is required for high fidelity dynamics in the vicinity of this planet.
pub struct MERCURY;

impl CelestialBody for MERCURY {
    fn gm() -> f64 {
        22_032.080486418
    }
    fn eq_radius() -> f64 {
        2439.7
    }
    fn flatenning() -> f64 {
        0.0
    }
}

impl NAIF for MERCURY {
    fn id() -> i32 {
        199
    }
}

/// Planet Venus as defined in [GMAT 2016a](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/GmatDefaults.hpp).
pub struct VENUS;

impl CelestialBody for VENUS {
    fn gm() -> f64 {
        324_858.59882646
    }
    fn eq_radius() -> f64 {
        6051.9
    }
    fn flatenning() -> f64 {
        0.0
    }
}

impl NAIF for VENUS {
    fn id() -> i32 {
        299
    }
}

/// Planet Earth as defined in [GMAT 2016a](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/GmatDefaults.hpp).
pub struct EARTH;

impl EARTH {
    /// Defines the semi major radius of the ellipsoid of Earth, as per WGS84, in km.
    pub fn semi_major_radius() -> f64 {
        6378.1370
    }

    /// The rotation rate of Earth, in radians per seconds; [source](http://hpiers.obspm.fr/eop-pc/models/constants.html).
    pub fn rotation_rate() -> f64 {
        7.292_115_146_706_4e-5
    }
}

impl CelestialBody for EARTH {
    fn gm() -> f64 {
        398_600.4415
    }
    fn eq_radius() -> f64 {
        6378.1363
    }
    fn flatenning() -> f64 {
        // From [EMG2008](http://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/egm08_wgs84.html)
        0.0033528106647474805
    }
}

impl NAIF for EARTH {
    fn id() -> i32 {
        399
    }
}

/// Planet Mars as defined in [GMAT 2016a](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/GmatDefaults.hpp).
pub struct MARS;

impl CelestialBody for MARS {
    fn gm() -> f64 {
        42_828.314258067
    }
    fn eq_radius() -> f64 {
        3.397e3
    }
    fn flatenning() -> f64 {
        0.00647630
    }
}

impl NAIF for MARS {
    fn id() -> i32 {
        499
    }
}

/// Planet Jupiter as defined in [GMAT 2016a](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/GmatDefaults.hpp).
pub struct JUPITER;

impl CelestialBody for JUPITER {
    fn gm() -> f64 {
        126_712_767.85780
    }
    fn eq_radius() -> f64 {
        7.1492e4
    }
    fn flatenning() -> f64 {
        0.06487439
    }
}

impl NAIF for JUPITER {
    fn id() -> i32 {
        599
    }
}

/// Planet Saturn as defined in [GMAT 2016a](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/GmatDefaults.hpp).
pub struct SATURN;

impl CelestialBody for SATURN {
    fn gm() -> f64 {
        37_940_626.061137
    }
    fn eq_radius() -> f64 {
        6.0268e4
    }
    fn flatenning() -> f64 {
        0.09796243
    }
}

impl NAIF for SATURN {
    fn id() -> i32 {
        699
    }
}

/// Planet Uranus as defined in [GMAT 2016a](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/GmatDefaults.hpp).
pub struct URANUS;

impl CelestialBody for URANUS {
    fn gm() -> f64 {
        5_794_549.0070719
    }
    fn eq_radius() -> f64 {
        2.5559e4
    }
    fn flatenning() -> f64 {
        0.02292734
    }
}

impl NAIF for URANUS {
    fn id() -> i32 {
        799
    }
}

/// Planet Neptune as defined in [GMAT 2016a](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/GmatDefaults.hpp).
pub struct NEPTUNE;

impl CelestialBody for NEPTUNE {
    fn gm() -> f64 {
        6_836_534.0638793
    }
    fn eq_radius() -> f64 {
        2.5269e4
    }
    fn flatenning() -> f64 {
        0.01856029
    }
}

impl NAIF for NEPTUNE {
    fn id() -> i32 {
        899
    }
}
