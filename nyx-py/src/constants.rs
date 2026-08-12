/*
 * ANISE Toolkit
 * Copyright (C) 2021-onward Christopher Rabotin <christopher.rabotin@gmail.com> et al. (cf. AUTHORS.md)
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *
 * Documentation: https://nyxspace.com/
 */

use anise::constants::SPEED_OF_LIGHT_KM_S;
use anise::constants::celestial_objects::*;
use anise::constants::orientations::*;
use anise::constants::usual_planetary_constants::MEAN_EARTH_ANGULAR_VELOCITY_DEG_S;
use anise::constants::usual_planetary_constants::MEAN_MOON_ANGULAR_VELOCITY_DEG_S;
use pyo3::prelude::*;

use anise::constants::frames::*;
use anise::frames::Frame;

#[pyclass]
#[pyo3(module = "anise.constants")]
struct Frames {}

#[pymethods]
impl Frames {
    #[classattr]
    const SSB_J2000: Frame = SSB_J2000;
    #[classattr]
    const MERCURY_J2000: Frame = MERCURY_J2000;
    #[classattr]
    const VENUS_J2000: Frame = VENUS_J2000;
    #[classattr]
    const EARTH_MOON_BARYCENTER_J2000: Frame = EARTH_MOON_BARYCENTER_J2000;
    #[classattr]
    const MARS_BARYCENTER_J2000: Frame = MARS_BARYCENTER_J2000;
    #[classattr]
    const JUPITER_BARYCENTER_J2000: Frame = JUPITER_BARYCENTER_J2000;
    #[classattr]
    const SATURN_BARYCENTER_J2000: Frame = SATURN_BARYCENTER_J2000;
    #[classattr]
    const URANUS_BARYCENTER_J2000: Frame = URANUS_BARYCENTER_J2000;
    #[classattr]
    const NEPTUNE_BARYCENTER_J2000: Frame = NEPTUNE_BARYCENTER_J2000;
    #[classattr]
    const PLUTO_BARYCENTER_J2000: Frame = PLUTO_BARYCENTER_J2000;
    #[classattr]
    const SUN_J2000: Frame = SUN_J2000;
    #[classattr]
    const MOON_J2000: Frame = MOON_J2000;
    #[classattr]
    const EARTH_J2000: Frame = EARTH_J2000;
    #[classattr]
    const EME2000: Frame = EME2000;
    #[classattr]
    const EARTH_ECLIPJ2000: Frame = EARTH_ECLIPJ2000;
    #[classattr]
    const GCRF: Frame = GCRF;
    #[classattr]
    const ICRF: Frame = ICRF;
    #[classattr]
    const IAU_MERCURY_FRAME: Frame = IAU_MERCURY_FRAME;
    #[classattr]
    const IAU_VENUS_FRAME: Frame = IAU_VENUS_FRAME;
    #[classattr]
    const IAU_EARTH_FRAME: Frame = IAU_EARTH_FRAME;
    #[classattr]
    const EARTH_ITRF93: Frame = EARTH_ITRF93;
    #[classattr]
    const EARTH_MOD_FRAME: Frame = EARTH_MOD_FRAME;
    #[classattr]
    const EARTH_TOD_FRAME: Frame = EARTH_TOD_FRAME;
    #[classattr]
    const EARTH_MOD_LEGACY_FRAME: Frame = EARTH_MOD_LEGACY_FRAME;
    #[classattr]
    const EARTH_TOD_LEGACY_FRAME: Frame = EARTH_TOD_LEGACY_FRAME;
    #[classattr]
    const EARTH_TEME_FRAME: Frame = EARTH_TEME_FRAME;
    #[classattr]
    const EARTH_TEME_LEGACY_FRAME: Frame = EARTH_TEME_LEGACY_FRAME;
    #[classattr]
    const MOON_ME_FRAME: Frame = MOON_ME_FRAME;
    #[classattr]
    const MOON_ME_DE421_FRAME: Frame = MOON_ME_DE421_FRAME;
    #[classattr]
    const MOON_ME_DE440_ME421_FRAME: Frame = MOON_ME_DE440_ME421_FRAME;
    #[classattr]
    const MOON_PA_FRAME: Frame = MOON_PA_FRAME;
    #[classattr]
    const MOON_PA_DE421_FRAME: Frame = MOON_PA_DE421_FRAME;
    #[classattr]
    const MOON_PA_DE440_FRAME: Frame = MOON_PA_DE440_FRAME;
    #[classattr]
    const MOON_MOD_FRAME: Frame = MOON_MOD_FRAME;
    #[classattr]
    const MOON_TOD_FRAME: Frame = MOON_TOD_FRAME;
    #[classattr]
    const IAU_MOON_FRAME: Frame = IAU_MOON_FRAME;
    #[classattr]
    const IAU_MARS_FRAME: Frame = IAU_MARS_FRAME;
    #[classattr]
    const IAU_JUPITER_FRAME: Frame = IAU_JUPITER_FRAME;
    #[classattr]
    const IAU_SATURN_FRAME: Frame = IAU_SATURN_FRAME;
    #[classattr]
    const IAU_NEPTUNE_FRAME: Frame = IAU_NEPTUNE_FRAME;
    #[classattr]
    const IAU_URANUS_FRAME: Frame = IAU_URANUS_FRAME;
    #[classattr]
    const MERCURY_INERTIAL_FRAME: Frame = MERCURY_INERTIAL_FRAME;
    #[classattr]
    const VENUS_INERTIAL_FRAME: Frame = VENUS_INERTIAL_FRAME;
    #[classattr]
    const MARS_INERTIAL_FRAME: Frame = MARS_INERTIAL_FRAME;
    #[classattr]
    const JUPITER_INERTIAL_FRAME: Frame = JUPITER_INERTIAL_FRAME;
    #[classattr]
    const SATURN_INERTIAL_FRAME: Frame = SATURN_INERTIAL_FRAME;
    #[classattr]
    const URANUS_INERTIAL_FRAME: Frame = URANUS_INERTIAL_FRAME;
    #[classattr]
    const NEPTUNE_INERTIAL_FRAME: Frame = NEPTUNE_INERTIAL_FRAME;
}

#[pyclass]
#[pyo3(module = "anise.constants")]
struct Orientations {}

#[pymethods]
impl Orientations {
    #[classattr]
    const J2000: i32 = J2000;
    #[classattr]
    const ECLIPJ2000: i32 = ECLIPJ2000;
    #[classattr]
    const ICRS: i32 = ICRS;
    #[classattr]
    const EARTH_MOD: i32 = EARTH_MOD;
    #[classattr]
    const EARTH_MOD_2000: i32 = EARTH_MOD_2000;
    #[classattr]
    const EARTH_MOD_1976: i32 = EARTH_MOD_1976;
    #[classattr]
    const EARTH_TOD: i32 = EARTH_TOD;
    #[classattr]
    const EARTH_TOD_2000A: i32 = EARTH_TOD_2000A;
    #[classattr]
    const EARTH_TOD_2000B: i32 = EARTH_TOD_2000B;
    #[classattr]
    const EARTH_TOD_1980: i32 = EARTH_TOD_1980;
    #[classattr]
    const EARTH_TEME: i32 = EARTH_TEME;
    #[classattr]
    const EARTH_TEME_LEGACY: i32 = EARTH_TEME_LEGACY;
    #[classattr]
    const IAU_MERCURY: i32 = IAU_MERCURY;
    #[classattr]
    const IAU_VENUS: i32 = IAU_VENUS;
    #[classattr]
    const IAU_EARTH: i32 = IAU_EARTH;
    #[classattr]
    const IAU_MOON: i32 = IAU_MOON;
    #[classattr]
    const MOON_ME: i32 = MOON_ME;
    #[classattr]
    const MOON_ME_DE421: i32 = MOON_ME_DE421;
    #[classattr]
    const MOON_ME_DE440_ME421: i32 = MOON_ME_DE440_ME421;
    #[classattr]
    const MOON_MOD: i32 = MOON_MOD;
    #[classattr]
    const MOON_TOD: i32 = MOON_TOD;
    #[classattr]
    const MOON_PA: i32 = MOON_PA;
    #[classattr]
    const MOON_PA_DE421: i32 = MOON_PA_DE421;
    #[classattr]
    const MOON_PA_DE440: i32 = MOON_PA_DE440;
    #[classattr]
    const ITRF93: i32 = ITRF93;
    #[classattr]
    const IAU_MARS: i32 = IAU_MARS;
    #[classattr]
    const IAU_JUPITER: i32 = IAU_JUPITER;
    #[classattr]
    const IAU_SATURN: i32 = IAU_SATURN;
    #[classattr]
    const IAU_NEPTUNE: i32 = IAU_NEPTUNE;
    #[classattr]
    const IAU_URANUS: i32 = IAU_URANUS;
}

#[pyclass]
#[pyo3(module = "anise.constants")]
struct CelestialObjects {}

#[pymethods]
impl CelestialObjects {
    #[classattr]
    const SOLAR_SYSTEM_BARYCENTER: i32 = SOLAR_SYSTEM_BARYCENTER;
    #[classattr]
    const MERCURY: i32 = MERCURY;
    #[classattr]
    const VENUS: i32 = VENUS;
    #[classattr]
    const EARTH_MOON_BARYCENTER: i32 = EARTH_MOON_BARYCENTER;
    #[classattr]
    const MARS_BARYCENTER: i32 = MARS_BARYCENTER;
    #[classattr]
    const JUPITER_BARYCENTER: i32 = JUPITER_BARYCENTER;
    #[classattr]
    const SATURN_BARYCENTER: i32 = SATURN_BARYCENTER;
    #[classattr]
    const URANUS_BARYCENTER: i32 = URANUS_BARYCENTER;
    #[classattr]
    const NEPTUNE_BARYCENTER: i32 = NEPTUNE_BARYCENTER;
    #[classattr]
    const PLUTO_BARYCENTER: i32 = PLUTO_BARYCENTER;
    #[classattr]
    const SUN: i32 = SUN;
    #[classattr]
    const MOON: i32 = MOON;
    #[classattr]
    const EARTH: i32 = EARTH;
    #[classattr]
    const MARS: i32 = MARS;
    #[classattr]
    const JUPITER: i32 = JUPITER;
    #[classattr]
    const SATURN: i32 = SATURN;
    #[classattr]
    const URANUS: i32 = URANUS;
    #[classattr]
    const NEPTUNE: i32 = NEPTUNE;
}

#[pyclass]
#[pyo3(module = "anise.constants")]
struct UsualConstants {}

#[pymethods]
impl UsualConstants {
    #[classattr]
    /// Mean angular velocity of the Earth in deg/s
    /// Source: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016 (confirmed by https://hpiers.obspm.fr/eop-pc/models/constants.html)
    const MEAN_EARTH_ANGULAR_VELOCITY_DEG_S: f64 = MEAN_EARTH_ANGULAR_VELOCITY_DEG_S;
    #[classattr]
    /// Mean angular velocity of the Moon in deg/s, computed from hifitime:
    /// ```py
    /// >>> moon_period = Unit.Day*27+Unit.Hour*7+Unit.Minute*43+Unit.Second*12
    /// >>> tau/moon_period.to_seconds()
    /// 2.661698975163682e-06
    /// ```
    /// Source: https://www.britannica.com/science/month#ref225844 via https://en.wikipedia.org/w/index.php?title=Lunar_day&oldid=1180701337
    const MEAN_MOON_ANGULAR_VELOCITY_DEG_S: f64 = MEAN_MOON_ANGULAR_VELOCITY_DEG_S;
    #[classattr]
    /// Speed of light in kilometers per second (km/s)
    const SPEED_OF_LIGHT_KM_S: f64 = SPEED_OF_LIGHT_KM_S;
}

#[pymodule]
pub(crate) fn constants(_py: Python, sm: &Bound<'_, PyModule>) -> PyResult<()> {
    sm.add_class::<CelestialObjects>()?;
    sm.add_class::<Frames>()?;
    sm.add_class::<Orientations>()?;
    sm.add_class::<UsualConstants>()?;

    Ok(())
}
