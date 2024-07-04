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

use anise::astro::orbit::ECC_EPSILON;
use anise::astro::PhysicsResult;
use anise::prelude::{Frame, Orbit};
use snafu::ResultExt;

use super::AstroError;
use crate::cosmic::AstroPhysicsSnafu;
use crate::linalg::{Vector3, U7};
use crate::md::StateParameter;
use crate::time::Epoch;
use crate::TimeTagged;
use hyperdual::linalg::norm;
use hyperdual::{Float, OHyperdual};
use std::f64::consts::PI;
use std::fmt;

/// Orbit defines an orbital state
///
/// Unless noted otherwise, algorithms are from GMAT 2016a [StateConversionUtil.cpp](https://github.com/ChristopherRabotin/GMAT/blob/37201a6290e7f7b941bc98ee973a527a5857104b/src/base/util/StateConversionUtil.cpp).
/// Regardless of the constructor used, this struct stores all the state information in Cartesian coordinates
/// as these are always non singular.
/// _Note:_ although not yet supported, this struct may change once True of Date or other nutation frames
/// are added to the toolkit.
#[derive(Copy, Clone, Debug)]
pub struct OrbitDual {
    /// in km
    pub x: OHyperdual<f64, U7>,
    /// in km
    pub y: OHyperdual<f64, U7>,
    /// in km
    pub z: OHyperdual<f64, U7>,
    /// in km/s
    pub vx: OHyperdual<f64, U7>,
    /// in km/s
    pub vy: OHyperdual<f64, U7>,
    /// in km/s
    pub vz: OHyperdual<f64, U7>,
    pub dt: Epoch,
    /// Frame contains everything we need to compute state information
    pub frame: Frame,
}

impl From<Orbit> for OrbitDual {
    /// Initialize a new OrbitDual from an orbit, no other initializers
    fn from(orbit: Orbit) -> Self {
        Self {
            x: OHyperdual::from_slice(&[orbit.radius_km.x, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0]),
            y: OHyperdual::from_slice(&[orbit.radius_km.y, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0]),
            z: OHyperdual::from_slice(&[orbit.radius_km.z, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0]),
            vx: OHyperdual::from_slice(&[orbit.velocity_km_s.x, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
            vy: OHyperdual::from_slice(&[orbit.velocity_km_s.y, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0]),
            vz: OHyperdual::from_slice(&[orbit.velocity_km_s.z, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
            dt: orbit.epoch,
            frame: orbit.frame,
        }
    }
}

/// A type which stores the partial of an element
#[derive(Copy, Clone, Debug)]
pub struct OrbitPartial {
    pub param: StateParameter,
    pub dual: OHyperdual<f64, U7>,
}

impl OrbitPartial {
    /// Returns the real value of this parameter
    pub fn real(&self) -> f64 {
        self.dual[0]
    }
    /// The partial of this parameter with respect to X
    pub fn wtr_x(&self) -> f64 {
        self.dual[1]
    }
    /// The partial of this parameter with respect to Y
    pub fn wtr_y(&self) -> f64 {
        self.dual[2]
    }
    /// The partial of this parameter with respect to Z
    pub fn wtr_z(&self) -> f64 {
        self.dual[3]
    }
    /// The partial of this parameter with respect to VX
    pub fn wtr_vx(&self) -> f64 {
        self.dual[4]
    }
    /// The partial of this parameter with respect to VY
    pub fn wtr_vy(&self) -> f64 {
        self.dual[5]
    }
    /// The partial of this parameter with respect to VZ
    pub fn wtr_vz(&self) -> f64 {
        self.dual[6]
    }
}

impl fmt::Display for OrbitPartial {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{:?} {}", self.param, self.dual)
    }
}

impl OrbitDual {
    pub fn partial_for(&self, param: StateParameter) -> Result<OrbitPartial, AstroError> {
        match param {
            StateParameter::X => Ok(OrbitPartial {
                dual: self.x,
                param: StateParameter::X,
            }),
            StateParameter::Y => Ok(OrbitPartial {
                dual: self.y,
                param: StateParameter::Y,
            }),
            StateParameter::Z => Ok(OrbitPartial {
                dual: self.z,
                param: StateParameter::Z,
            }),
            StateParameter::VX => Ok(OrbitPartial {
                dual: self.vx,
                param: StateParameter::VX,
            }),
            StateParameter::VY => Ok(OrbitPartial {
                dual: self.vy,
                param: StateParameter::VY,
            }),
            StateParameter::VZ => Ok(OrbitPartial {
                dual: self.vz,
                param: StateParameter::VZ,
            }),
            StateParameter::Rmag => Ok(self.rmag_km()),
            StateParameter::Vmag => Ok(self.vmag_km_s()),
            StateParameter::HX => Ok(self.hx()),
            StateParameter::HY => Ok(self.hy()),
            StateParameter::HZ => Ok(self.hz()),
            StateParameter::Hmag => Ok(self.hmag()),
            StateParameter::Energy => Ok(self.energy_km2_s2().context(AstroPhysicsSnafu)?),
            StateParameter::SMA => Ok(self.sma_km().context(AstroPhysicsSnafu)?),
            StateParameter::Eccentricity => Ok(self.ecc().context(AstroPhysicsSnafu)?),
            StateParameter::Inclination => Ok(self.inc_deg()),
            StateParameter::AoP => Ok(self.aop_deg().context(AstroPhysicsSnafu)?),
            StateParameter::AoL => Ok(self.aol_deg().context(AstroPhysicsSnafu)?),
            StateParameter::RAAN => Ok(self.raan_deg()),
            StateParameter::Periapsis => Ok(self.periapsis_km().context(AstroPhysicsSnafu)?),
            StateParameter::Apoapsis => Ok(self.apoapsis_km().context(AstroPhysicsSnafu)?),
            StateParameter::TrueLongitude => Ok(self.tlong_deg().context(AstroPhysicsSnafu)?),
            StateParameter::FlightPathAngle => Ok(self.fpa_deg().context(AstroPhysicsSnafu)?),
            StateParameter::MeanAnomaly => Ok(self.ma_deg().context(AstroPhysicsSnafu)?),
            StateParameter::EccentricAnomaly => Ok(self.ea_deg().context(AstroPhysicsSnafu)?),
            StateParameter::Height => Ok(self.height_km().context(AstroPhysicsSnafu)?),
            StateParameter::Latitude => Ok(self.latitude_deg().context(AstroPhysicsSnafu)?),
            StateParameter::Longitude => Ok(self.longitude_deg()),
            StateParameter::C3 => Ok(self.c3().context(AstroPhysicsSnafu)?),
            StateParameter::RightAscension => Ok(self.right_ascension_deg()),
            StateParameter::Declination => Ok(self.declination_deg()),
            StateParameter::HyperbolicAnomaly => self.hyperbolic_anomaly_deg(),
            StateParameter::SemiParameter => {
                Ok(self.semi_parameter_km().context(AstroPhysicsSnafu)?)
            }
            StateParameter::SemiMinorAxis => {
                Ok(self.semi_minor_axis_km().context(AstroPhysicsSnafu)?)
            }
            StateParameter::TrueAnomaly => Ok(self.ta_deg().context(AstroPhysicsSnafu)?),
            _ => Err(AstroError::PartialsUndefined),
        }
    }

    /// Returns the magnitude of the radius vector in km
    pub fn rmag_km(&self) -> OrbitPartial {
        OrbitPartial {
            param: StateParameter::Rmag,
            dual: (self.x.powi(2) + self.y.powi(2) + self.z.powi(2)).sqrt(),
        }
    }

    /// Returns the magnitude of the velocity vector in km/s
    pub fn vmag_km_s(&self) -> OrbitPartial {
        OrbitPartial {
            param: StateParameter::Vmag,
            dual: (self.vx.powi(2) + self.vy.powi(2) + self.vz.powi(2)).sqrt(),
        }
    }

    /// Returns the radius vector of this Orbit in [km, km, km]
    pub(crate) fn radius(&self) -> Vector3<OHyperdual<f64, U7>> {
        Vector3::new(self.x, self.y, self.z)
    }

    /// Returns the velocity vector of this Orbit in [km/s, km/s, km/s]
    pub(crate) fn velocity(&self) -> Vector3<OHyperdual<f64, U7>> {
        Vector3::new(self.vx, self.vy, self.vz)
    }

    /// Returns the orbital momentum vector
    pub(crate) fn hvec(&self) -> Vector3<OHyperdual<f64, U7>> {
        self.radius().cross(&self.velocity())
    }

    /// Returns the orbital momentum value on the X axis
    pub fn hx(&self) -> OrbitPartial {
        OrbitPartial {
            dual: self.hvec()[0],
            param: StateParameter::HX,
        }
    }

    /// Returns the orbital momentum value on the Y axis
    pub fn hy(&self) -> OrbitPartial {
        OrbitPartial {
            dual: self.hvec()[1],
            param: StateParameter::HY,
        }
    }

    /// Returns the orbital momentum value on the Z axis
    pub fn hz(&self) -> OrbitPartial {
        OrbitPartial {
            dual: self.hvec()[2],
            param: StateParameter::HZ,
        }
    }

    /// Returns the norm of the orbital momentum
    pub fn hmag(&self) -> OrbitPartial {
        OrbitPartial {
            dual: norm(&self.hvec()),
            param: StateParameter::Hmag,
        }
    }

    /// Returns the specific mechanical energy
    pub fn energy_km2_s2(&self) -> PhysicsResult<OrbitPartial> {
        Ok(OrbitPartial {
            dual: self.vmag_km_s().dual.powi(2) / OHyperdual::from(2.0)
                - OHyperdual::from(self.frame.mu_km3_s2()?) / self.rmag_km().dual,
            param: StateParameter::Energy,
        })
    }

    /// Returns the semi-major axis in km
    pub fn sma_km(&self) -> PhysicsResult<OrbitPartial> {
        Ok(OrbitPartial {
            dual: -OHyperdual::from(self.frame.mu_km3_s2()?)
                / (OHyperdual::from(2.0) * self.energy_km2_s2()?.dual),
            param: StateParameter::SMA,
        })
    }

    /// Returns the eccentricity vector (no unit)
    pub(crate) fn evec(&self) -> PhysicsResult<Vector3<OHyperdual<f64, U7>>> {
        let r = self.radius();
        let v = self.velocity();
        let hgm = OHyperdual::from(self.frame.mu_km3_s2()?);
        // Split up this operation because it doesn't seem to be implemented
        // ((norm(&v).powi(2) - hgm / norm(&r)) * r - (r.dot(&v)) * v) / hgm
        Ok(Vector3::new(
            ((norm(&v).powi(2) - hgm / norm(&r)) * r[0] - (r.dot(&v)) * v[0]) / hgm,
            ((norm(&v).powi(2) - hgm / norm(&r)) * r[1] - (r.dot(&v)) * v[1]) / hgm,
            ((norm(&v).powi(2) - hgm / norm(&r)) * r[2] - (r.dot(&v)) * v[2]) / hgm,
        ))
    }

    /// Returns the eccentricity (no unit)
    pub fn ecc(&self) -> PhysicsResult<OrbitPartial> {
        Ok(OrbitPartial {
            dual: norm(&self.evec()?),
            param: StateParameter::Eccentricity,
        })
    }

    /// Returns the inclination in degrees
    pub fn inc_deg(&self) -> OrbitPartial {
        OrbitPartial {
            dual: (self.hvec()[(2, 0)] / self.hmag().dual).acos().to_degrees(),
            param: StateParameter::Inclination,
        }
    }

    /// Returns the argument of periapsis in degrees
    pub fn aop_deg(&self) -> PhysicsResult<OrbitPartial> {
        let n = Vector3::new(
            OHyperdual::from(0.0),
            OHyperdual::from(0.0),
            OHyperdual::from(1.0),
        )
        .cross(&self.hvec());
        let aop = (n.dot(&self.evec()?) / (norm(&n) * self.ecc()?.dual)).acos();
        if aop.is_nan() {
            error!("AoP is NaN");
            Ok(OrbitPartial {
                dual: OHyperdual::from(0.0),
                param: StateParameter::AoP,
            })
        } else if self.evec()?[2].real() < 0.0 {
            Ok(OrbitPartial {
                dual: (OHyperdual::from(2.0 * PI) - aop).to_degrees(),
                param: StateParameter::AoP,
            })
        } else {
            Ok(OrbitPartial {
                dual: aop.to_degrees(),
                param: StateParameter::AoP,
            })
        }
    }

    /// Returns the right ascension of ther ascending node in degrees
    pub fn raan_deg(&self) -> OrbitPartial {
        let n = Vector3::new(
            OHyperdual::from(0.0),
            OHyperdual::from(0.0),
            OHyperdual::from(1.0),
        )
        .cross(&self.hvec());
        let raan = (n[(0, 0)] / norm(&n)).acos();
        if raan.is_nan() {
            warn!("RAAN is NaN");
            OrbitPartial {
                dual: OHyperdual::from(0.0),
                param: StateParameter::RAAN,
            }
        } else if n[(1, 0)] < 0.0 {
            OrbitPartial {
                dual: (OHyperdual::from(2.0 * PI) - raan).to_degrees(),
                param: StateParameter::RAAN,
            }
        } else {
            OrbitPartial {
                dual: raan.to_degrees(),
                param: StateParameter::RAAN,
            }
        }
    }

    /// Returns the true anomaly in degrees between 0 and 360.0
    pub fn ta_deg(&self) -> PhysicsResult<OrbitPartial> {
        if self.ecc()?.real() < ECC_EPSILON {
            warn!(
                "true anomaly ill-defined for circular orbit (e = {})",
                self.ecc()?
            );
        }
        let cos_nu = self.evec()?.dot(&self.radius()) / (self.ecc()?.dual * self.rmag_km().dual);
        if (cos_nu.real().abs() - 1.0).abs() < f64::EPSILON {
            // This bug drove me crazy when writing SMD in Go in 2017.
            if cos_nu > 1.0 {
                Ok(OrbitPartial {
                    dual: OHyperdual::from(180.0),
                    param: StateParameter::TrueAnomaly,
                })
            } else {
                Ok(OrbitPartial {
                    dual: OHyperdual::from(0.0),
                    param: StateParameter::TrueAnomaly,
                })
            }
        } else {
            let ta = cos_nu.acos();
            if ta.is_nan() {
                warn!("TA is NaN");
                Ok(OrbitPartial {
                    dual: OHyperdual::from(0.0),
                    param: StateParameter::TrueAnomaly,
                })
            } else if self.radius().dot(&self.velocity()) < 0.0 {
                Ok(OrbitPartial {
                    dual: (OHyperdual::from(2.0 * PI) - ta).to_degrees(),
                    param: StateParameter::TrueAnomaly,
                })
            } else {
                Ok(OrbitPartial {
                    dual: ta.to_degrees(),
                    param: StateParameter::TrueAnomaly,
                })
            }
        }
    }

    /// Returns the true longitude in degrees
    pub fn tlong_deg(&self) -> PhysicsResult<OrbitPartial> {
        // Angles already in degrees
        Ok(OrbitPartial {
            dual: self.aop_deg()?.dual + self.raan_deg().dual + self.ta_deg()?.dual,
            param: StateParameter::TrueLongitude,
        })
    }

    /// Returns the argument of latitude in degrees
    ///
    /// NOTE: If the orbit is near circular, the AoL will be computed from the true longitude
    /// instead of relying on the ill-defined true anomaly.
    pub fn aol_deg(&self) -> PhysicsResult<OrbitPartial> {
        if self.ecc()?.real() < ECC_EPSILON {
            Ok(OrbitPartial {
                dual: self.tlong_deg()?.dual - self.raan_deg().dual,
                param: StateParameter::AoL,
            })
        } else {
            Ok(OrbitPartial {
                dual: self.aop_deg()?.dual + self.ta_deg()?.dual,
                param: StateParameter::AoL,
            })
        }
    }

    /// Returns the radius of periapsis (or perigee around Earth), in kilometers.
    pub fn periapsis_km(&self) -> PhysicsResult<OrbitPartial> {
        Ok(OrbitPartial {
            dual: self.sma_km()?.dual * (OHyperdual::from(1.0) - self.ecc()?.dual),
            param: StateParameter::Periapsis,
        })
    }

    /// Returns the radius of apoapsis (or apogee around Earth), in kilometers.
    pub fn apoapsis_km(&self) -> PhysicsResult<OrbitPartial> {
        Ok(OrbitPartial {
            dual: self.sma_km()?.dual * (OHyperdual::from(1.0) + self.ecc()?.dual),
            param: StateParameter::Apoapsis,
        })
    }

    /// Returns the eccentric anomaly in degrees
    ///
    /// This is a conversion from GMAT's StateConversionUtil::TrueToEccentricAnomaly
    pub fn ea_deg(&self) -> PhysicsResult<OrbitPartial> {
        let (sin_ta, cos_ta) = self.ta_deg()?.dual.to_radians().sin_cos();
        let ecc_cos_ta = self.ecc()?.dual * cos_ta;
        let sin_ea = ((OHyperdual::from(1.0) - self.ecc()?.dual.powi(2)).sqrt() * sin_ta)
            / (OHyperdual::from(1.0) + ecc_cos_ta);
        let cos_ea = (self.ecc()?.dual + cos_ta) / (OHyperdual::from(1.0) + ecc_cos_ta);
        // The atan2 function is a bit confusing: https://doc.rust-lang.org/std/primitive.f64.html#method.atan2 .
        Ok(OrbitPartial {
            dual: sin_ea.atan2(cos_ea).to_degrees(),
            param: StateParameter::EccentricAnomaly,
        })
    }

    /// Returns the flight path angle in degrees
    pub fn fpa_deg(&self) -> PhysicsResult<OrbitPartial> {
        let nu = self.ta_deg()?.dual.to_radians();
        let ecc = self.ecc()?.dual;
        let denom =
            (OHyperdual::from(1.0) + OHyperdual::from(2.0) * ecc * nu.cos() + ecc.powi(2)).sqrt();
        let sin_fpa = ecc * nu.sin() / denom;
        let cos_fpa = OHyperdual::from(1.0) + ecc * nu.cos() / denom;
        Ok(OrbitPartial {
            dual: sin_fpa.atan2(cos_fpa).to_degrees(),
            param: StateParameter::FlightPathAngle,
        })
    }

    /// Returns the mean anomaly in degrees
    pub fn ma_deg(&self) -> PhysicsResult<OrbitPartial> {
        if self.ecc()?.real() < 1.0 {
            Ok(OrbitPartial {
                dual: (self.ea_deg()?.dual.to_radians()
                    - self.ecc()?.dual * self.ea_deg()?.dual.to_radians().sin())
                .to_degrees(),
                param: StateParameter::MeanAnomaly,
            })
        } else if self.ecc()?.real() > 1.0 {
            info!("computing the hyperbolic anomaly");
            // From GMAT's TrueToHyperbolicAnomaly
            Ok(OrbitPartial {
                dual: ((self.ta_deg()?.dual.to_radians().sin()
                    * (self.ecc()?.dual.powi(2) - OHyperdual::from(1.0)))
                .sqrt()
                    / (OHyperdual::from(1.0)
                        + self.ecc()?.dual * self.ta_deg()?.dual.to_radians().cos()))
                .asinh()
                .to_degrees(),
                param: StateParameter::MeanAnomaly,
            })
        } else {
            error!("parabolic orbit: setting mean anomaly to 0.0");
            Ok(OrbitPartial {
                dual: OHyperdual::from(0.0),
                param: StateParameter::MeanAnomaly,
            })
        }
    }

    /// Returns the semi parameter (or semilatus rectum)
    pub fn semi_parameter_km(&self) -> PhysicsResult<OrbitPartial> {
        Ok(OrbitPartial {
            dual: self.sma_km()?.dual * (OHyperdual::from(1.0) - self.ecc()?.dual.powi(2)),
            param: StateParameter::SemiParameter,
        })
    }

    /// Returns the geodetic longitude (λ) in degrees. Value is between 0 and 360 degrees.
    ///
    /// Although the reference is not Vallado, the math from Vallado proves to be equivalent.
    /// Reference: G. Xu and Y. Xu, "GPS", DOI 10.1007/978-3-662-50367-6_2, 2016
    pub fn longitude_deg(&self) -> OrbitPartial {
        OrbitPartial {
            dual: self.y.atan2(self.x).to_degrees(),
            param: StateParameter::Longitude,
        }
    }

    /// Returns the geodetic latitude (φ) in degrees. Value is between -180 and +180 degrees.
    ///
    /// Reference: Vallado, 4th Ed., Algorithm 12 page 172.
    pub fn latitude_deg(&self) -> PhysicsResult<OrbitPartial> {
        let flattening = self.frame.flattening()?;
        let eps = 1e-12;
        let max_attempts = 20;
        let mut attempt_no = 0;
        let r_delta = (self.x.powi(2) + self.y.powi(2)).sqrt();
        let mut latitude = (self.z / self.rmag_km().dual).asin();
        let e2 = OHyperdual::from(flattening * (2.0 - flattening));
        loop {
            attempt_no += 1;
            let c_earth = OHyperdual::from(self.frame.semi_major_radius_km()?)
                / ((OHyperdual::from(1.0) - e2 * (latitude).sin().powi(2)).sqrt());
            let new_latitude = (self.z + c_earth * e2 * (latitude).sin()).atan2(r_delta);
            if (latitude - new_latitude).abs() < eps {
                return Ok(OrbitPartial {
                    dual: new_latitude.to_degrees(),
                    param: StateParameter::Latitude,
                });
            } else if attempt_no >= max_attempts {
                error!(
                    "geodetic latitude failed to converge -- error = {}",
                    (latitude - new_latitude).abs()
                );
                return Ok(OrbitPartial {
                    dual: new_latitude.to_degrees(),
                    param: StateParameter::Latitude,
                });
            }
            latitude = new_latitude;
        }
    }

    /// Returns the geodetic height in km.
    ///
    /// Reference: Vallado, 4th Ed., Algorithm 12 page 172.
    pub fn height_km(&self) -> PhysicsResult<OrbitPartial> {
        let flattening = self.frame.flattening()?;
        let semi_major_radius = self.frame.semi_major_radius_km()?;
        let e2 = OHyperdual::from(flattening * (2.0 - flattening));
        let latitude = self.latitude_deg()?.dual.to_radians();
        let sin_lat = latitude.sin();
        if (latitude - OHyperdual::from(1.0)).abs() < 0.1 {
            // We are near poles, let's use another formulation.
            let s_earth0: f64 = semi_major_radius * (1.0 - flattening).powi(2);
            let s_earth = OHyperdual::from(s_earth0)
                / ((OHyperdual::from(1.0) - e2 * sin_lat.powi(2)).sqrt());
            Ok(OrbitPartial {
                dual: self.z / latitude.sin() - s_earth,
                param: StateParameter::Height,
            })
        } else {
            let c_earth = OHyperdual::from(semi_major_radius)
                / ((OHyperdual::from(1.0) - e2 * sin_lat.powi(2)).sqrt());
            let r_delta = (self.x.powi(2) + self.y.powi(2)).sqrt();
            Ok(OrbitPartial {
                dual: r_delta / latitude.cos() - c_earth,
                param: StateParameter::Height,
            })
        }
    }

    /// Returns the right ascension of this orbit in degrees
    #[allow(clippy::eq_op)]
    pub fn right_ascension_deg(&self) -> OrbitPartial {
        OrbitPartial {
            dual: (self.y.atan2(self.x)).to_degrees(),
            param: StateParameter::RightAscension,
        }
    }

    /// Returns the declination of this orbit in degrees
    #[allow(clippy::eq_op)]
    pub fn declination_deg(&self) -> OrbitPartial {
        OrbitPartial {
            dual: (self.z / self.rmag_km().dual).asin().to_degrees(),
            param: StateParameter::Declination,
        }
    }

    /// Returns the semi minor axis in km, includes code for a hyperbolic orbit
    pub fn semi_minor_axis_km(&self) -> PhysicsResult<OrbitPartial> {
        if self.ecc()?.real() <= 1.0 {
            Ok(OrbitPartial {
                dual: ((self.sma_km()?.dual * self.ecc()?.dual).powi(2)
                    - self.sma_km()?.dual.powi(2))
                .sqrt(),
                param: StateParameter::SemiMinorAxis,
            })
        } else {
            Ok(OrbitPartial {
                dual: self.hmag().dual.powi(2)
                    / (OHyperdual::from(self.frame.mu_km3_s2()?)
                        * (self.ecc()?.dual.powi(2) - OHyperdual::from(1.0)).sqrt()),
                param: StateParameter::SemiMinorAxis,
            })
        }
    }

    /// Returns the velocity declination of this orbit in degrees
    pub fn velocity_declination(&self) -> OrbitPartial {
        OrbitPartial {
            dual: (self.vz / self.vmag_km_s().dual).asin().to_degrees(),
            param: StateParameter::VelocityDeclination,
        }
    }

    /// Returns the $C_3$ of this orbit
    pub fn c3(&self) -> PhysicsResult<OrbitPartial> {
        Ok(OrbitPartial {
            dual: -OHyperdual::from(self.frame.mu_km3_s2()?) / self.sma_km()?.dual,
            param: StateParameter::C3,
        })
    }

    /// Returns the hyperbolic anomaly in degrees between 0 and 360.0
    pub fn hyperbolic_anomaly_deg(&self) -> Result<OrbitPartial, AstroError> {
        let ecc = self.ecc().context(AstroPhysicsSnafu)?;
        if ecc.real() <= 1.0 {
            Err(AstroError::NotHyperbolic)
        } else {
            let (sin_ta, cos_ta) = self
                .ta_deg()
                .context(AstroPhysicsSnafu)?
                .dual
                .to_radians()
                .sin_cos();
            let sinh_h = (sin_ta * (ecc.dual.powi(2) - OHyperdual::from(1.0)).sqrt())
                / (OHyperdual::from(1.0) + ecc.dual * cos_ta);
            Ok(OrbitPartial {
                dual: sinh_h.asinh().to_degrees(),
                param: StateParameter::HyperbolicAnomaly,
            })
        }
    }
}

impl TimeTagged for OrbitDual {
    fn epoch(&self) -> Epoch {
        self.dt
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.dt = epoch
    }
}
