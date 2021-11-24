/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use super::serde_derive::Deserialize;
use super::EpochFormat;
use crate::cosmic::{Cosm, Frame, Orbit};
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::StateParameter;
use crate::od::estimate::NavSolution;
use crate::{NyxError, State};
use std::cmp::PartialEq;
use std::collections::HashMap;
use std::fmt;
use std::str::FromStr;
use std::sync::Arc;

#[derive(Deserialize)]
pub struct OutputSerde {
    pub filename: String,
    /// If not specified, the standard
    pub headers: Option<Vec<String>>,
}

impl OutputSerde {
    pub fn to_state_formatter(&self, cosm: Arc<Cosm>) -> Result<StateFormatter, NyxError> {
        match &self.headers {
            Some(hdr) => StateFormatter::from_headers(
                hdr.iter().map(|s| s.as_str()).collect::<Vec<&str>>(),
                self.filename.clone(),
                cosm,
            ),
            None => Ok(StateFormatter::default(self.filename.clone(), cosm)),
        }
    }

    pub fn to_nav_sol_formatter(&self, cosm: Arc<Cosm>) -> Result<NavSolutionFormatter, NyxError> {
        match &self.headers {
            Some(hdr) => {
                NavSolutionFormatter::from_headers(hdr.to_vec(), self.filename.clone(), cosm)
            }
            None => Ok(NavSolutionFormatter::default(self.filename.clone(), cosm)),
        }
    }
}

/// Allowed headers, with an optional frame.
/// TODO: Support units
#[allow(non_camel_case_types)]
#[derive(Clone, Debug, PartialEq)]
pub struct StateHeader {
    /// Stores either the state paramater or the epoch
    pub param: StateParameter,
    pub frame_name: Option<String>,
    pub epoch_fmt: Option<EpochFormat>,
}

impl From<StateParameter> for StateHeader {
    fn from(param: StateParameter) -> Self {
        StateHeader {
            param,
            frame_name: None,
            epoch_fmt: if param == StateParameter::Epoch {
                Some(EpochFormat::GregorianUtc)
            } else {
                None
            },
        }
    }
}

impl fmt::Display for StateHeader {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, fh: &mut fmt::Formatter) -> fmt::Result {
        let fmtd = match self.param {
            StateParameter::X
            | StateParameter::Y
            | StateParameter::Z
            | StateParameter::ApoapsisRadius
            | StateParameter::PeriapsisRadius
            | StateParameter::GeodeticHeight
            | StateParameter::SemiMinorAxis
            | StateParameter::SemiParameter
            | StateParameter::SMA
            | StateParameter::Rmag => {
                format!("{:?} (km)", self.param)
            }
            StateParameter::VX | StateParameter::VY | StateParameter::VZ | StateParameter::Vmag => {
                format!("{:?} (km/s)", self.param)
            }
            StateParameter::AoL
            | StateParameter::AoP
            | StateParameter::Declination
            | StateParameter::EccentricAnomaly
            | StateParameter::GeodeticLatitude
            | StateParameter::GeodeticLongitude
            | StateParameter::Inclination
            | StateParameter::MeanAnomaly
            | StateParameter::RightAscension
            | StateParameter::RAAN
            | StateParameter::TrueAnomaly
            | StateParameter::TrueLongitude => {
                format!("{:?} (deg)", self.param)
            }
            _ => format!("{:?}", self.param),
        };
        write!(fh, "{}", fmtd)?;
        if let Some(frame) = &self.frame_name {
            write!(fh, ":{}", frame)?;
        } else if let Some(epoch_fmt) = self.epoch_fmt {
            write!(fh, ":{:?}", epoch_fmt)?;
        }
        Ok(())
    }
}

impl Serialize for StateHeader {
    /// NOTE: This is not part of unit testing because there is no deseralization of Orbit (yet)
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        serializer.serialize_str(&format!("{}", self))
    }
}

/// Allowed headers, with an optional frame.
/// TODO: Support units
#[allow(non_camel_case_types)]
#[derive(Clone, Debug, PartialEq)]
pub enum NavSolutionHeader {
    /// The epoch in the specified format
    Epoch(EpochFormat),
    /// Headers for the estimated state
    EstimatedState(Vec<StateHeader>),
    /// Headers for the nominal state
    NominalState(Vec<StateHeader>),
    /// Orbit deviation X (km)
    Delta_x,
    /// Orbit deviation Y (km)
    Delta_y,
    /// Orbit deviation Z (km)
    Delta_z,
    /// Orbit deviation VX (km/s)
    Delta_vx,
    /// Orbit deviation VY (km/s)
    Delta_vy,
    /// Orbit deviation VZ (km/s)
    Delta_vz,
    /// Covariance matrix [1,1]
    Cx_x { frame: Option<String> },
    /// Covariance matrix [2,1]
    Cy_x { frame: Option<String> },
    /// Covariance matrix [2,2]
    Cy_y { frame: Option<String> },
    /// Covariance matrix [3,1]
    Cz_x { frame: Option<String> },
    /// Covariance matrix [3,2]
    Cz_y { frame: Option<String> },
    /// Covariance matrix [3,3]
    Cz_z { frame: Option<String> },
    /// Covariance matrix [4,1]
    Cx_dot_x { frame: Option<String> },
    /// Covariance matrix [4,2]
    Cx_dot_y { frame: Option<String> },
    /// Covariance matrix [4,3]
    Cx_dot_z { frame: Option<String> },
    /// Covariance matrix [4,4]
    Cx_dot_x_dot { frame: Option<String> },
    /// Covariance matrix [5,1]
    Cy_dot_x { frame: Option<String> },
    /// Covariance matrix [5,2]
    Cy_dot_y { frame: Option<String> },
    /// Covariance matrix [5,3]
    Cy_dot_z { frame: Option<String> },
    /// Covariance matrix [5,4]
    Cy_dot_x_dot { frame: Option<String> },
    /// Covariance matrix [5,5]
    Cy_dot_y_dot { frame: Option<String> },
    /// Covariance matrix [6,1]
    Cz_dot_x { frame: Option<String> },
    /// Covariance matrix [6,2]
    Cz_dot_y { frame: Option<String> },
    /// Covariance matrix [6,3]
    Cz_dot_z { frame: Option<String> },
    /// Covariance matrix [6,4]
    Cz_dot_x_dot { frame: Option<String> },
    /// Covariance matrix [6,5]
    Cz_dot_y_dot { frame: Option<String> },
    /// Covariance matrix [6,6]
    Cz_dot_z_dot { frame: Option<String> },
    /// Boolean specifying whether this is a prediction or not
    Prediction,
    /// Norm of the position items of the covariance (Cx_x, Cy_y, Cz_z)
    Covar_pos,
    /// Norm of the velocity items of the covartiance (Cx_dot_x_dot, Cy_dot_y_dot, Cz_dot_z_dot)
    Covar_vel,
}

impl fmt::Display for NavSolutionHeader {
    fn fmt(&self, fh: &mut fmt::Formatter) -> fmt::Result {
        match self {
            NavSolutionHeader::Epoch(efmt) => write!(fh, "Epoch:{:?}", efmt),
            NavSolutionHeader::EstimatedState(hdr) => {
                let mut seq = Vec::with_capacity(hdr.len());
                for element in hdr {
                    seq.push(format!("Estimate:{}", element));
                }
                write!(fh, "{}", seq.join(","))
            }
            NavSolutionHeader::NominalState(hdr) => {
                let mut seq = Vec::with_capacity(hdr.len());
                for element in hdr {
                    seq.push(format!("Nominal:{}", element));
                }
                write!(fh, "{}", seq.join(","))
            }
            NavSolutionHeader::Delta_x => write!(fh, "delta_x"),
            NavSolutionHeader::Delta_y => write!(fh, "delta_y"),
            NavSolutionHeader::Delta_z => write!(fh, "delta_z"),
            NavSolutionHeader::Delta_vx => write!(fh, "delta_vx"),
            NavSolutionHeader::Delta_vy => write!(fh, "delta_vy"),
            NavSolutionHeader::Delta_vz => write!(fh, "delta_vz"),
            NavSolutionHeader::Cx_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cx_x:{}", f)
                } else {
                    write!(fh, "cx_x")
                }
            }
            NavSolutionHeader::Cy_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_x:{}", f)
                } else {
                    write!(fh, "cy_x")
                }
            }
            NavSolutionHeader::Cy_y { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_y:{}", f)
                } else {
                    write!(fh, "cy_y")
                }
            }
            NavSolutionHeader::Cz_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_x:{}", f)
                } else {
                    write!(fh, "cz_x")
                }
            }
            NavSolutionHeader::Cz_y { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_y:{}", f)
                } else {
                    write!(fh, "cz_y")
                }
            }
            NavSolutionHeader::Cz_z { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_z:{}", f)
                } else {
                    write!(fh, "cz_z")
                }
            }
            NavSolutionHeader::Cx_dot_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cx_dot_x:{}", f)
                } else {
                    write!(fh, "cx_dot_x")
                }
            }
            NavSolutionHeader::Cx_dot_y { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cx_dot_y:{}", f)
                } else {
                    write!(fh, "cx_dot_y")
                }
            }
            NavSolutionHeader::Cx_dot_z { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cx_dot_z:{}", f)
                } else {
                    write!(fh, "cx_dot_z")
                }
            }
            NavSolutionHeader::Cx_dot_x_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cx_dot_x_dot:{}", f)
                } else {
                    write!(fh, "cx_dot_x_dot")
                }
            }
            NavSolutionHeader::Cy_dot_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_dot_x:{}", f)
                } else {
                    write!(fh, "cy_dot_x")
                }
            }
            NavSolutionHeader::Cy_dot_y { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_dot_y:{}", f)
                } else {
                    write!(fh, "cy_dot_y")
                }
            }
            NavSolutionHeader::Cy_dot_z { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_dot_z:{}", f)
                } else {
                    write!(fh, "cy_dot_z")
                }
            }
            NavSolutionHeader::Cy_dot_x_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_dot_x_dot:{}", f)
                } else {
                    write!(fh, "cy_dot_x_dot")
                }
            }
            NavSolutionHeader::Cy_dot_y_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_dot_y_dot:{}", f)
                } else {
                    write!(fh, "cy_dot_y_dot")
                }
            }
            NavSolutionHeader::Cz_dot_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_x:{}", f)
                } else {
                    write!(fh, "cz_dot_x")
                }
            }
            NavSolutionHeader::Cz_dot_y { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_y:{}", f)
                } else {
                    write!(fh, "cz_dot_y")
                }
            }
            NavSolutionHeader::Cz_dot_z { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_z:{}", f)
                } else {
                    write!(fh, "cz_dot_z")
                }
            }
            NavSolutionHeader::Cz_dot_x_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_x_dot:{}", f)
                } else {
                    write!(fh, "cz_dot_x_dot")
                }
            }
            NavSolutionHeader::Cz_dot_y_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_y_dot:{}", f)
                } else {
                    write!(fh, "cz_dot_y_dot")
                }
            }
            NavSolutionHeader::Cz_dot_z_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_z_dot:{}", f)
                } else {
                    write!(fh, "cz_dot_z_dot")
                }
            }
            NavSolutionHeader::Prediction => write!(fh, "prediction"),
            NavSolutionHeader::Covar_pos => write!(fh, "covar_position"),
            NavSolutionHeader::Covar_vel => write!(fh, "covar_velocity"),
        }
    }
}

impl Serialize for NavSolutionHeader {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self {
            NavSolutionHeader::EstimatedState(hdr) => {
                let mut seq = serializer.serialize_seq(Some(hdr.len()))?;
                for element in hdr {
                    seq.serialize_element(&format!("Estimate:{}", element))?;
                }
                seq.end()
            }
            NavSolutionHeader::NominalState(hdr) => {
                let mut seq = serializer.serialize_seq(Some(hdr.len()))?;
                for element in hdr {
                    seq.serialize_element(&format!("Nominal:{}", element))?;
                }
                seq.end()
            }
            _ => serializer.serialize_str(&format!("{}", self)),
        }
    }
}

/// A formatter for states
#[derive(Clone)]
pub struct StateFormatter {
    pub filename: String,
    pub headers: Vec<StateHeader>,
    frames: HashMap<String, Frame>,
    cosm: Arc<Cosm>,
}

impl StateFormatter {
    /// ```
    /// extern crate nyx_space as nyx;
    /// use nyx::io::formatter::StateFormatter;
    /// use nyx::cosmic::Cosm;
    ///
    /// let cosm = Cosm::de438();
    /// // In this case, we're initializing the formatter to output the AoL and the eccentric anomaly in the EME2000 frame.
    /// let hdrs = vec!["AoL", "ea:eme2000"];
    /// StateFormatter::from_headers(hdrs, "nope".to_string(), cosm);
    /// ```
    pub fn from_headers(
        headers: Vec<&str>,
        filename: String,
        cosm: Arc<Cosm>,
    ) -> Result<Self, NyxError> {
        let mut frames = HashMap::new();
        let mut hdrs = Vec::with_capacity(20);
        // Rebuild the header tokens
        for hdr in &headers {
            let splt: Vec<&str> = hdr.split(':').collect();

            match splt[0].to_lowercase().as_str() {
                "epoch" => {
                    let epoch_fmt = if splt.len() == 2 {
                        match EpochFormat::from_str(splt[1]) {
                            Ok(e) => e,
                            Err(_) => {
                                return Err(NyxError::LoadingError(format!(
                                    "Unknown epoch format {}",
                                    splt[1]
                                )))
                            }
                        }
                    } else {
                        EpochFormat::GregorianUtc
                    };

                    let hdr = StateHeader {
                        param: StateParameter::Epoch,
                        frame_name: None,
                        epoch_fmt: Some(epoch_fmt),
                    };

                    hdrs.push(hdr);
                }
                _ => {
                    let frame_name = if splt.len() == 2 {
                        Some(splt[1].to_owned())
                    } else {
                        None
                    };

                    let param = StateParameter::from_str(splt[0].to_lowercase().as_str())?;

                    let hdr = StateHeader {
                        param,
                        frame_name,
                        epoch_fmt: None,
                    };

                    hdrs.push(hdr);
                }
            }

            if splt.len() == 2 && splt[0].to_lowercase() != "epoch" {
                // Get the frame
                match splt[1].to_lowercase().as_str() {
                    "ric" => frames.insert(splt[1].to_string(), Frame::RIC),
                    "rcn" => frames.insert(splt[1].to_string(), Frame::RCN),
                    "vnc" => frames.insert(splt[1].to_string(), Frame::VNC),
                    _ => frames.insert(splt[1].to_string(), cosm.try_frame(splt[1])?),
                };
            }
        }

        Ok(Self {
            filename,
            headers: hdrs,
            frames,
            cosm,
        })
    }

    /// Default headers are [Epoch (GregorianTai), X, Y, Z, VX, VY, VZ], where position is in km and velocity in km/s.
    pub fn default(filename: String, cosm: Arc<Cosm>) -> Self {
        Self {
            filename,
            headers: vec![
                From::from(StateParameter::Epoch),
                From::from(StateParameter::X),
                From::from(StateParameter::Y),
                From::from(StateParameter::Z),
                From::from(StateParameter::VX),
                From::from(StateParameter::VY),
                From::from(StateParameter::VZ),
            ],
            frames: HashMap::new(),
            cosm,
        }
    }

    pub fn fmt(&self, orbit: &Orbit) -> Vec<String> {
        // Start by computing the state in all of the frames needed
        let mut mapped = HashMap::new();
        for (name, frame) in &self.frames {
            let mapped_orbit = match frame {
                Frame::RCN | Frame::RIC | Frame::VNC => {
                    orbit.with_position_rotated_by(orbit.dcm_from_traj_frame(*frame).unwrap())
                }
                _ => self.cosm.frame_chg(orbit, *frame),
            };
            mapped.insert(name.to_lowercase(), mapped_orbit);
        }
        let mut formatted = Vec::with_capacity(self.headers.len());

        for hdr in &self.headers {
            // Grab the state in the other frame if needed
            let orbit = if hdr.frame_name.is_some() {
                &mapped[&hdr.frame_name.as_ref().unwrap().to_lowercase()]
            } else {
                orbit
            };

            formatted.push(match hdr.param {
                StateParameter::Epoch => hdr.epoch_fmt.as_ref().unwrap().format(orbit.dt),
                StateParameter::AoL => format!("{:.16}", orbit.aol()),
                StateParameter::AoP => format!("{:.16}", orbit.aop()),
                StateParameter::Apoapsis => format!("{:.16}", orbit.ta()),
                StateParameter::C3 => format!("{:.16}", orbit.c3()),
                StateParameter::Declination => format!("{:.16}", orbit.declination()),
                StateParameter::ApoapsisRadius => format!("{:.16}", orbit.apoapsis()),
                StateParameter::EccentricAnomaly => format!("{:.16}", orbit.ea()),
                StateParameter::Eccentricity => format!("{:.16}", orbit.ecc()),
                StateParameter::Energy => format!("{:.16}", orbit.energy()),
                StateParameter::GeodeticHeight => format!("{:.16}", orbit.geodetic_height()),
                StateParameter::GeodeticLatitude => format!("{:.16}", orbit.geodetic_latitude()),
                StateParameter::GeodeticLongitude => format!("{:.16}", orbit.geodetic_longitude()),
                StateParameter::FlightPathAngle => format!("{:.16}", orbit.fpa()),
                StateParameter::Hmag => format!("{:.16}", orbit.hmag()),
                StateParameter::HX => format!("{:.16}", orbit.hx()),
                StateParameter::HY => format!("{:.16}", orbit.hy()),
                StateParameter::HZ => format!("{:.16}", orbit.hz()),
                StateParameter::HyperbolicAnomaly => {
                    format!("{:.16}", orbit.hyperbolic_anomaly().unwrap())
                }
                StateParameter::Inclination => format!("{:.16}", orbit.inc()),
                StateParameter::MeanAnomaly => format!("{:.16}", orbit.ma()),
                StateParameter::Periapsis => format!("{:.16}", orbit.ta()),
                StateParameter::PeriapsisRadius => format!("{:.16}", orbit.periapsis()),
                StateParameter::Period => format!("{:.16}", orbit.period().in_seconds()),
                StateParameter::RightAscension => format!("{:.16}", orbit.right_ascension()),
                StateParameter::RAAN => format!("{:.16}", orbit.raan()),
                StateParameter::Rmag => format!("{:.16}", orbit.rmag()),
                StateParameter::SemiParameter => format!("{:.16}", orbit.semi_parameter()),
                StateParameter::SemiMinorAxis => format!("{:.16}", orbit.semi_minor_axis()),
                StateParameter::SMA => format!("{:.16}", orbit.sma()),
                StateParameter::TrueAnomaly => format!("{:.16}", orbit.ta()),
                StateParameter::TrueLongitude => format!("{:.16}", orbit.tlong()),
                StateParameter::VelocityDeclination => {
                    format!("{:.16}", orbit.velocity_declination())
                }
                StateParameter::Vmag => format!("{:.16}", orbit.vmag()),
                StateParameter::X => format!("{:.16}", orbit.x),
                StateParameter::Y => format!("{:.16}", orbit.y),
                StateParameter::Z => format!("{:.16}", orbit.z),
                StateParameter::VX => format!("{:.16}", orbit.vx),
                StateParameter::VY => format!("{:.16}", orbit.vy),
                StateParameter::VZ => format!("{:.16}", orbit.vz),
                StateParameter::FuelMass => {
                    unimplemented!("No fuel for an orbit, only for spacecraft!")
                }
                StateParameter::Custom { .. } => {
                    unimplemented!("Cannot format custom state parameters yet")
                }
                _ => unimplemented!("{:?} cannot yet be formatted (yet)", hdr.param),
            });
        }

        formatted
    }
}

/// A formatter for navigation solution
pub struct NavSolutionFormatter {
    pub filename: String,
    pub headers: Vec<NavSolutionHeader>,
    pub estimated_headers: StateFormatter,
    pub nominal_headers: StateFormatter,
    frames: HashMap<String, Frame>,
    cosm: Arc<Cosm>,
}

impl NavSolutionFormatter {
    /// ```
    /// extern crate nyx_space as nyx;
    /// use nyx::io::formatter::NavSolutionFormatter;
    /// use nyx::cosmic::Cosm;
    ///
    /// let cosm = Cosm::de438();
    /// // In this case, we're initializing the formatter to output the AoL and the eccentric anomaly in the EME2000 frame.
    /// let hdrs = vec!["estimate:AoL".to_string(), "nominal:ea:eme2000".to_string(), "delta_x".to_string()];
    /// NavSolutionFormatter::from_headers(hdrs, "nope".to_string(), cosm);
    /// ```
    pub fn from_headers(
        headers: Vec<String>,
        filename: String,
        cosm: Arc<Cosm>,
    ) -> Result<Self, NyxError> {
        let mut frames = HashMap::new();
        let mut hdrs = Vec::with_capacity(40);
        let mut est_hdrs = Vec::with_capacity(20);
        let mut nom_hdrs = Vec::with_capacity(20);
        // Rebuild the header tokens
        for hdr in &headers {
            let lowered = hdr.to_lowercase();
            let splt: Vec<&str> = lowered.split(':').collect();

            // Check the frames for the nominal and estimated state
            // The frames for the covariance are checked blow
            let frame_name = if splt.len() == 3 {
                // Check that the frame is valid
                let name = splt[2].to_owned();
                // Get the frame
                frames.insert(name.clone(), cosm.try_frame(&name)?);
                Some(name)
            } else if splt.len() == 2
                && splt[0].to_lowercase() != "epoch"
                && splt[0].to_lowercase() != "nominal"
                && splt[0].to_lowercase() != "estimate"
            {
                // See if this is a frame name
                match splt[1].to_lowercase().as_str() {
                    "ric" => frames.insert(splt[1].to_string(), Frame::RIC),
                    "rcn" => frames.insert(splt[1].to_string(), Frame::RCN),
                    "vnc" => frames.insert(splt[1].to_string(), Frame::VNC),
                    _ => {
                        if let Ok(frame) = cosm.try_frame(splt[1]) {
                            frames.insert(splt[1].to_string(), frame)
                        } else {
                            None
                        }
                    }
                };
                Some(splt[1].to_string())
            } else {
                None
            };

            match splt[0].to_lowercase().as_str() {
                "epoch" => {
                    hdrs.push(NavSolutionHeader::Epoch(if splt.len() == 2 {
                        EpochFormat::from_str(splt[1]).unwrap()
                    } else {
                        EpochFormat::GregorianUtc
                    }));
                }
                "delta_x" => hdrs.push(NavSolutionHeader::Delta_x),
                "delta_y" => hdrs.push(NavSolutionHeader::Delta_y),
                "delta_z" => hdrs.push(NavSolutionHeader::Delta_z),
                "delta_vx" => hdrs.push(NavSolutionHeader::Delta_vx),
                "delta_vy" => hdrs.push(NavSolutionHeader::Delta_vy),
                "delta_vz" => hdrs.push(NavSolutionHeader::Delta_vz),
                "cx_x" => hdrs.push(NavSolutionHeader::Cx_x { frame: frame_name }),
                "cy_x" => hdrs.push(NavSolutionHeader::Cy_x { frame: frame_name }),
                "cy_y" => hdrs.push(NavSolutionHeader::Cy_y { frame: frame_name }),
                "cz_x" => hdrs.push(NavSolutionHeader::Cz_x { frame: frame_name }),
                "cz_y" => hdrs.push(NavSolutionHeader::Cz_y { frame: frame_name }),
                "cz_z" => hdrs.push(NavSolutionHeader::Cz_z { frame: frame_name }),
                "cx_dot_x" => hdrs.push(NavSolutionHeader::Cx_dot_x { frame: frame_name }),
                "cx_dot_y" => hdrs.push(NavSolutionHeader::Cx_dot_y { frame: frame_name }),
                "cx_dot_z" => hdrs.push(NavSolutionHeader::Cx_dot_z { frame: frame_name }),
                "cx_dot_x_dot" => hdrs.push(NavSolutionHeader::Cx_dot_x_dot { frame: frame_name }),
                "cy_dot_x" => hdrs.push(NavSolutionHeader::Cy_dot_x { frame: frame_name }),
                "cy_dot_y" => hdrs.push(NavSolutionHeader::Cy_dot_y { frame: frame_name }),
                "cy_dot_z" => hdrs.push(NavSolutionHeader::Cy_dot_z { frame: frame_name }),
                "cy_dot_x_dot" => hdrs.push(NavSolutionHeader::Cy_dot_x_dot { frame: frame_name }),
                "cy_dot_y_dot" => hdrs.push(NavSolutionHeader::Cy_dot_y_dot { frame: frame_name }),
                "cz_dot_x" => hdrs.push(NavSolutionHeader::Cz_dot_x { frame: frame_name }),
                "cz_dot_y" => hdrs.push(NavSolutionHeader::Cz_dot_y { frame: frame_name }),
                "cz_dot_z" => hdrs.push(NavSolutionHeader::Cz_dot_z { frame: frame_name }),
                "cz_dot_x_dot" => hdrs.push(NavSolutionHeader::Cz_dot_x_dot { frame: frame_name }),
                "cz_dot_y_dot" => hdrs.push(NavSolutionHeader::Cz_dot_y_dot { frame: frame_name }),
                "cz_dot_z_dot" => hdrs.push(NavSolutionHeader::Cz_dot_z_dot { frame: frame_name }),
                "prediction" => hdrs.push(NavSolutionHeader::Prediction),
                "covar_position" => hdrs.push(NavSolutionHeader::Covar_pos),
                "covar_velocity" => hdrs.push(NavSolutionHeader::Covar_vel),
                "estimate" | "nominal" => {
                    let param = StateParameter::from_str(splt[1].to_lowercase().as_str())
                        .expect("Unknown paramater");

                    let state_hdr = StateHeader {
                        param,
                        frame_name,
                        epoch_fmt: None,
                    };

                    if splt[0] == "estimate" {
                        est_hdrs.push(state_hdr);
                    } else {
                        nom_hdrs.push(state_hdr);
                    }
                }
                _ => {
                    return Err(NyxError::CustomError(format!(
                        "unknown header `{}`",
                        splt[0]
                    )))
                }
            }
        }

        // Add the nominal and estimate headers (needed to add the header row)
        hdrs.push(NavSolutionHeader::EstimatedState(est_hdrs.clone()));
        hdrs.push(NavSolutionHeader::NominalState(nom_hdrs.clone()));

        Ok(Self {
            filename,
            headers: hdrs,
            nominal_headers: StateFormatter {
                filename: "file_should_not_exist".to_owned(),
                headers: nom_hdrs,
                frames: frames.clone(),
                cosm: cosm.clone(),
            },
            estimated_headers: StateFormatter {
                filename: "file_should_not_exist".to_owned(),
                headers: est_hdrs,
                frames: frames.clone(),
                cosm: cosm.clone(),
            },
            frames,
            cosm,
        })
    }

    /// Default headers are [Epoch (GregorianTai), X, Y, Z, VX, VY, VZ], where position is in km and velocity in km/s.
    pub fn default(filename: String, cosm: Arc<Cosm>) -> Self {
        let est_hdrs = vec![
            From::from(StateParameter::X),
            From::from(StateParameter::Y),
            From::from(StateParameter::Z),
            From::from(StateParameter::VX),
            From::from(StateParameter::VY),
            From::from(StateParameter::VZ),
        ];
        Self {
            filename,
            headers: vec![
                NavSolutionHeader::Epoch(EpochFormat::GregorianTai),
                NavSolutionHeader::Delta_x,
                NavSolutionHeader::Delta_y,
                NavSolutionHeader::Delta_z,
                NavSolutionHeader::Delta_vx,
                NavSolutionHeader::Delta_vy,
                NavSolutionHeader::Delta_vz,
                NavSolutionHeader::Cx_x { frame: None },
                NavSolutionHeader::Cy_y { frame: None },
                NavSolutionHeader::Cz_z { frame: None },
                NavSolutionHeader::Cx_dot_x_dot { frame: None },
                NavSolutionHeader::Cy_dot_y_dot { frame: None },
                NavSolutionHeader::Cz_dot_z_dot { frame: None },
                NavSolutionHeader::EstimatedState(est_hdrs.clone()),
                NavSolutionHeader::Prediction,
            ],
            nominal_headers: StateFormatter {
                filename: "file_should_not_exist".to_owned(),
                headers: Vec::new(),
                frames: HashMap::new(),
                cosm: cosm.clone(),
            },
            estimated_headers: StateFormatter {
                filename: "file_should_not_exist".to_owned(),
                headers: est_hdrs,
                frames: HashMap::new(),
                cosm: cosm.clone(),
            },
            frames: HashMap::new(),
            cosm,
        }
    }

    pub fn fmt<T: State, S: NavSolution<T>>(&self, sol: &S) -> Vec<String>
    where
        DefaultAllocator: Allocator<f64, <T as State>::Size>
            + Allocator<f64, <T as State>::Size, <T as State>::Size>
            + Allocator<f64, <T as State>::VecLength>,
    {
        // Start by computing the state in all of the frames needed
        let mut dcms = HashMap::new();
        let orbit = sol.orbital_state();
        for (name, frame) in &self.frames {
            let mapping_dcm = match frame {
                Frame::RCN | Frame::RIC | Frame::VNC => {
                    orbit.dcm6x6_from_traj_frame(*frame).unwrap()
                }
                _ => self
                    .cosm
                    .try_dcm_from_to(&orbit.frame, frame, orbit.dt)
                    .unwrap(),
            };
            dcms.insert(name.to_lowercase(), mapping_dcm);
        }

        let mut formatted = Vec::new();

        for hdr in &self.headers {
            match hdr {
                NavSolutionHeader::EstimatedState(_) => {
                    // The formatter is already initialized
                    for fmtval in self.estimated_headers.fmt(&sol.orbital_state()) {
                        formatted.push(fmtval);
                    }
                }
                NavSolutionHeader::NominalState(_) => {
                    // The formatter is already initialized
                    for fmtval in self.nominal_headers.fmt(&sol.expected_state()) {
                        formatted.push(fmtval);
                    }
                }
                NavSolutionHeader::Epoch(efmt) => formatted.push(efmt.format(sol.epoch())),
                NavSolutionHeader::Delta_x => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[0]))
                }
                NavSolutionHeader::Delta_y => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[1]))
                }
                NavSolutionHeader::Delta_z => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[2]))
                }
                NavSolutionHeader::Delta_vx => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[3]))
                }
                NavSolutionHeader::Delta_vy => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[4]))
                }
                NavSolutionHeader::Delta_vz => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[5]))
                }
                NavSolutionHeader::Cx_x { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(0, 0)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(0, 0))),
                },
                NavSolutionHeader::Cy_x { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(1, 0)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(1, 0))),
                },
                NavSolutionHeader::Cy_y { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(1, 1)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(1, 1))),
                },
                NavSolutionHeader::Cz_x { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(2, 0)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(2, 0))),
                },
                NavSolutionHeader::Cz_y { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(2, 1)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(2, 1))),
                },
                NavSolutionHeader::Cz_z { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(2, 2)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(2, 2))),
                },
                NavSolutionHeader::Cx_dot_x { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(3, 0)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(3, 0))),
                },
                NavSolutionHeader::Cx_dot_y { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(3, 1)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(3, 1))),
                },
                NavSolutionHeader::Cx_dot_z { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(3, 2)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(3, 2))),
                },
                NavSolutionHeader::Cx_dot_x_dot { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(3, 3)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(3, 3))),
                },
                NavSolutionHeader::Cy_dot_x { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(4, 0)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(4, 0))),
                },
                NavSolutionHeader::Cy_dot_y { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(4, 1)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(4, 1))),
                },
                NavSolutionHeader::Cy_dot_z { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(4, 2)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(4, 2))),
                },
                NavSolutionHeader::Cy_dot_x_dot { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(4, 3)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(4, 3))),
                },
                NavSolutionHeader::Cy_dot_y_dot { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(4, 4)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(4, 4))),
                },
                NavSolutionHeader::Cz_dot_x { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(5, 0)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(5, 0))),
                },
                NavSolutionHeader::Cz_dot_y { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(5, 1)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(5, 1))),
                },
                NavSolutionHeader::Cz_dot_z { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(5, 2)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(5, 2))),
                },
                NavSolutionHeader::Cz_dot_x_dot { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(5, 3)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(5, 3))),
                },
                NavSolutionHeader::Cz_dot_y_dot { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(5, 4)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(5, 4))),
                },
                NavSolutionHeader::Cz_dot_z_dot { frame } => match frame {
                    Some(frame_name) => {
                        // Rotate the covariance
                        let covar = dcms[frame_name] * sol.covar() * dcms[frame_name].transpose();
                        formatted.push(format!("{:.16e}", covar[(5, 5)]))
                    }
                    None => formatted.push(format!("{:.16e}", sol.covar_ij(5, 5))),
                },
                NavSolutionHeader::Prediction => formatted.push(format!("{}", sol.predicted())),
                NavSolutionHeader::Covar_pos => formatted.push(format!(
                    "{:.16e}",
                    (sol.covar_ij(0, 0).powi(2)
                        + sol.covar_ij(1, 1).powi(2)
                        + sol.covar_ij(2, 2).powi(2))
                    .sqrt()
                )),
                NavSolutionHeader::Covar_vel => formatted.push(format!(
                    "{:.16e}",
                    (sol.covar_ij(3, 3).powi(2)
                        + sol.covar_ij(4, 4).powi(2)
                        + sol.covar_ij(5, 5).powi(2))
                    .sqrt()
                )),
            };
        }

        formatted
    }
}
