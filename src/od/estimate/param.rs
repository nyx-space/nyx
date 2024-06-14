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

use crate::md::StateParameter;

use super::EpochFormat;
use std::cmp::PartialEq;
use std::fmt;

/// Allowed headers, with an optional frame.
/// TODO: Support units
#[allow(non_camel_case_types)]
#[derive(Clone, Debug, PartialEq)]
pub enum EstimateParameter {
    /// The epoch in the specified format
    Epoch(EpochFormat),
    /// Parameters of the estimated state
    EstimatedState(Vec<StateParameter>),
    /// Parameters of the nominal state
    NominalState(Vec<StateParameter>),
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
    /// Norm of the velocity items of the covariance (Cx_dot_x_dot, Cy_dot_y_dot, Cz_dot_z_dot)
    Covar_vel,
}

impl fmt::Display for EstimateParameter {
    fn fmt(&self, fh: &mut fmt::Formatter) -> fmt::Result {
        match self {
            EstimateParameter::Epoch(efmt) => write!(fh, "Epoch:{efmt:?}"),
            EstimateParameter::EstimatedState(hdr) => {
                let mut seq = Vec::with_capacity(hdr.len());
                for element in hdr {
                    seq.push(format!("Estimate:{element}"));
                }
                write!(fh, "{}", seq.join(","))
            }
            EstimateParameter::NominalState(hdr) => {
                let mut seq = Vec::with_capacity(hdr.len());
                for element in hdr {
                    seq.push(format!("Nominal:{element}"));
                }
                write!(fh, "{}", seq.join(","))
            }
            EstimateParameter::Delta_x => write!(fh, "delta_x"),
            EstimateParameter::Delta_y => write!(fh, "delta_y"),
            EstimateParameter::Delta_z => write!(fh, "delta_z"),
            EstimateParameter::Delta_vx => write!(fh, "delta_vx"),
            EstimateParameter::Delta_vy => write!(fh, "delta_vy"),
            EstimateParameter::Delta_vz => write!(fh, "delta_vz"),
            EstimateParameter::Cx_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cx_x:{f}")
                } else {
                    write!(fh, "cx_x")
                }
            }
            EstimateParameter::Cy_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_x:{f}")
                } else {
                    write!(fh, "cy_x")
                }
            }
            EstimateParameter::Cy_y { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_y:{f}")
                } else {
                    write!(fh, "cy_y")
                }
            }
            EstimateParameter::Cz_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_x:{f}")
                } else {
                    write!(fh, "cz_x")
                }
            }
            EstimateParameter::Cz_y { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_y:{f}")
                } else {
                    write!(fh, "cz_y")
                }
            }
            EstimateParameter::Cz_z { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_z:{f}")
                } else {
                    write!(fh, "cz_z")
                }
            }
            EstimateParameter::Cx_dot_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cx_dot_x:{f}")
                } else {
                    write!(fh, "cx_dot_x")
                }
            }
            EstimateParameter::Cx_dot_y { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cx_dot_y:{f}")
                } else {
                    write!(fh, "cx_dot_y")
                }
            }
            EstimateParameter::Cx_dot_z { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cx_dot_z:{f}")
                } else {
                    write!(fh, "cx_dot_z")
                }
            }
            EstimateParameter::Cx_dot_x_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cx_dot_x_dot:{f}")
                } else {
                    write!(fh, "cx_dot_x_dot")
                }
            }
            EstimateParameter::Cy_dot_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_dot_x:{f}")
                } else {
                    write!(fh, "cy_dot_x")
                }
            }
            EstimateParameter::Cy_dot_y { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_dot_y:{f}")
                } else {
                    write!(fh, "cy_dot_y")
                }
            }
            EstimateParameter::Cy_dot_z { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_dot_z:{f}")
                } else {
                    write!(fh, "cy_dot_z")
                }
            }
            EstimateParameter::Cy_dot_x_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_dot_x_dot:{f}")
                } else {
                    write!(fh, "cy_dot_x_dot")
                }
            }
            EstimateParameter::Cy_dot_y_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cy_dot_y_dot:{f}")
                } else {
                    write!(fh, "cy_dot_y_dot")
                }
            }
            EstimateParameter::Cz_dot_x { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_x:{f}")
                } else {
                    write!(fh, "cz_dot_x")
                }
            }
            EstimateParameter::Cz_dot_y { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_y:{f}")
                } else {
                    write!(fh, "cz_dot_y")
                }
            }
            EstimateParameter::Cz_dot_z { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_z:{f}")
                } else {
                    write!(fh, "cz_dot_z")
                }
            }
            EstimateParameter::Cz_dot_x_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_x_dot:{f}")
                } else {
                    write!(fh, "cz_dot_x_dot")
                }
            }
            EstimateParameter::Cz_dot_y_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_y_dot:{f}")
                } else {
                    write!(fh, "cz_dot_y_dot")
                }
            }
            EstimateParameter::Cz_dot_z_dot { frame } => {
                if let Some(f) = frame {
                    write!(fh, "cz_dot_z_dot:{f}")
                } else {
                    write!(fh, "cz_dot_z_dot")
                }
            }
            EstimateParameter::Prediction => write!(fh, "prediction"),
            EstimateParameter::Covar_pos => write!(fh, "covar_position"),
            EstimateParameter::Covar_vel => write!(fh, "covar_velocity"),
        }
    }
}
