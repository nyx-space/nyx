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

use anise::analysis::prelude::OrbitalElement;
use arrow::datatypes::{DataType, Field};
use core::fmt;

use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Common state parameters
#[allow(non_camel_case_types, clippy::upper_case_acronyms)]
#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize)]
pub enum StateParameter {
    Element(OrbitalElement),
    /// B-Plane B⋅R
    BdotR,
    /// B-Plane B⋅T
    BdotT,
    /// B-Plane LTOF
    BLTOF,
    /// Coefficient of drag
    Cd,
    /// Coefficient of reflectivity
    Cr,
    /// Dry mass (kg)
    DryMass,
    /// The epoch of the state
    Epoch,
    /// Return the guidance mode of the spacecraft
    GuidanceMode,
    /// Specific impulse (isp) in seconds
    Isp,
    /// prop mass in kilograms
    PropMass,
    /// Thrust (Newtons)
    Thrust,
    /// Total mass
    TotalMass,
}

impl StateParameter {
    /// Returns the default event finding precision in the unit of that parameter
    pub fn default_event_precision(&self) -> f64 {
        match self {
            Self::Element(el) => {
                if el == &OrbitalElement::Period {
                    1e-1
                } else {
                    1e-3
                }
            }
            Self::BdotR | Self::BdotT => 1e-3,

            // Special
            Self::DryMass | Self::PropMass => 1e-3,
            _ => unimplemented!("{self} cannot be used for targeting"),
        }
    }

    /// Returns whether this parameter is of the B-Plane kind
    pub const fn is_b_plane(&self) -> bool {
        matches!(&self, Self::BdotR | Self::BdotT | Self::BLTOF)
    }

    /// Returns whether this is an orbital parameter
    pub const fn is_orbital(&self) -> bool {
        !self.is_for_spacecraft() && !matches!(self, Self::Element(..))
    }

    /// Returns whether this parameter is only applicable to a spacecraft state
    pub const fn is_for_spacecraft(&self) -> bool {
        matches!(
            &self,
            Self::DryMass
                | Self::PropMass
                | Self::Cr
                | Self::Cd
                | Self::Isp
                | Self::GuidanceMode
                | Self::Thrust
        )
    }

    pub const fn unit(&self) -> &'static str {
        match self {
            Self::Element(e) => e.unit(),
            Self::BdotR | Self::BdotT => "km",

            Self::DryMass | Self::PropMass => "kg",
            Self::Isp => "isp",
            Self::Thrust => "N",
            _ => "",
        }
    }
}

impl StateParameter {
    /// Returns the parquet field of this parameter
    pub(crate) fn to_field(self, more_meta: Option<Vec<(String, String)>>) -> Field {
        self.to_field_generic(false, more_meta)
    }

    /// Returns the parquet field of this parameter
    pub(crate) fn to_cov_field(self, more_meta: Option<Vec<(String, String)>>) -> Field {
        self.to_field_generic(true, more_meta)
    }

    /// Returns the parquet field of this parameter
    fn to_field_generic(self, is_sigma: bool, more_meta: Option<Vec<(String, String)>>) -> Field {
        let mut meta = HashMap::new();
        meta.insert("unit".to_string(), self.unit().to_string());
        if let Some(more_data) = more_meta {
            for (k, v) in more_data {
                meta.insert(k, v);
            }
        }

        Field::new(
            if is_sigma {
                format!("Sigma {self}")
            } else {
                format!("{self}")
            },
            if self == Self::GuidanceMode {
                DataType::Utf8
            } else {
                DataType::Float64
            },
            false,
        )
        .with_metadata(meta)
    }
}

impl fmt::Display for StateParameter {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let repr = match *self {
            Self::Element(e) => return write!(f, "{e}"),
            Self::BLTOF => "BLToF",
            Self::BdotR => "BdotR",
            Self::BdotT => "BdotT",
            Self::Cd => "cd",
            Self::Cr => "cr",
            Self::DryMass => "dry_mass",
            Self::Epoch => "epoch",
            Self::GuidanceMode => "guidance_mode",
            Self::Isp => "isp",
            Self::PropMass => "prop_mass",
            Self::Thrust => "thrust",
            Self::TotalMass => "total_mass",
        };
        let unit = if self.unit().is_empty() {
            String::new()
        } else {
            format!(" ({})", self.unit())
        };
        write!(f, "{repr}{unit}")
    }
}
