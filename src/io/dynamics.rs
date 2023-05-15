/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2023 Christopher Rabotin <christopher.rabotin@gmail.com>

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

#[cfg(feature = "python")]
use pyo3::prelude::*;
use serde::{Deserialize, Serialize};

use super::{frames_from_str, frames_to_str, ConfigRepr};

use crate::cosmic::{Bodies, Frame};

#[derive(Debug, Deserialize, Serialize)]
pub struct HarmonicsSerde {
    pub frame: String,
    pub coeffs: String,
    pub degree: usize,
    pub order: usize,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct SrpSerde {
    pub phi: Option<f64>,
    #[serde(serialize_with = "frames_to_str", deserialize_with = "frames_from_str")]
    pub shadows: Vec<Frame>,
}

/// A representation of spacecraft dynamics that need to be used in Python with the spacecraft Propagator class.
#[derive(Debug, Deserialize, Serialize)]
#[cfg_attr(feature = "python", pyclass)]
pub struct DynamicsSerde {
    pub point_masses: Vec<Bodies>,
    pub harmonics: Option<Vec<HarmonicsSerde>>,
    pub srp: Option<SrpSerde>,
}

impl ConfigRepr for DynamicsSerde {}

#[test]
fn test_serde() {
    use crate::cosmic::Cosm;
    use crate::io::Configurable;
    use crate::md::prelude::SpacecraftDynamics;
    use std::collections::HashMap;

    let yaml = "
lofi:
  point_masses:
    - Sun
    - Earth
hifi:
  point_masses:
    - Sun
    - Earth
    - Luna
  harmonics:
    - frame: IAU Earth
      coeffs: data/JGM3.cof.gz
      degree: 10
      order: 10
  srp:
    phi: 1367.0
    shadows:
      - Sun J2000
      - Moon J2000
";

    let cosm = Cosm::de438();

    let mut dynamics_serde: HashMap<String, DynamicsSerde> = serde_yaml::from_str(yaml).unwrap();

    // Access the "hifi" dynamics
    let hifi_dynamics =
        SpacecraftDynamics::from_config(dynamics_serde.remove("hifi").unwrap(), cosm.clone())
            .unwrap();
    println!("hifi dynamics: {}", hifi_dynamics);

    // Access the "lofi" dynamics
    let lofi_dynamics =
        SpacecraftDynamics::from_config(dynamics_serde.remove("lofi").unwrap(), cosm).unwrap();
    println!("lofi dynamics: {}", lofi_dynamics);
}
