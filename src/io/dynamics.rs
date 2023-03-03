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

// use std::collections::HashMap;

// use serde_derive::{Deserialize, Serialize};
// #[derive(Debug, Serialize, Deserialize, PartialEq)]
// pub struct DynamicsSerde {
//     #[serde(transparent)]
//     inner: HashMap<String, String>,
// }

use serde::{Deserialize, Serialize};
// use std::collections::HashMap;

#[derive(Debug, Deserialize, Serialize)]
pub struct Harmonics {
    frame: String,
    coeffs: String,
    degree: usize,
    order: usize,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Srp {
    phi: Option<f64>,
    #[serde(flatten)]
    shadows: Vec<String>,
}

#[derive(Debug, Deserialize, Serialize)]
pub struct Dynamics {
    point_masses: Vec<String>,
    harmonics: Vec<Harmonics>,
    srp: Option<Srp>,
}

// fn main() -> Result<(), Box<dyn std::error::Error>> {
//     let yaml = "
// lofi:
//   point_masses:
//     - Sun
//     - Earth
// hifi:
//   point_masses:
//     - Sun
//     - Earth
//     - Moon
//   harmonics:
//     - frame: IAU Earth
//       coeffs: data/JGM3.cof.gz
//       degree: 10
//       order: 10
//   srp:
//     phi: 1367.0
//     shadows:
//       - Sun
//       - Moon
// ";

//     let mut dynamics: HashMap<String, Dynamics> = serde_yaml::from_str(yaml)?;

//     // Access the "hifi" dynamics
//     let hifi_dynamics = dynamics.remove("hifi").unwrap();
//     println!("hifi dynamics: {:?}", hifi_dynamics);

//     // Access the "lofi" dynamics
//     let lofi_dynamics = dynamics.remove("lofi").unwrap();
//     println!("lofi dynamics: {:?}", lofi_dynamics);

//     Ok(())
// }
