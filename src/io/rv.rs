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

use super::serde_derive::Deserialize;

#[derive(Deserialize, PartialEq)]
#[serde(rename_all = "camelCase")]
pub enum Distribution {
    Normal { mean: f64, std_dev: f64 },
    Cauchy { median: f64, scale: f64 },
    Exponential { lambda: f64 },
    Poisson { lambda: f64 },
}

#[test]
fn test_deser_distr() {
    extern crate toml;

    let _std_norm: Distribution = toml::from_str(
        r#"[normal]
        mean = 0.0
        std_dev = 1.0"#,
    )
    .unwrap();

    let _cauchy: Distribution = toml::from_str(
        r#"[cauchy]
        median = 0.0
        scale = 1.0"#,
    )
    .unwrap();

    let _my_exp: Distribution = toml::from_str(
        r#"[exponential]
        lambda = 0.5"#,
    )
    .unwrap();

    let _fish: Distribution = toml::from_str(
        r#"[poisson]
        lambda = 10.0
    "#,
    )
    .unwrap();
}

#[test]
fn test_deser_distr_multi() {
    use std::collections::HashMap;
    extern crate toml;

    #[derive(Deserialize)]
    struct MapRv {
        pub rvs: HashMap<String, Distribution>,
    }

    let _as_map: MapRv = toml::from_str(
        r#"rvs = { one = { normal = { mean = 0.0, std_dev = 0.2 } }, two = { poisson = { lambda = 10.0 } } }"#,
    )
    .unwrap();
    println!("{:?}", _as_map.rvs.keys());
}
