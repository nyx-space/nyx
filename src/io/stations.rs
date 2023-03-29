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

use super::ConfigError;
use super::{ConfigRepr, Configurable};
use crate::md::ui::Cosm;
use crate::od::prelude::GroundStation;
use serde_derive::{Deserialize, Serialize};
use std::fmt::Debug;
use std::sync::Arc;

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct StationSerde {
    pub name: String,
    pub frame: Option<String>,
    pub elevation_mask_deg: f64,
    pub range_noise_km: f64,
    pub range_rate_noise_km_s: f64,
    pub latitude_deg: Option<f64>,
    pub longitude_deg: Option<f64>,
    pub height_km: Option<f64>,
    pub inherit: Option<String>,
}

impl ConfigRepr for StationSerde {}

impl Configurable for GroundStation {
    type IntermediateRepr = StationSerde;

    fn from_config(cfg: Self::IntermediateRepr, cosm: Arc<Cosm>) -> Result<Self, ConfigError>
    where
        Self: Sized,
    {
        let iau_earth = cosm.frame("IAU Earth");
        let gs = if let Some(base) = &cfg.inherit {
            match base.to_lowercase().as_str() {
                "dss13" => GroundStation::dss13_goldstone(
                    cfg.elevation_mask_deg,
                    cfg.range_noise_km,
                    cfg.range_rate_noise_km_s,
                    iau_earth,
                ),
                "dss34" => GroundStation::dss34_canberra(
                    cfg.elevation_mask_deg,
                    cfg.range_noise_km,
                    cfg.range_rate_noise_km_s,
                    iau_earth,
                ),
                "dss65" => GroundStation::dss65_madrid(
                    cfg.elevation_mask_deg,
                    cfg.range_noise_km,
                    cfg.range_rate_noise_km_s,
                    iau_earth,
                ),
                _ => {
                    return Err(ConfigError::InvalidConfig(format!(
                        "unknown built-in station {base}"
                    )))
                }
            }
        } else {
            let frame_name = cfg
                .frame
                .clone()
                .or_else(|| Some("IAU Earth".to_string()))
                .unwrap();

            GroundStation::from_noise_values(
                cfg.name.clone(),
                cfg.elevation_mask_deg,
                cfg.latitude_deg.unwrap(),
                cfg.longitude_deg.unwrap(),
                cfg.height_km.unwrap(),
                cfg.range_noise_km,
                cfg.range_rate_noise_km_s,
                cosm.try_frame(&frame_name).map_or_else(
                    |e| {
                        error!("{e} - using {iau_earth}");
                        iau_earth
                    },
                    |f| f,
                ),
            )
        };

        Ok(gs)
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, ConfigError> {
        Ok(StationSerde {
            name: self.name.clone(),
            frame: Some(self.frame.to_string()),
            elevation_mask_deg: self.elevation_mask_deg,
            range_noise_km: self.range_noise.mean(),
            range_rate_noise_km_s: self.range_rate_noise.mean(),
            latitude_deg: Some(self.latitude_deg),
            longitude_deg: Some(self.longitude_deg),
            height_km: Some(self.height_km),
            inherit: None,
        })
    }
}

#[test]
fn test_load_single() {
    use std::env;
    use std::path::PathBuf;

    // Get the path to the root directory of the current Cargo project
    let test_data: PathBuf = [
        env::var("CARGO_MANIFEST_DIR").unwrap(),
        "data".to_string(),
        "tests".to_string(),
        "config".to_string(),
        "one_ground_station.yaml".to_string(),
    ]
    .iter()
    .collect();

    assert!(test_data.exists(), "Could not find the test data");

    let gs = StationSerde::load(test_data).unwrap();

    dbg!(&gs);

    let expected_gs = StationSerde {
        name: "Demo ground station".to_string(),
        frame: Some("IAU Earth".to_string()),
        elevation_mask_deg: 5.0,
        range_noise_km: 1e-3,
        range_rate_noise_km_s: 1e-5,
        latitude_deg: Some(2.3522),
        longitude_deg: Some(48.8566),
        height_km: Some(0.4),
        inherit: None,
    };

    assert_eq!(expected_gs, gs);
}

#[test]
fn test_load_many() {
    use crate::cosmic::Cosm;
    use std::env;
    use std::path::PathBuf;

    // Get the path to the root directory of the current Cargo project

    let test_file: PathBuf = [
        env::var("CARGO_MANIFEST_DIR").unwrap(),
        "data".to_string(),
        "tests".to_string(),
        "config".to_string(),
        "many_ground_stations.yaml".to_string(),
    ]
    .iter()
    .collect();

    let stations = StationSerde::load_many(test_file).unwrap();

    dbg!(&stations);

    let expected = vec![
        StationSerde {
            name: "Demo ground station".to_string(),
            frame: Some("IAU Earth".to_string()),
            elevation_mask_deg: 5.0,
            range_noise_km: 1e-3,
            range_rate_noise_km_s: 1e-5,
            latitude_deg: Some(2.3522),
            longitude_deg: Some(48.8566),
            height_km: Some(0.4),
            inherit: None,
        },
        StationSerde {
            name: "Inherited".to_string(),
            frame: None,
            elevation_mask_deg: 5.0,
            range_noise_km: 1e-3,
            range_rate_noise_km_s: 1e-5,
            inherit: Some("DSS34".to_string()),
            latitude_deg: None,
            longitude_deg: None,
            height_km: None,
        },
    ];

    assert_eq!(expected, stations);

    let cosm = Cosm::de438();

    // Convert into a ground station
    let expected_name = "Demo ground station".to_string();
    let expected_el_mask = 5.0;
    let expected_lat = 2.3522;
    let expected_long = 48.8566;
    let expected_height = 0.4;
    let expected_frame = cosm.frame("IAU Earth");

    let gs = GroundStation::from_config(stations[0].clone(), cosm.clone()).unwrap();
    dbg!(&gs);
    assert_eq!(expected_name, gs.name);
    assert_eq!(expected_el_mask, gs.elevation_mask_deg);
    assert_eq!(expected_lat, gs.latitude_deg);
    assert_eq!(expected_long, gs.longitude_deg);
    assert_eq!(expected_frame, gs.frame);
    assert_eq!(expected_height, gs.height_km);
}
