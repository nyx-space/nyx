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

use std::sync::Arc;

use super::serde_derive::{Deserialize, Serialize};
use super::{Configurable, ParsingError};
use crate::od::ui::GroundStation;

#[derive(Serialize, Deserialize)]
pub struct StationSerde {
    pub name: String,
    pub frame_name: Option<String>,
    pub elevation: f64,
    pub range_noise: f64,
    pub range_rate_noise: f64,
    pub latitude: Option<f64>,
    pub longitude: Option<f64>,
    pub height: Option<f64>,
    pub inherit: Option<String>,
}

impl<'a> Configurable<'a> for GroundStation {
    type IntermediateRepr = StationSerde;

    fn from_config(
        cfg: &Self::IntermediateRepr,
        cosm: Arc<super::odp::Cosm>,
    ) -> Result<Self, ParsingError>
    where
        Self: Sized,
    {
        let gs = if let Some(base) = &cfg.inherit {
            match base.to_lowercase().as_str() {
                "dss13" => GroundStation::dss13_goldstone(
                    cfg.elevation,
                    cfg.range_noise,
                    cfg.range_rate_noise,
                    cosm.clone(), // Dammit it.
                ),
                "dss34" => GroundStation::dss34_canberra(
                    cfg.elevation,
                    cfg.range_noise,
                    cfg.range_rate_noise,
                    cosm.clone(),
                ),
                "dss65" => GroundStation::dss65_madrid(
                    cfg.elevation,
                    cfg.range_noise,
                    cfg.range_rate_noise,
                    cosm.clone(),
                ),
                _ => return Err(ParsingError::OD(format!("unknown base station `{base}`",))),
            }
        } else {
            let frame_name = cfg
                .frame_name
                .clone()
                .or_else(|| Some("IAU Earth".to_string()))
                .unwrap();

            GroundStation::from_noise_values(
                cfg.name.clone(),
                cfg.elevation,
                cfg.latitude.unwrap(),
                cfg.longitude.unwrap(),
                cfg.height.unwrap(),
                cfg.range_noise,
                cfg.range_rate_noise,
                cosm.frame(&frame_name),
                cosm.clone(),
            )
        };

        Ok(gs)
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, ParsingError> {
        todo!()
    }
}
