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
