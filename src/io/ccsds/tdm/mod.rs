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

/*
 * This file defines how to parse a XML as TDM.
 * This is coded from CCSDS 503.0-B-2, published in June 2020
 */

extern crate serde;
extern crate serde_derive;
extern crate yaserde;
extern crate yaserde_derive;

use self::serde_derive::Deserialize;
use yaserde_derive::{YaDeserialize, YaSerialize};
// use self::serde_with::{CommaSeparator, SpaceSeparator, StringWithSeparator};

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
pub struct Tdm {
    #[yaserde(child)]
    id: String,
    #[yaserde(child)]
    version: f64,
    #[yaserde(child)]
    header: Header,
    #[yaserde(child)]
    body: Body,
}

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]

struct Header {
    #[yaserde(child)]
    comment: Vec<String>,
    #[yaserde(child)]
    creation_date: String,
    #[yaserde(child)]
    originator: String,
}

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
struct Body {
    #[yaserde(child)]
    segment: Vec<Segment>,
}

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
struct Segment {
    #[yaserde(child)]
    metadata: Metadata,
    // data: Data,
}

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
struct Metadata {
    #[yaserde(child, rename = "COMMENT")]
    comment: Vec<String>,
    #[yaserde(child, rename = "TRACK_ID")]
    track_id: Option<String>,
    data_types: Option<String>, // Check Table 3.5 for valid values
    #[yaserde(child, rename = "TIME_SYSTEM")]
    time_system: TimeSystem,
    // time_system: TimeSystemTag,
    // start_time: Option<String>,
    // stop_time: Option<String>,
    // participant_1: Option<String>,
    // participant_2: Option<String>,
    // participant_3: Option<String>,
    // participant_4: Option<String>,
    // participant_5: Option<String>,
    // #[serde(rename = "$value")]
    // mode: Option<TrackingMode>,
    // // #[serde(with = "StringWithSeparator::<CommaSeparator>")]
    // path: Option<Vec<u8>>,
    // path_1: Option<Vec<u8>>,
    // path_2: Option<Vec<u8>>,
    // ephemeris_1: Option<String>,
    // ephemeris_2: Option<String>,
    // ephemeris_3: Option<String>,
    // ephemeris_4: Option<String>,
    // ephemeris_5: Option<String>,
    // transmit_band: Option<Band>,
    // receive_band: Option<Band>,
    // turnaround_denominator: Option<i32>,
    // turnaround_numerator: Option<i32>,
    // timetag_ref: Option<TimetagRef>,
    // integration_interval: Option<f64>,
    // integration_ref: Option<IntegrationRef>,
    // freq_offset: Option<f64>,
    // range_mode: Option<RangeMode>,
    // range_modulus: Option<f64>,
    // range_units: Option<RangeUnit>,
    // angle_type: Option<String>,
    // reference_frame: Option<String>,
    // interpolation: Option<Interpolation>,
    // interpolation_degree: Option<u32>,
    // doppler_count_bias: Option<f64>,
    // doppler_count_scale: Option<u32>,
    // transmit_delay_1: Option<f64>,
    // transmit_delay_2: Option<f64>,
    // transmit_delay_3: Option<f64>,
    // transmit_delay_4: Option<f64>,
    // transmit_delay_5: Option<f64>,
    // receive_delay_1: Option<f64>,
    // receive_delay_2: Option<f64>,
    // receive_delay_3: Option<f64>,
    // receive_delay_4: Option<f64>,
    // receive_delay_5: Option<f64>,
    // data_quality: Option<DataQuality>,
    // correction_angle_1: Option<f64>,
    // correction_angle_2: Option<f64>,
    // correction_doppler: Option<f64>,
    // correction_mag: Option<f64>,
    // correction_range: Option<f64>,
    // correction_rcs: Option<f64>,
    // correction_receive: Option<f64>,
    // correction_transmit: Option<f64>,
    // correction_aberration_yearly: Option<f64>,
    // correction_aberration_diurnal: Option<f64>,
    // corrections_applied: Option<YesNo>,
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
enum TimeSystem {
    Utc,
    Tai,
    Gps,
    Sclk,
}

impl Default for TimeSystem {
    fn default() -> Self {
        TimeSystem::Utc
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
#[allow(non_camel_case_types)]
enum TrackingMode {
    Sequential,
    Single_Diff,
}

impl Default for TrackingMode {
    fn default() -> Self {
        Self::Sequential
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
#[allow(clippy::upper_case_acronyms)]
enum Band {
    S,
    X,
    Ka,
    L,
    UHF,
    GREEN,
}

impl Default for Band {
    fn default() -> Self {
        Self::S
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
enum TimetagRef {
    Transmit,
    Receive,
}

impl Default for TimetagRef {
    fn default() -> Self {
        Self::Transmit
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
enum IntegrationRef {
    Start,
    Middle,
    End,
}

impl Default for IntegrationRef {
    fn default() -> Self {
        Self::Start
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
#[allow(non_camel_case_types)]
enum RangeMode {
    Coherent,
    Constant,
    One_Way,
}

impl Default for RangeMode {
    fn default() -> Self {
        Self::Coherent
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
#[allow(non_camel_case_types)]
enum RangeUnit {
    km,
    s,
    RU,
}

impl Default for RangeUnit {
    fn default() -> Self {
        Self::km
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
#[allow(clippy::upper_case_acronyms)]
enum AngleType {
    /// Azimuth, Elevation (local horizontal)
    AZEL,
    /// Right ascension, declination (must be referenced to inertial frame)
    RADEC,
    /// x-east, y-north
    XEYN,
    /// x-south, y-east
    XSYE,
}

impl Default for AngleType {
    fn default() -> Self {
        Self::AZEL
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
enum Interpolation {
    Hermite,
    Lagrange,
    Linear,
}

impl Default for Interpolation {
    fn default() -> Self {
        Self::Hermite
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
enum DataQuality {
    Raw,
    Validated,
    Degraded,
}

impl Default for DataQuality {
    fn default() -> Self {
        Self::Raw
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
enum YesNo {
    Yes,
    No,
}

impl Default for YesNo {
    fn default() -> Self {
        Self::Yes
    }
}
