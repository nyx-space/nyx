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
use crate::NyxError;
use yaserde_derive::{YaDeserialize, YaSerialize};
// use self::serde_with::{CommaSeparator, SpaceSeparator, StringWithSeparator};

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
pub struct Tdm {
    #[yaserde(attribute)]
    pub id: String,
    #[yaserde(attribute)]
    pub version: f64,
    #[yaserde(child)]
    pub header: Header,
    #[yaserde(child)]
    pub body: Body,
}

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
pub struct Header {
    #[yaserde(child, rename = "COMMENT")]
    comment: Vec<String>,
    #[yaserde(child, rename = "CREATION_DATE")]
    pub creation_date: String,
    #[yaserde(child, rename = "ORIGINATOR")]
    pub originator: String,
}

impl Header {
    pub fn comments(&self) -> String {
        self.comment.join(" ")
    }
}

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
pub struct Body {
    #[yaserde(child)]
    pub segment: Vec<Segment>,
}

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
pub struct Segment {
    #[yaserde(child)]
    pub metadata: Metadata,
    // data: Data,
}

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
pub struct Metadata {
    #[yaserde(child, rename = "COMMENT")]
    comment: Vec<String>,
    #[yaserde(child, rename = "TRACK_ID")]
    pub track_id: Option<String>,
    #[yaserde(child, rename = "DATA_TYPES")]
    data_types: Option<String>, // Check Table 3.5 for valid values
    #[yaserde(child, rename = "TIME_SYSTEM")]
    pub time_system: TimeSystem,
    #[yaserde(child, rename = "START_TIME")]
    pub start_time: Option<String>,
    #[yaserde(child, rename = "STOP_TIME")]
    pub stop_time: Option<String>,
    #[yaserde(child, rename = "PARTICIPANT_1")]
    pub participant_1: Option<String>,
    #[yaserde(child, rename = "PARTICIPANT_2")]
    pub participant_2: Option<String>,
    #[yaserde(child, rename = "PARTICIPANT_3")]
    pub participant_3: Option<String>,
    #[yaserde(child, rename = "PARTICIPANT_4")]
    pub participant_4: Option<String>,
    #[yaserde(child, rename = "PARTICIPANT_5")]
    pub participant_5: Option<String>,
    #[yaserde(child, rename = "MODE")]
    pub mode: Option<TrackingMode>,
    #[yaserde(child, rename = "PATH")]
    pub path: Option<String>,
    #[yaserde(child, rename = "PATH_1")]
    pub path_1: Option<String>,
    #[yaserde(child, rename = "PATH_2")]
    pub path_2: Option<String>,
    #[yaserde(child, rename = "EPHEMERIS_1")]
    pub ephemeris_1: Option<String>,
    #[yaserde(child, rename = "EPHEMERIS_2")]
    pub ephemeris_2: Option<String>,
    #[yaserde(child, rename = "EPHEMERIS_3")]
    pub ephemeris_3: Option<String>,
    #[yaserde(child, rename = "EPHEMERIS_4")]
    pub ephemeris_4: Option<String>,
    #[yaserde(child, rename = "EPHEMERIS_5")]
    pub ephemeris_5: Option<String>,
    #[yaserde(child, rename = "TRANSMIT_BAND")]
    pub transmit_band: Option<Band>,
    #[yaserde(child, rename = "RECEIVE_BAND")]
    pub receive_band: Option<Band>,
    #[yaserde(child, rename = "TURNAROUND_DENOMINATOR")]
    pub turnaround_denominator: Option<i32>,
    #[yaserde(child, rename = "TURNAROUND_NUMERATOR")]
    pub turnaround_numerator: Option<i32>,
    #[yaserde(child, rename = "TIMETAG_REF")]
    pub timetag_ref: Option<TimetagRef>,
    #[yaserde(child, rename = "INTEGRATION_INTERVAL")]
    pub integration_interval: Option<f64>,
    #[yaserde(child, rename = "INTEGRATION_REF")]
    pub integration_ref: Option<IntegrationRef>,
    #[yaserde(child, rename = "FREQ_OFFSET")]
    pub freq_offset: Option<f64>,
    #[yaserde(child, rename = "RANGE_MODE")]
    pub range_mode: Option<RangeMode>,
    #[yaserde(child, rename = "RANGE_MODULUS")]
    pub range_modulus: Option<f64>,
    #[yaserde(child, rename = "RANGE_UNITS")]
    pub range_units: Option<RangeUnit>,
    #[yaserde(child, rename = "ANGLE_TYPE")]
    pub angle_type: Option<String>,
    #[yaserde(child, rename = "REFERENCE_FRAME")]
    pub reference_frame: Option<String>,
    #[yaserde(child, rename = "INTERPOLATION")]
    pub interpolation: Option<Interpolation>,
    #[yaserde(child, rename = "INTERPOLATION_DEGREE")]
    pub interpolation_degree: Option<u32>,
    #[yaserde(child, rename = "DOPPLER_COUNT_BIAS")]
    pub doppler_count_bias: Option<f64>,
    #[yaserde(child, rename = "DOPPLER_COUNT_SCALE")]
    pub doppler_count_scale: Option<u32>,
    #[yaserde(child, rename = "TRANSMIT_DELAY_1")]
    pub transmit_delay_1: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_DELAY_2")]
    pub transmit_delay_2: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_DELAY_3")]
    pub transmit_delay_3: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_DELAY_4")]
    pub transmit_delay_4: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_DELAY_5")]
    pub transmit_delay_5: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_DELAY_1")]
    pub receive_delay_1: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_DELAY_2")]
    pub receive_delay_2: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_DELAY_3")]
    pub receive_delay_3: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_DELAY_4")]
    pub receive_delay_4: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_DELAY_5")]
    pub receive_delay_5: Option<f64>,
    #[yaserde(child, rename = "DATA_QUALITY")]
    pub data_quality: Option<DataQuality>,
    #[yaserde(child, rename = "CORRECTION_ANGLE_1")]
    pub correction_angle_1: Option<f64>,
    #[yaserde(child, rename = "CORRECTION_ANGLE_2")]
    pub correction_angle_2: Option<f64>,
    #[yaserde(child, rename = "CORRECTION_DOPPLER")]
    pub correction_doppler: Option<f64>,
    #[yaserde(child, rename = "CORRECTION_MAG")]
    pub correction_mag: Option<f64>,
    #[yaserde(child, rename = "CORRECTION_RANGE")]
    pub correction_range: Option<f64>,
    #[yaserde(child, rename = "CORRECTION_RCS")]
    pub correction_rcs: Option<f64>,
    #[yaserde(child, rename = "CORRECTION_RECEIVE")]
    pub correction_receive: Option<f64>,
    #[yaserde(child, rename = "CORRECTION_TRANSMIT")]
    pub correction_transmit: Option<f64>,
    #[yaserde(child, rename = "CORRECTION_ABERRATION_YEARLY")]
    pub correction_aberration_yearly: Option<f64>,
    #[yaserde(child, rename = "CORRECTION_ABERRATION_DIURNAL")]
    pub correction_aberration_diurnal: Option<f64>,
    #[yaserde(child, rename = "CORRECTIONS_APPLIED")]
    pub corrections_applied: Option<YesNo>,
}

impl Metadata {
    pub fn comments(&self) -> String {
        self.comment.join(" ")
    }

    pub fn participant(&self, n: usize) -> Result<Participant, NyxError> {
        if n == 0 || n > 5 {
            Err(NyxError::CCSDS(
                "Valid participants numbered 1 through 5".to_string(),
            ))
        } else {
            match n {
                1 => match &self.participant_1 {
                    Some(participant) => Ok(Participant {
                        name: participant.clone(),
                        ephemeris: self.ephemeris_1.clone(),
                        transmit_delay: self.transmit_delay_1.unwrap_or(0.0),
                        receive_delay: self.receive_delay_1.unwrap_or(0.0),
                    }),
                    None => Err(NyxError::CCSDS(format!("No participant #{}", n))),
                },
                2 => match &self.participant_2 {
                    Some(participant) => Ok(Participant {
                        name: participant.clone(),
                        ephemeris: self.ephemeris_2.clone(),
                        transmit_delay: self.transmit_delay_2.unwrap_or(0.0),
                        receive_delay: self.receive_delay_2.unwrap_or(0.0),
                    }),
                    None => Err(NyxError::CCSDS(format!("No participant #{}", n))),
                },
                3 => match &self.participant_3 {
                    Some(participant) => Ok(Participant {
                        name: participant.clone(),
                        ephemeris: self.ephemeris_3.clone(),
                        transmit_delay: self.transmit_delay_3.unwrap_or(0.0),
                        receive_delay: self.receive_delay_3.unwrap_or(0.0),
                    }),
                    None => Err(NyxError::CCSDS(format!("No participant #{}", n))),
                },
                4 => match &self.participant_4 {
                    Some(participant) => Ok(Participant {
                        name: participant.clone(),
                        ephemeris: self.ephemeris_4.clone(),
                        transmit_delay: self.transmit_delay_4.unwrap_or(0.0),
                        receive_delay: self.receive_delay_4.unwrap_or(0.0),
                    }),
                    None => Err(NyxError::CCSDS(format!("No participant #{}", n))),
                },
                5 => match &self.participant_5 {
                    Some(participant) => Ok(Participant {
                        name: participant.clone(),
                        ephemeris: self.ephemeris_5.clone(),
                        transmit_delay: self.transmit_delay_5.unwrap_or(0.0),
                        receive_delay: self.receive_delay_5.unwrap_or(0.0),
                    }),
                    None => Err(NyxError::CCSDS(format!("No participant #{}", n))),
                },
                _ => unreachable!(),
            }
        }
    }
}

pub struct Participant {
    pub name: String,
    pub ephemeris: Option<String>,
    pub transmit_delay: f64,
    pub receive_delay: f64,
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
pub enum TimeSystem {
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
pub enum TrackingMode {
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
pub enum Band {
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
pub enum TimetagRef {
    Transmit,
    Receive,
}

impl Default for TimetagRef {
    fn default() -> Self {
        Self::Transmit
    }
}

#[derive(Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
pub enum IntegrationRef {
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
pub enum RangeMode {
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
pub enum RangeUnit {
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
pub enum AngleType {
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
pub enum Interpolation {
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
pub enum DataQuality {
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
pub enum YesNo {
    Yes,
    No,
}

impl Default for YesNo {
    fn default() -> Self {
        Self::Yes
    }
}
