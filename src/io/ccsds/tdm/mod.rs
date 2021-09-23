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
use crate::{time::Epoch, NyxError};
use std::convert::TryFrom;
use std::str::FromStr;
use yaserde_derive::{YaDeserialize, YaSerialize};

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
    #[yaserde(child)]
    pub data: Data,
}

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
pub struct Data {
    #[yaserde(child, rename = "COMMENT")]
    comment: Vec<String>,
    #[yaserde(child, rename = "observation")]
    pub observations: Vec<TdmObservation>,
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

impl Data {
    pub fn comments(&self) -> String {
        self.comment.join(" ")
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

#[derive(Default, Debug, Deserialize, PartialEq, YaSerialize, YaDeserialize)]
#[yaserde(flatten)]
#[allow(non_snake_case)]
pub struct TdmObservation {
    #[yaserde(child, rename = "EPOCH")]
    epoch: String,
    #[yaserde(child, rename = "ANGLE_1")]
    angle_1_deg: Option<f64>,
    #[yaserde(child, rename = "ANGLE_2")]
    angle_2_deg: Option<f64>,
    #[yaserde(child, rename = "CARRIER_POWER")]
    carrier_power_dbw: Option<f64>,
    #[yaserde(child, rename = "CLOCK_BIAS")]
    clock_bias_s: Option<f64>,
    #[yaserde(child, rename = "CLOCK_DRIFT")]
    clock_drift_ss: Option<f64>,
    #[yaserde(child, rename = "COMMENT")]
    comment: Option<String>,
    #[yaserde(child, rename = "DOPPLER_COUNT")]
    doppler_count: Option<f64>,
    #[yaserde(child, rename = "DOPPLER_INSTANTANEOUS")]
    doppler_instantaneous_kms: Option<f64>,
    #[yaserde(child, rename = "DOPPLER_INTEGRATED")]
    doppler_integrated_kms: Option<f64>,
    #[yaserde(child, rename = "DOR")]
    dor_s: Option<f64>,
    #[yaserde(child, rename = "MAG")]
    mag: Option<f64>,
    /// Carrier power to noise spectral density ratio
    #[yaserde(child, rename = "PC_N0")]
    pc_n0_dbhz: Option<f64>,
    /// Ranging power to noise spectral density ratio
    #[yaserde(child, rename = "PR_N0")]
    pr_n0_dbhz: Option<f64>,
    #[yaserde(child, rename = "PRESSURE")]
    pressure_hpa: Option<f64>,
    /// Range is specified in either km, s, or RU (hence the "xyz" notation for the units)
    #[yaserde(child, rename = "RANGE")]
    range_xyz: Option<f64>,
    /// Radar cross section
    #[yaserde(child, rename = "RCS")]
    rcs_m2: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_FREQ_1")]
    receive_freq_1_hz: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_FREQ_2")]
    receive_freq_2_hz: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_FREQ_3")]
    receive_freq_3_hz: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_FREQ_4")]
    receive_freq_4_hz: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_FREQ_5")]
    receive_freq_5_hz: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_FREQ")]
    receive_freq_hz: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_PHASE_CT_1")]
    receive_phase_ct_1: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_PHASE_CT_2")]
    receive_phase_ct_2: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_PHASE_CT_3")]
    receive_phase_ct_3: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_PHASE_CT_4")]
    receive_phase_ct_4: Option<f64>,
    #[yaserde(child, rename = "RECEIVE_PHASE_CT_5")]
    receive_phase_ct_5: Option<f64>,
    #[yaserde(child, rename = "RHUMIDITY")]
    rhumidity_prct: Option<f64>,
    /// Slant Total Electron Count
    #[yaserde(child, rename = "STEC")]
    stec_tecu: Option<f64>,
    #[yaserde(child, rename = "TEMPERATURE")]
    temperature_k: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_FREQ_1")]
    transmit_freq_1_hz: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_FREQ_2")]
    transmit_freq_2_hz: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_FREQ_3")]
    transmit_freq_3_hz: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_FREQ_4")]
    transmit_freq_4_hz: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_FREQ_5")]
    transmit_freq_5_hz: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_FREQ_RATE_1")]
    transmit_freq_rate_1_hzs: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_FREQ_RATE_2")]
    transmit_freq_rate_2_hzs: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_FREQ_RATE_3")]
    transmit_freq_rate_3_hzs: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_FREQ_RATE_4")]
    transmit_freq_rate_4_hzs: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_FREQ_RATE_5")]
    transmit_freq_rate_5_hzs: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_PHASE_CT_1")]
    transmit_phase_ct_1: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_PHASE_CT_2")]
    transmit_phase_ct_2: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_PHASE_CT_3")]
    transmit_phase_ct_3: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_PHASE_CT_4")]
    transmit_phase_ct_4: Option<f64>,
    #[yaserde(child, rename = "TRANSMIT_PHASE_CT_5")]
    transmit_phase_ct_5: Option<f64>,
    /// Dry zenith delay through the  troposphere  measured  at  the  timetag
    #[yaserde(child, rename = "TROPO_DRY")]
    tropo_dry_m: Option<f64>,
    /// Wet zenith delay through the  troposphere  measured  at  the  timetag
    #[yaserde(child, rename = "TROPO_WET")]
    tropo_wet_m: Option<f64>,
    #[yaserde(child, rename = "VLBI_DELAY")]
    vlbi_delay_s: Option<f64>,
}

/// A parsed observation value from a TDM
#[derive(Debug, Copy, Clone, PartialEq)]
#[allow(non_camel_case_types)]
pub enum ObservationValue {
    Angle1_deg(f64),
    Angle2_deg(f64),
    CarrierPower_dbw(f64),
    ClockBias_s(f64),
    ClockDrift_ss(f64),
    DopplerCount(f64),
    DopplerInstantaneous_kms(f64),
    DopplerIntegrated_kms(f64),
    Dor_s(f64),
    Mag(f64),
    /// Carrier power to noise spectral density ratio
    PcN0_dbhz(f64),
    /// Ranging power to noise spectral density ratio
    PrN0_dbhz(f64),
    Pressure_hpa(f64),
    /// Range is specified in either km, s, or RU (hence the "xyz" notation for the units)
    Range_xyz(f64),
    /// Radar cross section
    Rcs_m2(f64),
    ReceiveFreq1_hz(f64),
    ReceiveFreq2_hz(f64),
    ReceiveFreq3_hz(f64),
    ReceiveFreq4_hz(f64),
    ReceiveFreq5_hz(f64),
    ReceiveFreq_hz(f64),
    ReceivePhase_ct_1(f64),
    ReceivePhase_ct_2(f64),
    ReceivePhase_ct_3(f64),
    ReceivePhase_ct_4(f64),
    ReceivePhase_ct_5(f64),
    Rhumidity_prct(f64),
    Stec_tecu(f64),
    /// Slant Total Electron Count
    Temperature_k(f64),
    TransmitFreq1_hz(f64),
    TransmitFreq2_hz(f64),
    TransmitFreq3_hz(f64),
    TransmitFreq4_hz(f64),
    TransmitFreq5_hz(f64),
    TransmitFreqRate1_hzs(f64),
    TransmitFreqRate2_hzs(f64),
    TransmitFreqRate3_hzs(f64),
    TransmitFreqRate4_hzs(f64),
    TransmitFreqRate5_hzs(f64),
    TransmitPhaseCt_1(f64),
    TransmitPhaseCt_2(f64),
    TransmitPhaseCt_3(f64),
    TransmitPhaseCt_4(f64),
    TransmitPhaseCt_5(f64),
    TropoDry_m(f64),
    /// Dry zenith delay through the  troposphere  measured  at  the  timetag
    TropoWet_m(f64),
    /// Wet zenith delay through the  troposphere  measured  at  the  timetag
    VlbiDelay_s(f64),
}

/// A parsed observation from a TDM, with its epoch and value
#[derive(Debug)]
pub struct Observation {
    pub epoch: Epoch,
    pub value: ObservationValue,
}

impl TryFrom<&TdmObservation> for Observation {
    type Error = NyxError;

    fn try_from(raw_obs: &TdmObservation) -> Result<Self, Self::Error> {
        let epoch = match Epoch::from_str(&raw_obs.epoch) {
            Ok(e) => e,
            Err(_) => return Err(NyxError::CCSDS("Could not decode epoch".to_string())),
        };

        let value = if raw_obs.angle_1_deg.is_some() {
            ObservationValue::Angle1_deg(raw_obs.angle_1_deg.unwrap())
        } else if raw_obs.angle_2_deg.is_some() {
            ObservationValue::Angle2_deg(raw_obs.angle_2_deg.unwrap())
        } else if raw_obs.carrier_power_dbw.is_some() {
            ObservationValue::CarrierPower_dbw(raw_obs.carrier_power_dbw.unwrap())
        } else if raw_obs.clock_bias_s.is_some() {
            ObservationValue::ClockBias_s(raw_obs.clock_bias_s.unwrap())
        } else if raw_obs.clock_drift_ss.is_some() {
            ObservationValue::ClockDrift_ss(raw_obs.clock_drift_ss.unwrap())
        } else if raw_obs.doppler_count.is_some() {
            ObservationValue::DopplerCount(raw_obs.doppler_count.unwrap())
        } else if raw_obs.doppler_instantaneous_kms.is_some() {
            ObservationValue::DopplerInstantaneous_kms(raw_obs.doppler_instantaneous_kms.unwrap())
        } else if raw_obs.doppler_integrated_kms.is_some() {
            ObservationValue::DopplerIntegrated_kms(raw_obs.doppler_integrated_kms.unwrap())
        } else if raw_obs.dor_s.is_some() {
            ObservationValue::Dor_s(raw_obs.dor_s.unwrap())
        } else if raw_obs.mag.is_some() {
            ObservationValue::Mag(raw_obs.mag.unwrap())
        } else if raw_obs.pc_n0_dbhz.is_some() {
            ObservationValue::PcN0_dbhz(raw_obs.pc_n0_dbhz.unwrap())
        } else if raw_obs.pr_n0_dbhz.is_some() {
            ObservationValue::PrN0_dbhz(raw_obs.pr_n0_dbhz.unwrap())
        } else if raw_obs.pressure_hpa.is_some() {
            ObservationValue::Pressure_hpa(raw_obs.pressure_hpa.unwrap())
        } else if raw_obs.range_xyz.is_some() {
            ObservationValue::Range_xyz(raw_obs.range_xyz.unwrap())
        } else if raw_obs.rcs_m2.is_some() {
            ObservationValue::Rcs_m2(raw_obs.rcs_m2.unwrap())
        } else if raw_obs.receive_freq_1_hz.is_some() {
            ObservationValue::ReceiveFreq1_hz(raw_obs.receive_freq_1_hz.unwrap())
        } else if raw_obs.receive_freq_2_hz.is_some() {
            ObservationValue::ReceiveFreq2_hz(raw_obs.receive_freq_2_hz.unwrap())
        } else if raw_obs.receive_freq_3_hz.is_some() {
            ObservationValue::ReceiveFreq3_hz(raw_obs.receive_freq_3_hz.unwrap())
        } else if raw_obs.receive_freq_4_hz.is_some() {
            ObservationValue::ReceiveFreq4_hz(raw_obs.receive_freq_4_hz.unwrap())
        } else if raw_obs.receive_freq_5_hz.is_some() {
            ObservationValue::ReceiveFreq5_hz(raw_obs.receive_freq_5_hz.unwrap())
        } else if raw_obs.receive_freq_hz.is_some() {
            ObservationValue::ReceiveFreq_hz(raw_obs.receive_freq_hz.unwrap())
        } else if raw_obs.rhumidity_prct.is_some() {
            ObservationValue::Rhumidity_prct(raw_obs.rhumidity_prct.unwrap())
        } else if raw_obs.stec_tecu.is_some() {
            ObservationValue::Stec_tecu(raw_obs.stec_tecu.unwrap())
        } else if raw_obs.temperature_k.is_some() {
            ObservationValue::Temperature_k(raw_obs.temperature_k.unwrap())
        } else if raw_obs.transmit_freq_1_hz.is_some() {
            ObservationValue::TransmitFreq1_hz(raw_obs.transmit_freq_1_hz.unwrap())
        } else if raw_obs.transmit_freq_2_hz.is_some() {
            ObservationValue::TransmitFreq2_hz(raw_obs.transmit_freq_2_hz.unwrap())
        } else if raw_obs.transmit_freq_3_hz.is_some() {
            ObservationValue::TransmitFreq3_hz(raw_obs.transmit_freq_3_hz.unwrap())
        } else if raw_obs.transmit_freq_4_hz.is_some() {
            ObservationValue::TransmitFreq4_hz(raw_obs.transmit_freq_4_hz.unwrap())
        } else if raw_obs.transmit_freq_5_hz.is_some() {
            ObservationValue::TransmitFreq5_hz(raw_obs.transmit_freq_5_hz.unwrap())
        } else if raw_obs.transmit_phase_ct_1.is_some() {
            ObservationValue::TransmitPhaseCt_1(raw_obs.transmit_phase_ct_1.unwrap())
        } else if raw_obs.transmit_phase_ct_2.is_some() {
            ObservationValue::TransmitPhaseCt_2(raw_obs.transmit_phase_ct_2.unwrap())
        } else if raw_obs.transmit_phase_ct_3.is_some() {
            ObservationValue::TransmitPhaseCt_3(raw_obs.transmit_phase_ct_3.unwrap())
        } else if raw_obs.transmit_phase_ct_4.is_some() {
            ObservationValue::TransmitPhaseCt_4(raw_obs.transmit_phase_ct_4.unwrap())
        } else if raw_obs.transmit_phase_ct_5.is_some() {
            ObservationValue::TransmitPhaseCt_5(raw_obs.transmit_phase_ct_5.unwrap())
        } else if raw_obs.tropo_dry_m.is_some() {
            ObservationValue::TropoDry_m(raw_obs.tropo_dry_m.unwrap())
        } else if raw_obs.tropo_wet_m.is_some() {
            ObservationValue::TropoWet_m(raw_obs.tropo_wet_m.unwrap())
        } else if raw_obs.vlbi_delay_s.is_some() {
            ObservationValue::VlbiDelay_s(raw_obs.vlbi_delay_s.unwrap())
        } else if raw_obs.transmit_freq_rate_1_hzs.is_some() {
            ObservationValue::TransmitFreqRate1_hzs(raw_obs.transmit_freq_rate_1_hzs.unwrap())
        } else if raw_obs.transmit_freq_rate_2_hzs.is_some() {
            ObservationValue::TransmitFreqRate2_hzs(raw_obs.transmit_freq_rate_2_hzs.unwrap())
        } else if raw_obs.transmit_freq_rate_3_hzs.is_some() {
            ObservationValue::TransmitFreqRate3_hzs(raw_obs.transmit_freq_rate_3_hzs.unwrap())
        } else if raw_obs.transmit_freq_rate_4_hzs.is_some() {
            ObservationValue::TransmitFreqRate4_hzs(raw_obs.transmit_freq_rate_4_hzs.unwrap())
        } else if raw_obs.transmit_freq_rate_5_hzs.is_some() {
            ObservationValue::TransmitFreqRate5_hzs(raw_obs.transmit_freq_rate_5_hzs.unwrap())
        } else {
            return Err(NyxError::CCSDS(
                "No observation value specified".to_string(),
            ));
        };

        Ok(Observation { epoch, value })
    }
}
