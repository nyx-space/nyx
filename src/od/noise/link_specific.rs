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

use super::{StochasticNoise, WhiteNoise};
use anise::constants::SPEED_OF_LIGHT_KM_S;
use hifitime::Duration;
use serde_derive::{Deserialize, Serialize};
use std::f64::consts::TAU;

/// Signal power to noise density for stochastic modeling, typical values.
/// IMPORTANT: The S/N0 will always be lower or equal to the Carrier Power to noise density (C/N0)
#[derive(Copy, Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum SN0 {
    /// 65 dB-Hz
    Strong,
    /// 50 dB-Hz
    #[default]
    Average,
    /// 40 dB-Hz
    Poor,
    /// Manual value provided in dB-Hz, converted to Hertz automatically
    ManualDbHz(f64),
}

impl SN0 {
    /// Note that this returns the data in Hertz not dB-Hz
    pub(crate) fn value_hz(self) -> f64 {
        match self {
            Self::Strong => 10.0_f64.powf(6.5),
            Self::Average => 10.0_f64.powi(5),
            Self::Poor => 10.0_f64.powi(4),
            Self::ManualDbHz(value) => 10.0_f64.powf(value / 10.0),
        }
    }
}

/// Carrier power to noise density for stochastic modeling, typical values.
/// IMPORTANT: The C/N0 will always be greater or equal to the Signal Power to noise density (C/N0) because the S/N0 for a ranging
/// tone is the power dedicated to the ranging within the uplink (cf. "modulation index" in the DESCANSO monograph).
#[derive(Copy, Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum CN0 {
    /// 70 dB-Hz
    Strong,
    /// 55 dB-Hz
    #[default]
    Average,
    /// 45 dB-Hz
    Poor,
    /// Manual value provided in dB-Hz, converted to Hertz automatically
    ManualDbHz(f64),
}

impl CN0 {
    /// Note that this returns the data in Hertz not dB-Hz
    pub(crate) fn value_hz(self) -> f64 {
        match self {
            Self::Strong => 10.0_f64.powi(7),
            Self::Average => 10.0_f64.powf(5.5),
            Self::Poor => 10.0_f64.powf(4.5),
            Self::ManualDbHz(value) => 10.0_f64.powf(value / 10.0),
        }
    }
}

/// Carrier frequency helper enum, typical values.
pub enum CarrierFreq {
    /// 2.2 GHz
    SBand,
    /// 8.4 GHz
    XBand,
    /// 32 Ghz
    KaBand,
    ManualHz(f64),
}

impl CarrierFreq {
    pub(crate) fn value_hz(self) -> f64 {
        match self {
            Self::SBand => 2.2e9,
            Self::XBand => 8.4e9,
            Self::KaBand => 32e9,
            Self::ManualHz(value) => value,
        }
    }
}

/// An enum helper with typical chip rates.
#[derive(Copy, Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub enum ChipRate {
    /// 1 kchip/s -- basically emergency ranging
    Lowest,
    /// 100 kchip/s -- could be used for weaker links
    Low,
    /// 1 Mchip/s -- typical of xGEO/cislunar missions
    #[default]
    StandardT4B,
    /// 10 Mchip/s -- high-precision scientific missions (e.g. gravity modeling)
    High,
    /// 25 Mchip/s -- highly specialized missions
    VeryHigh,
    /// Provide your own chip rate depending on the ground station configuration
    ManualHz(f64),
}

impl ChipRate {
    pub(crate) fn value_chip_s(self) -> f64 {
        match self {
            Self::Lowest => 1e3,
            Self::Low => 1e5,
            Self::StandardT4B => 1e6,
            Self::High => 1e7,
            Self::VeryHigh => 2.5e7,
            Self::ManualHz(value) => value,
        }
    }
}

impl StochasticNoise {
    /// Constructs a high precision zero-mean range noise model (accounting for clock error and thermal error) from
    /// the Allan deviation of the clock, integration time, chip rate (depends on the ranging code), and
    /// signal-power-to-noise-density ratio (S/N₀).
    ///
    /// NOTE: The Allan Deviation should be provided given the integration time. For example, if the integration time
    /// is one second, the Allan Deviation should be the deviation over one second.
    ///
    /// IMPORTANT: These do NOT include atmospheric noises, which add up to ~10 cm one-sigma.
    pub fn from_hardware_range_km(
        allan_deviation: f64,
        integration_time: Duration,
        chip_rate: ChipRate,
        s_n0: SN0,
    ) -> Self {
        // Compute the thermal noise.
        let sigma_thermal_km =
            SPEED_OF_LIGHT_KM_S / (TAU * chip_rate.value_chip_s() * (2.0 * s_n0.value_hz()).sqrt());
        // Compute the clock noise.
        let sigma_clock_km =
            (SPEED_OF_LIGHT_KM_S * allan_deviation * integration_time.to_seconds())
                / (3.0_f64.sqrt());

        Self {
            white_noise: Some(WhiteNoise::constant_white_noise(
                (sigma_clock_km.powi(2) + sigma_thermal_km.powi(2)).sqrt(),
            )),
            bias: None,
        }
    }

    pub fn from_hardware_doppler_km_s(
        allan_deviation: f64,
        integration_time: Duration,
        carrier: CarrierFreq,
        c_n0: CN0,
    ) -> Self {
        // Compute the thermal noise
        let sigma_thermal_km_s = SPEED_OF_LIGHT_KM_S
            / (TAU
                * carrier.value_hz()
                * (2.0 * c_n0.value_hz() * integration_time.to_seconds()).sqrt());

        // Compute the clock noise.
        let sigma_clock_km_s = SPEED_OF_LIGHT_KM_S * allan_deviation;

        Self {
            white_noise: Some(WhiteNoise::constant_white_noise(
                (sigma_clock_km_s.powi(2) + sigma_thermal_km_s.powi(2)).sqrt(),
            )),
            bias: None,
        }
    }
}

#[cfg(test)]
mod link_noise {
    use super::{CarrierFreq, ChipRate, StochasticNoise, CN0, SN0};
    use hifitime::Unit;
    #[test]
    fn nasa_dsac() {
        // The DSAC has an Allan Dev of 1e-14 over one day.
        // Gemini claims that such a good clock likely has the same deviation over 60 seconds (hitting the flicker floor).
        // But worst case scenario, its AD is the square root of the ratio of one day over 1 minute, or 38 times worse.

        for (case_num, allan_dev) in [1e-14, 3.8e-13].iter().copied().enumerate() {
            println!("AD = {allan_dev:e}");

            let range_dsac_no_flicker = StochasticNoise::from_hardware_range_km(
                allan_dev,
                Unit::Minute * 1,
                ChipRate::StandardT4B,
                SN0::Average,
            );

            let range_sigma_m = range_dsac_no_flicker.white_noise.unwrap().sigma * 1e3;

            println!("range sigma = {range_sigma_m:.3e} m",);

            assert!(range_sigma_m.abs() < 1.1e-1);

            let doppler_dsac_no_flicker = StochasticNoise::from_hardware_doppler_km_s(
                allan_dev,
                Unit::Minute * 1,
                CarrierFreq::XBand,
                CN0::Average,
            );

            let doppler_sigma_m_s = doppler_dsac_no_flicker.white_noise.unwrap().sigma * 1e3;

            println!("doppler sigma = {doppler_sigma_m_s:.3e} m/s");

            match case_num {
                0 => assert!(doppler_sigma_m_s < 3.2e-6),
                1 => assert!(doppler_sigma_m_s < 1.2e-4),
                _ => unreachable!(),
            };
        }
    }
}
