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

use crate::io::watermark::pq_writer;
use arrow::array::{ArrayRef, Float64Array, UInt32Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use der::{Decode, Encode, Reader};
use hifitime::{Epoch, TimeSeries, TimeUnits};
use parquet::arrow::ArrowWriter;

use rand::rngs::SysRng;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;
use serde::{Deserialize, Serialize};
use std::error::Error;
use std::fmt::Display;
use std::fs::File;
use std::ops::{Mul, MulAssign};
use std::path::Path;
use std::sync::Arc;

pub mod gauss_markov;
pub mod link_specific;
pub mod white;

#[cfg(feature = "python")]
use hifitime::Duration;
#[cfg(feature = "python")]
use pyo3::exceptions::PyValueError;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::PyType;

pub use gauss_markov::GaussMarkov;
pub use white::WhiteNoise;

/// Trait for any kind of stochastic modeling, developing primarily for synthetic orbit determination measurements.
pub trait Stochastics {
    /// Return the variance of this stochastic noise model at a given time.
    fn covariance(&self, epoch: Epoch) -> f64;

    /// Returns a new sample of these stochastics
    fn sample<R: Rng>(&mut self, epoch: Epoch, rng: &mut R) -> f64;
}

/// Stochastic noise modeling used primarily for synthetic orbit determination measurements.
///
/// This implementation distinguishes between the white noise model and the bias model. It also includes a constant offset.
#[derive(Copy, Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
#[cfg_attr(feature = "python", pyclass(from_py_object, get_all, set_all))]
pub struct StochasticNoise {
    pub white_noise: Option<WhiteNoise>,
    pub bias: Option<GaussMarkov>,
}

impl StochasticNoise {
    /// Zero noise stochastic process.
    pub const ZERO: Self = Self {
        white_noise: None,
        bias: None,
    };

    /// The minimum stochastic noise process with a zero mean white noise of 1e-6.
    pub const MIN: Self = Self {
        white_noise: Some(WhiteNoise {
            mean: 0.0,
            sigma: 1e-6,
        }),
        bias: None,
    };

    /// Default stochastic process of the Deep Space Network, as per DESCANSO Chapter 3, Table 3-3.
    /// Using the instrument bias as the white noise value, zero constant bias.
    pub fn default_range_km() -> Self {
        Self {
            white_noise: Some(WhiteNoise {
                sigma: 2.0e-3, // 2 m
                ..Default::default()
            }),
            // Until Nyx can support bias estimation the default bias is None, cf. https://github.com/nyx-space/nyx/issues/326
            // bias: Some(GaussMarkov::default_range_km()),
            bias: None,
        }
    }

    /// Default stochastic process of the Deep Space Network, using as per DESCANSO Chapter 3, Table 3-3 for the GM process.
    pub fn default_doppler_km_s() -> Self {
        Self {
            white_noise: Some(WhiteNoise {
                sigma: 3e-6, // 3 mm/s
                ..Default::default()
            }),
            // Until Nyx can support bias estimation the default bias is None, cf. https://github.com/nyx-space/nyx/issues/326
            // bias: Some(GaussMarkov::default_doppler_km_s()),
            bias: None,
        }
    }

    /// Default stochastic process for an angle measurement (azimuth or elevation)
    /// Using the instrument bias as the white noise value, zero constant bias.
    pub fn default_angle_deg() -> Self {
        Self {
            white_noise: Some(WhiteNoise {
                sigma: 1.0e-2, // 0.01 deg
                ..Default::default()
            }),
            // Until Nyx can support bias estimation the default bias is None, cf. https://github.com/nyx-space/nyx/issues/326
            // bias: Some(GaussMarkov::default_range_km()),
            bias: None,
        }
    }

    /// Sample these stochastics
    pub fn sample<R: Rng>(&mut self, epoch: Epoch, rng: &mut R) -> f64 {
        let mut sample = 0.0;
        if let Some(wn) = &mut self.white_noise {
            sample += wn.sample(epoch, rng)
        }
        if let Some(gm) = &mut self.bias {
            sample += gm.sample(epoch, rng);
        }
        sample
    }

    /// Executes a hardcoded 24-hour Monte Carlo simulation of the stochastic model, exporting the time history to a Parquet file.
    ///
    /// # Warning: Hardcoded Time Series & Diagnostic Data Gaps
    /// This method does *not* accept a user-defined tracking schedule or time series. It inherently evaluates the stochastic process
    /// over a strict 24-hour period, beginning at the exact system clock moment of method execution, utilizing a 1-minute step size.
    ///
    /// Furthermore, users will observe exactly 1,082 samples per simulation run, rather than the 1,441 samples expected from a
    /// continuous 24-hour 1-minute cadence. The simulation intentionally drops all epochs strictly greater than +6 hours and
    /// strictly less than +12 hours from the start time. This hardcoded artifact is designed to demonstrate variance bounds
    /// expansion in the absence of measurements (e.g., simulating a tracking dropout for a Gauss-Markov bias).
    ///
    /// # Algorithm
    /// 1. Establish `start` as the system clock time at invocation.
    /// 2. Construct an inclusive time series from `start` to `start + 24 hours` at 1-minute intervals.
    /// 3. For each configured run, seed a PRNG (`Pcg64Mcg`) using system entropy.
    /// 4. Evaluate the process covariance and sample the stochastic noise at each epoch.
    /// 5. Discard all epochs inside the `(start + 6h, start + 12h)` open interval.
    /// 6. Export the remaining 1,082 samples per run to an Apache Arrow RecordBatch and write to disk via Parquet.
    pub fn simulate<P: AsRef<Path>>(
        self,
        path: P,
        runs: Option<u32>,
        unit: Option<String>,
    ) -> Result<Vec<StochasticState>, Box<dyn Error>> {
        let num_runs = runs.unwrap_or(25);

        let start = Epoch::now().unwrap();
        let (step, end) = (1.minutes(), start + 1.days());

        let capacity = ((end - start).to_seconds() / step.to_seconds()).ceil() as usize;

        let mut samples = Vec::with_capacity(capacity);

        for run in 0..num_runs {
            let mut rng = Pcg64Mcg::try_from_rng(&mut SysRng).unwrap();

            let mut mdl = self;
            for epoch in TimeSeries::inclusive(start, end, step) {
                if epoch > start + 6.hours() && epoch < start + 12.hours() {
                    // Skip to see how the variance changes.
                    continue;
                }
                let variance = mdl.covariance(epoch);
                let sample = mdl.sample(epoch, &mut rng);
                samples.push(StochasticState {
                    run,
                    dt_s: (epoch - start).to_seconds(),
                    sample,
                    variance,
                });
            }
        }

        let bias_unit = match unit {
            Some(unit) => format!("({unit})"),
            None => "(unitless)".to_string(),
        };

        // Build the parquet file
        let hdrs = vec![
            Field::new("Run", DataType::UInt32, false),
            Field::new("Delta Time (s)", DataType::Float64, false),
            Field::new(format!("Bias {bias_unit}"), DataType::Float64, false),
            Field::new(format!("Variance {bias_unit}"), DataType::Float64, false),
        ];

        let schema = Arc::new(Schema::new(hdrs));
        let record = vec![
            Arc::new(UInt32Array::from(
                samples.iter().map(|s| s.run).collect::<Vec<u32>>(),
            )) as ArrayRef,
            Arc::new(Float64Array::from(
                samples.iter().map(|s| s.dt_s).collect::<Vec<f64>>(),
            )) as ArrayRef,
            Arc::new(Float64Array::from(
                samples.iter().map(|s| s.sample).collect::<Vec<f64>>(),
            )) as ArrayRef,
            Arc::new(Float64Array::from(
                samples.iter().map(|s| s.variance).collect::<Vec<f64>>(),
            )) as ArrayRef,
        ];

        let props = pq_writer(None);

        let file = File::create(path)?;
        let mut writer = ArrowWriter::try_new(file, schema.clone(), props).unwrap();

        let batch = RecordBatch::try_new(schema, record)?;
        writer.write(&batch)?;
        writer.close()?;

        Ok(samples)
    }

    fn available_data(&self) -> u8 {
        let mut bits: u8 = 0;

        if self.white_noise.is_some() {
            bits |= 1 << 0;
        }

        if self.bias.is_some() {
            bits |= 1 << 1;
        }

        bits
    }
}

#[cfg_attr(feature = "python", pymethods)]
impl StochasticNoise {
    #[cfg(feature = "python")]
    #[pyo3(signature=(white_noise=None, bias=None, name=None))]
    #[new]
    fn py_new(
        white_noise: Option<WhiteNoise>,
        bias: Option<GaussMarkov>,
        name: Option<String>,
    ) -> PyResult<Self> {
        if let Some(name) = name {
            match name.to_ascii_lowercase().as_str() {
                "range" => Ok(Self::default_range_km()),
                "doppler" => Ok(Self::default_doppler_km_s()),
                "angles" => Ok(Self::default_angle_deg()),
                _ => Err(PyValueError::new_err(format!(
                    "name must be `range`, `doppler`, or `angles` (received `{name}`)"
                ))),
            }
        } else {
            Ok(Self { white_noise, bias })
        }
    }

    /// Return the covariance of these stochastics at a given time.
    pub fn covariance(&self, epoch: Epoch) -> f64 {
        let mut variance = 0.0;
        if let Some(wn) = &self.white_noise {
            variance += wn.covariance(epoch);
        }
        if let Some(gm) = &self.bias {
            variance += gm.covariance(epoch);
        }
        variance
    }

    /// Executes a hardcoded 24-hour Monte Carlo simulation of the stochastic model, exporting the time history to a Parquet file.
    ///
    /// # Warning: Hardcoded Time Series & Diagnostic Data Gaps
    /// This method does *not* accept a user-defined tracking schedule or time series. It inherently evaluates the stochastic process
    /// over a strict 24-hour period, beginning at the exact system clock moment of method execution, utilizing a 1-minute step size.
    ///
    /// Furthermore, users will observe exactly 1,082 samples per simulation run, rather than the 1,441 samples expected from a
    /// continuous 24-hour 1-minute cadence. The simulation intentionally drops all epochs strictly greater than +6 hours and
    /// strictly less than +12 hours from the start time. This hardcoded artifact is designed to demonstrate variance bounds
    /// expansion in the absence of measurements (e.g., simulating a tracking dropout for a Gauss-Markov bias).
    ///
    /// # Algorithm
    /// 1. Establish `start` as the system clock time at invocation.
    /// 2. Construct an inclusive time series from `start` to `start + 24 hours` at 1-minute intervals.
    /// 3. For each configured run, seed a PRNG (`Pcg64Mcg`) using system entropy.
    /// 4. Evaluate the process covariance and sample the stochastic noise at each epoch.
    /// 5. Discard all epochs inside the `(start + 6h, start + 12h)` open interval.
    /// 6. Export the remaining 1,082 samples per run to an Apache Arrow RecordBatch and write to disk via Parquet.
    ///
    /// :param path: The filesystem path for the output Parquet file.
    /// :type path: str
    /// :param runs: The number of Monte Carlo runs. Defaults to 25 if not provided.
    /// :type runs: int | None
    /// :param unit: An optional string appended to the Parquet column headers for plotting clarity.
    /// :type unit: str | None
    /// :rtype: list[StochasticState]
    /// :raises Exception: If the underlying Apache Arrow RecordBatch fails to allocate or write to the specified filesystem path.
    #[cfg(feature = "python")]
    #[pyo3(name = "simulate")]
    fn py_simulate(
        &self,
        path: &str,
        runs: Option<u32>,
        unit: Option<String>,
    ) -> PyResult<Vec<StochasticState>> {
        self.simulate(path, runs, unit)
            .map_err(|e| PyValueError::new_err(e.to_string()))
    }

    /// Constructs a high precision zero-mean range noise model (accounting for clock error and thermal error) from
    /// the Allan deviation of the clock, integration time, chip rate (depends on the ranging code), and
    /// signal-power-to-noise-density ratio (S/N₀).
    ///
    /// NOTE: The Allan Deviation should be provided given the integration time. For example, if the integration time
    /// is one second, the Allan Deviation should be the deviation over one second.
    ///
    /// IMPORTANT: These do NOT include atmospheric noises, which add up to ~10 cm one-sigma.
    #[cfg(feature = "python")]
    #[pyo3(name = "from_hardware_range_km")]
    #[classmethod]
    fn py_from_hardware_range_km(
        _cls: &Bound<'_, PyType>,

        allan_deviation: f64,
        integration_time: Duration,
        chip_rate: link_specific::ChipRate,
        s_n0: link_specific::SN0,
    ) -> Self {
        Self::from_hardware_range_km(allan_deviation, integration_time, chip_rate, s_n0)
    }

    #[cfg(feature = "python")]
    #[pyo3(name = "from_hardware_doppler_km_s")]
    #[classmethod]
    fn py_from_hardware_doppler_km_s(
        _cls: &Bound<'_, PyType>,
        allan_deviation: f64,
        integration_time: Duration,
        carrier: link_specific::CarrierFreq,
        c_n0: link_specific::CN0,
    ) -> Self {
        Self::from_hardware_doppler_km_s(allan_deviation, integration_time, carrier, c_n0)
    }

    fn __str__(&self) -> String {
        format!("{self}")
    }

    fn __repr__(&self) -> String {
        format!("{self} @ {self:p}")
    }
}

impl Display for StochasticNoise {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match (self.white_noise, self.bias) {
            (Some(wn), None) => write!(f, "Stochastics with {wn:?}"),
            (None, Some(bias)) => write!(f, "Stochastics with bias {bias}"),
            (None, None) => write!(f, "Noiseless stochastics"),
            (Some(wn), Some(bias)) => write!(f, "Stochastics with {wn:?} and bias {bias}"),
        }
    }
}

impl Mul<f64> for StochasticNoise {
    type Output = Self;

    fn mul(mut self, rhs: f64) -> Self::Output {
        if let Some(wn) = &mut self.white_noise {
            *wn *= rhs;
        }
        if let Some(gm) = &mut self.bias {
            *gm *= rhs;
        }

        self
    }
}

impl MulAssign<f64> for StochasticNoise {
    fn mul_assign(&mut self, rhs: f64) {
        *self = *self * rhs;
    }
}

impl Encode for StochasticNoise {
    fn encoded_len(&self) -> der::Result<der::Length> {
        let flags = self.available_data();
        flags.encoded_len()? + self.white_noise.encoded_len()? + self.bias.encoded_len()?
    }

    fn encode(&self, encoder: &mut impl der::Writer) -> der::Result<()> {
        let flags = self.available_data();

        flags.encode(encoder)?;
        self.white_noise.encode(encoder)?;
        self.bias.encode(encoder)
    }
}

impl<'a> Decode<'a> for StochasticNoise {
    fn decode<R: Reader<'a>>(decoder: &mut R) -> der::Result<Self> {
        let flags: u8 = decoder.decode()?;

        let white_noise = if flags & (1 << 0) != 0 {
            Some(decoder.decode()?)
        } else {
            None
        };

        let bias = if flags & (1 << 1) != 0 {
            Some(decoder.decode()?)
        } else {
            None
        };

        Ok(Self { white_noise, bias })
    }
}

#[derive(Copy, Clone, Debug)]
#[cfg_attr(feature = "python", pyclass(from_py_object, get_all))]
pub struct StochasticState {
    pub run: u32,
    pub dt_s: f64,
    pub sample: f64,
    pub variance: f64,
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl StochasticState {
    fn __str__(&self) -> String {
        format!("{self:?}")
    }
    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}

#[cfg(test)]
mod ut_stochastics {
    use std::path::PathBuf;

    use super::{StochasticNoise, white::WhiteNoise};

    #[test]
    fn test_simulate_zero() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "../data",
            "04_output",
            "stochastics_zero.parquet",
        ]
        .iter()
        .collect();

        let noise = StochasticNoise::default();

        let rslts = noise.simulate(path, None, None).unwrap();
        assert!(!rslts.is_empty());
        assert!(rslts.iter().map(|rslt| rslt.sample).sum::<f64>().abs() < f64::EPSILON);
    }

    #[test]
    fn test_simulate_constant() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "../data",
            "04_output",
            "stochastics_constant.parquet",
        ]
        .iter()
        .collect();

        let noise = StochasticNoise {
            white_noise: Some(WhiteNoise {
                mean: 15.0,
                sigma: 2.0,
            }),
            ..Default::default()
        };

        noise.simulate(path, None, None).unwrap();
    }

    #[test]
    fn test_simulate_dsn_range() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "../data",
            "04_output",
            "stochastics_dsn_range.parquet",
        ]
        .iter()
        .collect();

        let noise = StochasticNoise::default_range_km();

        noise
            .simulate(path, None, Some("kilometer".to_string()))
            .unwrap();
    }

    #[test]
    fn test_simulate_dsn_range_gm_only() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "../data",
            "04_output",
            "stochastics_dsn_range_gm_only.parquet",
        ]
        .iter()
        .collect();

        let mut noise = StochasticNoise::default_range_km();
        noise.white_noise = None;

        noise
            .simulate(path, None, Some("kilometer".to_string()))
            .unwrap();
    }
}
