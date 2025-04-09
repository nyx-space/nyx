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
use hifitime::{Epoch, TimeSeries, TimeUnits};
use parquet::arrow::ArrowWriter;
use rand::{Rng, SeedableRng};
use rand_pcg::Pcg64Mcg;
use serde_derive::{Deserialize, Serialize};
use std::error::Error;
use std::fs::File;
use std::ops::{Mul, MulAssign};
use std::path::Path;
use std::sync::Arc;

pub mod gauss_markov;
pub mod white;

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
            bias: Some(GaussMarkov::default_range_km()),
        }
    }

    /// Default stochastic process of the Deep Space Network, using as per DESCANSO Chapter 3, Table 3-3 for the GM process.
    pub fn default_doppler_km_s() -> Self {
        Self {
            white_noise: Some(WhiteNoise {
                sigma: 3e-6, // 3 mm/s
                ..Default::default()
            }),
            bias: Some(GaussMarkov::default_doppler_km_s()),
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
            bias: Some(GaussMarkov::default_range_km()),
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

    /// Simulate the configured stochastic model and store the bias in a parquet file.
    /// Python: call as `simulate(path, runs=25, unit=None)` where the path is the output Parquet file, runs is the number of runs, and unit is the unit of the bias, reflected only in the headers of the parquet file.
    ///
    /// The unit is only used in the headers of the parquet file.
    ///
    /// This will simulate the model with "runs" different seeds, sampling the process 500 times for a duration of 5 times the time constant.
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
            let mut rng = Pcg64Mcg::from_os_rng();

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
}

impl Mul<f64> for StochasticNoise {
    type Output = Self;

    fn mul(mut self, rhs: f64) -> Self::Output {
        if let Some(mut wn) = &mut self.white_noise {
            wn *= rhs;
        }
        if let Some(mut gm) = &mut self.bias {
            gm *= rhs;
        }

        self
    }
}

impl MulAssign<f64> for StochasticNoise {
    fn mul_assign(&mut self, rhs: f64) {
        *self = *self * rhs;
    }
}

pub struct StochasticState {
    pub run: u32,
    pub dt_s: f64,
    pub sample: f64,
    pub variance: f64,
}

#[cfg(test)]
mod ut_stochastics {
    use std::path::PathBuf;

    use super::{white::WhiteNoise, StochasticNoise};

    #[test]
    fn test_simulate_zero() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "data",
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
            "data",
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
            "data",
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
            "data",
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
