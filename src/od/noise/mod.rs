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
use rand_distr::Normal;
use rand_pcg::Pcg64Mcg;
use serde_derive::{Deserialize, Serialize};
use std::error::Error;
use std::fs::File;
use std::path::Path;
use std::sync::Arc;

pub mod gauss_markov;
pub mod walk;

pub use gauss_markov::GaussMarkov;
use walk::RandomWalk;

/// Trait for any kind of stochastic modeling, developing primarily for synthetic orbit determination measurements.
pub trait Stochastics {
    /// Return the variance of this stochastic noise model at a given time.
    fn variance(&mut self, epoch: Epoch) -> f64;

    /// A sampler based on the variance.
    fn sample<R: Rng>(&mut self, epoch: Epoch, rng: &mut R) -> f64 {
        rng.sample(Normal::new(0.0, self.variance(epoch)).unwrap())
    }
}

/// Stochastic noise modeling used primarily for synthetic orbit determination measurements.
///
/// This implementation distinguishes between the white noise model and the bias model. It also includes a constant offset.
#[derive(Copy, Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct StochasticNoise {
    pub white_noise: Option<RandomWalk>,
    pub bias: Option<GaussMarkov>,
    pub constant: Option<f64>,
}

impl StochasticNoise {
    /// Zero noise stochastic process.
    pub const ZERO: Self = Self {
        white_noise: None,
        bias: None,
        constant: None,
    };

    /// Default stochastic process of the Deep Space Network, as per DESCANSO Chapter 3, Table 3-3.
    /// Using the instrument bias as the white noise value, zero constant bias.
    pub fn default_range_km() -> Self {
        Self {
            white_noise: Some(RandomWalk {
                process_noise_per_s: 2.0e-3, // 2 m
                ..Default::default()
            }),
            bias: Some(GaussMarkov::default_range_km()),
            constant: None,
        }
    }

    /// Default stochastic process of the Deep Space Network, using as per DESCANSO Chapter 3, Table 3-3 for the GM process.
    pub fn default_doppler_km_s() -> Self {
        Self {
            white_noise: Some(RandomWalk {
                process_noise_per_s: 0.3e-6, // 3 mm/s
                ..Default::default()
            }),
            bias: Some(GaussMarkov::default_doppler_km_s()),
            constant: None,
        }
    }

    /// Sample these stochastics
    pub fn sample<R: Rng>(&mut self, epoch: Epoch, rng: &mut R) -> f64 {
        rng.sample(Normal::new(0.0, self.variance(epoch)).unwrap()) + self.constant.unwrap_or(0.0)
    }

    /// Return the variance of these stochastics at a given time.
    pub fn variance(&mut self, epoch: Epoch) -> f64 {
        let mut variance = 0.0;
        if let Some(wn) = &mut self.white_noise {
            variance += wn.variance(epoch);
        }
        variance
    }

    /// Simulate a Gauss Markov model and store the bias in a parquet file.
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
    ) -> Result<(), Box<dyn Error>> {
        struct StochasticState {
            run: u32,
            dt_s: f64,
            sample: f64,
            variance: f64,
        }

        let num_runs = runs.unwrap_or(25);

        let mut samples = Vec::with_capacity(num_runs as usize);
        let start = Epoch::now().unwrap();
        // Run until the steady state is clearly visible.
        let (step, end) = if let Some(gm) = self.bias {
            (gm.tau / 100, start + gm.tau * 5)
        } else {
            (1.minutes(), start + 5.days())
        };

        for run in 0..num_runs {
            let mut rng = Pcg64Mcg::from_entropy();

            let mut mdl = self;
            for epoch in TimeSeries::inclusive(start, end, step) {
                let variance = mdl.variance(epoch);
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

        Ok(())
    }
}

#[cfg(test)]
mod ut_stochastics {
    use std::path::PathBuf;

    use super::{walk::RandomWalk, StochasticNoise};

    #[test]
    fn test_simulate_zero() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "output_data",
            "stochastics_zero.parquet",
        ]
        .iter()
        .collect();

        let noise = StochasticNoise::default();

        noise.simulate(path, None, None).unwrap();
    }

    #[test]
    fn test_simulate_constant() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "output_data",
            "stochastics_constant.parquet",
        ]
        .iter()
        .collect();

        let noise = StochasticNoise {
            constant: Some(1.0),
            ..Default::default()
        };

        noise.simulate(path, None, None).unwrap();
    }

    #[test]
    fn test_simulate_wn_only() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "output_data",
            "stochastics_wn.parquet",
        ]
        .iter()
        .collect();

        let noise = StochasticNoise {
            white_noise: Some(RandomWalk {
                process_noise_per_s: 2.1,
                ..Default::default()
            }),
            ..Default::default()
        };

        noise
            .simulate(path, None, Some("meter".to_string()))
            .unwrap();
    }

    #[test]
    fn test_simulate_wn_constant() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "output_data",
            "stochastics_wn_and_constant.parquet",
        ]
        .iter()
        .collect();

        let noise = StochasticNoise {
            white_noise: Some(RandomWalk {
                process_noise_per_s: 2.1,
                ..Default::default()
            }),
            constant: Some(0.5),
            ..Default::default()
        };

        noise
            .simulate(path, None, Some("meter".to_string()))
            .unwrap();
    }

    #[test]
    fn test_simulate_dsn_range() {
        let path: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "output_data",
            "stochastics_dsn_range.parquet",
        ]
        .iter()
        .collect();

        let noise = StochasticNoise::default_range_km();

        noise
            .simulate(path, None, Some("kilometer".to_string()))
            .unwrap();
    }
}
