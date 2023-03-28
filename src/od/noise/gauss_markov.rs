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

use crate::io::watermark::pq_writer;
use crate::io::{duration_from_str, duration_to_str, ConfigError};
use crate::NyxError;
use arrow::array::{ArrayRef, Float64Array, UInt32Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use hifitime::{Duration, Epoch, TimeSeries, TimeUnits};
use parquet::arrow::ArrowWriter;
#[cfg(feature = "python")]
use pyo3::prelude::*;
use rand::{Rng, SeedableRng};
use rand_distr::Normal;
use rand_pcg::Pcg64Mcg;
use serde_derive::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::sync::Arc;

/// A first order Gauss-Markov process for modeling biases as described in section 5.2.4 of the NASA Best Practices for Navigation Filters (D'Souza et al.).
///
/// The model is as follows, where b(t) is the bias at epoch t, and σ_b and σ_d are the standard deviations of the zero-mean white noise of the bias and drift, respectively.
/// The time constant τ is only used as a time constant on the bias, and is not dependent on the absolute epoch, only on the original drift.
/// ω_n^2 is the natural frequency of the drift Gauss-Markov process, and ζ is the damping ratio of the drift Gauss-Markov process.
///
/// \dot{b(t)} = -1/τ b(t) + d(t) + σ_b(t)
///
/// \dot{d(t)} = -ω_n^2 b(t) - 2ζω_n d(t) + σ_d(t)
///
/// # Implementation notes
///
/// 1. There is no epoch in this Gauss Markov process. The time constant τ is only used as a time constant on the bias.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(text_signature = "(tau, sigma, state_state)"))]
pub struct GaussMarkov {
    /// The time constant, tau gives the correlation time, or the time over which the intensity of the time correlation will fade to 1/ e of its prior value 2.
    #[serde(
        serialize_with = "duration_to_str",
        deserialize_with = "duration_from_str"
    )]
    pub tau: Duration,
    /// Standard deviation (or covariance) of the zero-mean white noise of the bias, and updated at each sample.
    pub sigma: f64,
    /// The steady-state bias on which this process converges as t approaches infinity, noted `q`, and sometimes called the "constant".
    pub steady_state: f64,
    /// Latest bias, if unset prior to the first call, then one will be generated from a zero mean normal distribution with standard deviation `bias_sigma`.
    #[serde(skip)]
    pub bias: Option<f64>,
    /// Initial epoch, unset before first sample.
    #[serde(skip)]
    pub epoch: Option<Epoch>,
}

impl fmt::Display for GaussMarkov {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "First order Gauss-Markov process with τ = {}, σ = {}, q = {}: bias = {:?}",
            self.tau, self.sigma, self.steady_state, self.bias
        )
    }
}

impl GaussMarkov {
    /// Create a new first order Gauss-Markov process.
    /// # Arguments
    /// * `tau` - The time constant, tau gives the correlation time, or the time over which the intensity of the time correlation will fade to 1/ e of its prior value 2.
    /// * `sigma` - Standard deviation (or covariance) of the zero-mean white noise of the bias, and updated at each sample.
    /// * `steady_state` - The steady-state bias on which this process converges as t approaches infinity, noted `q`, and sometimes called the "constant".
    pub fn new(tau: Duration, sigma: f64, steady_state: f64) -> Result<Self, ConfigError> {
        if tau <= Duration::ZERO {
            return Err(ConfigError::InvalidConfig(format!(
                "tau must be positive but got {tau}"
            )));
        }
        if sigma <= 0.0 {
            return Err(ConfigError::InvalidConfig(format!(
                "sigma must be positive but got {sigma}"
            )));
        }

        Ok(Self {
            tau,
            sigma,
            steady_state,
            bias: None,
            epoch: None,
        })
    }

    pub fn next_bias<R: Rng>(&mut self, epoch: Epoch, rng: &mut R) -> f64 {
        // Compute the delta time in seconds between the very first epoch and the sample epoch.
        let dt_s = (match self.epoch {
            None => {
                self.epoch = Some(epoch);
                Duration::ZERO
            }
            Some(prev_epoch) => epoch - prev_epoch,
        })
        .to_seconds();

        // If there is no bias, generate one using the standard deviation of the bias
        if self.bias.is_none() {
            self.bias = Some(rng.sample(Normal::new(0.0, self.sigma).unwrap()));
        }

        // For some reason, there is a different expression for the decay and anti decay.
        let decay = (-2.0 * dt_s / self.tau.to_seconds()).exp();
        let anti_decay = 1.0 - decay;

        dbg!(decay, anti_decay);

        // let ss_stddev = (self.steady_state * self.tau.to_seconds() / 2.0) * anti_decay;
        // The steady state contribution. This is the bias that the process will converge to as t approaches infinity.
        let ss_contrib = rng.sample(Normal::new(0.0, self.sigma).unwrap());
        let ss = self.steady_state * anti_decay * ss_contrib;

        let bias = self.bias.unwrap() * decay + ss;

        self.bias = Some(bias);

        // Return the new bias
        bias
    }

    /// Zero noise Gauss-Markov process.
    pub const ZERO: Self = Self {
        tau: Duration::MAX,
        sigma: 0.0,
        steady_state: 0.0,
        bias: None,
        epoch: None,
    };

    /// Typical noise on the ranging data from a non-high-precision ground station.
    pub fn default_range_km() -> Self {
        Self {
            tau: 5.minutes(),
            sigma: 50.0e-3,     // 50 m
            steady_state: 5e-3, // 5 m
            bias: None,
            epoch: None,
        }
    }

    /// Typical noise on the Doppler data from a non-high-precision ground station.
    pub fn default_doppler_km_s() -> Self {
        Self {
            tau: 20.minutes(),
            sigma: 50.0e-6,      // 50 cm/s
            steady_state: 75e-6, // 7.5 cm/s
            bias: None,
            epoch: None,
        }
    }

    /// Example noise on the ranging data from a high-precision ground station.
    pub fn high_precision_range_km() -> Self {
        Self {
            tau: 12.hours(),
            sigma: 50.0e-6,     // 5 cm
            steady_state: 1e-4, // 0.1 m
            bias: Some(-4e-3),  // -0.4 m
            epoch: None,
        }
    }

    /// Example noise on the Doppler data from a high-precision ground station.
    pub fn high_precision_doppler_km_s() -> Self {
        Self {
            tau: 12.hours(),
            sigma: 50.0e-6,       // 5 cm/s
            steady_state: 1.5e-6, // 0.15 cm/s
            bias: Some(-40e-6),   // -4 cm/s
            epoch: None,
        }
    }
}

#[cfg_attr(feature = "python", pymethods)]
impl GaussMarkov {
    /// Simulate a Gauss Markov model and store the bias and drift in a parquet file.
    ///
    /// The unit is only used in the headers of the parquet file.
    ///
    /// This will simulate the model with 100 different seeds, sampling the process 1140 times.
    /// This corresponds to one sample per minute for 24 hours.
    /// TODO: Add text_signature
    pub fn simulate(&self, path: String, unit: Option<String>) -> Result<(), NyxError> {
        struct GMState {
            run: u32,
            dt_s: f64,
            bias: f64,
        }

        let num_runs: u32 = 2;

        let mut samples = Vec::with_capacity(num_runs as usize);
        let start = Epoch::now().unwrap();
        // Run until the steady state is clearly visible.
        let end = start + self.tau * 3;
        let step = self.tau / 100;

        for run in 0..num_runs {
            let mut rng = Pcg64Mcg::from_entropy();

            let mut gm = self.clone();
            for epoch in TimeSeries::inclusive(start, end, step) {
                let bias = gm.next_bias(epoch, &mut rng);
                samples.push(GMState {
                    run,
                    dt_s: (epoch - start).to_seconds(),
                    bias,
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
        ];

        let schema = Arc::new(Schema::new(hdrs));
        let mut record = Vec::new();

        record.push(Arc::new(UInt32Array::from(
            samples.iter().map(|s| s.run).collect::<Vec<u32>>(),
        )) as ArrayRef);
        record.push(Arc::new(Float64Array::from(
            samples.iter().map(|s| s.dt_s).collect::<Vec<f64>>(),
        )) as ArrayRef);
        record.push(Arc::new(Float64Array::from(
            samples.iter().map(|s| s.bias).collect::<Vec<f64>>(),
        )) as ArrayRef);

        let mut metadata = HashMap::new();
        metadata.insert("Purpose".to_string(), "Gauss Markov simulation".to_string());

        let props = pq_writer(Some(metadata));

        let file = File::create(&path).map_err(|e| NyxError::CustomError(e.to_string()))?;
        let mut writer = ArrowWriter::try_new(file, schema.clone(), props).unwrap();

        let batch = RecordBatch::try_new(schema, record)
            .map_err(|e| NyxError::CustomError(e.to_string()))?;
        writer
            .write(&batch)
            .map_err(|e| NyxError::CustomError(e.to_string()))?;
        writer
            .close()
            .map_err(|e| NyxError::CustomError(e.to_string()))?;

        Ok(())
    }

    #[cfg(feature = "python")]
    pub fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    #[cfg(feature = "python")]
    pub fn __str__(&self) -> String {
        format!("{self}")
    }

    #[cfg(feature = "python")]
    #[new]
    fn py_new(tau: Duration, sigma: f64, steady_state: f64) -> Result<Self, ConfigError> {
        Self::new(tau, sigma, steady_state)
    }

    #[cfg(feature = "python")]
    #[getter]
    fn get_tau(&self) -> Duration {
        self.tau
    }

    #[cfg(feature = "python")]
    #[getter]
    fn get_bias(&self) -> Option<f64> {
        self.bias
    }

    /// Initializes a new Gauss Markov process for the provided kind of model.
    ///
    /// Available models are: `Range`, `Doppler`, `RangeHP`, `Doppler HP` (HP stands for high precision).
    #[cfg(feature = "python")]
    #[staticmethod]
    fn from_default(kind: String) -> Result<Self, NyxError> {
        match kind.as_str() {
            "Range" => Ok(Self::default_range_km()),
            "Doppler" => Ok(Self::default_doppler_km_s()),
            "RangeHP" => Ok(Self::high_precision_range_km()),
            "DopplerHP" => Ok(Self::high_precision_doppler_km_s()),
            _ => Err(NyxError::CustomError(format!(
                "No default Gauss Markov model for `{kind}`"
            ))),
        }
    }
}

#[test]
fn fogm_test() {
    use hifitime::TimeUnits;
    use rstats::{triangmat::Vecops, Stats};

    let mut gm = GaussMarkov::new(24.hours(), 0.1, 0.0).unwrap();

    let mut biases = Vec::with_capacity(1000);
    let epoch = Epoch::now().unwrap();

    let mut rng = Pcg64Mcg::new(0);
    for seconds in 0..1000 {
        biases.push(gm.next_bias(epoch + seconds.seconds(), &mut rng));
    }

    // Result was inspected visually with the test_gauss_markov.py Python script
    // I'm not sure how to correctly test this and open to ideas.
    let min_max = biases.minmax();

    assert_eq!(biases.amean().unwrap(), 0.024688771669655347);
    assert_eq!(min_max.min, 9.019181340919657e-7);
    assert_eq!(min_max.max, 0.0947757142410022);
}

#[test]
fn zero_noise_test() {
    use rstats::{triangmat::Vecops, Stats};

    let mut gm = GaussMarkov::ZERO;

    let mut biases = Vec::with_capacity(1000);
    let epoch = Epoch::now().unwrap();

    let mut rng = Pcg64Mcg::new(0);
    for seconds in 0..1000 {
        biases.push(gm.next_bias(epoch + seconds.seconds(), &mut rng));
    }

    let min_max = biases.minmax();

    assert_eq!(biases.amean().unwrap(), 0.0);
    assert_eq!(min_max.min, 0.0);
    assert_eq!(min_max.max, 0.0);
}

#[test]
fn serde_test() {
    use serde_yaml;

    // Note that we set the initial bias to zero because it is not serialized.
    let gm = GaussMarkov::new(Duration::MAX, 0.1, 0.0).unwrap();
    let serialized = serde_yaml::to_string(&gm).unwrap();
    println!("{serialized}");
    let gm_deser: GaussMarkov = serde_yaml::from_str(&serialized).unwrap();
    assert_eq!(gm_deser, gm);
}
