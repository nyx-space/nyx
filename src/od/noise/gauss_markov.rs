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
use crate::io::{duration_from_str, duration_to_str};
use crate::linalg::{Matrix2, Vector2};
use crate::NyxError;
use arrow::array::{ArrayRef, Float64Array, UInt32Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use hifitime::{Duration, TimeUnits};
use parquet::arrow::ArrowWriter;
#[cfg(feature = "python")]
use pyo3::prelude::*;
use rand::{Rng, SeedableRng};
use rand_distr::{Distribution, Normal};
use rand_pcg::Pcg64Mcg;
use serde_derive::{Deserialize, Serialize};
use std::collections::HashMap;
use std::fmt;
use std::fs::File;
use std::sync::Arc;

/// A Gauss-Markov process for modeling biases.
///
/// If the natural frequency of the process is zero, then the process is a first-order Gauss Markov process as described in section 5.2.4 of the NASA Best Practices for Navigation Filters (D'Souza et al.).
/// Otherwise, it is a bias and drift couple first- and second-order Gauss Markov process as described in section 5.3.3 of the NASA Best Practices for Navigation Filters (D'Souza et al.).
/// It is up to the caller to ensure that the units at initialization match the units used where the model is applied.
///
/// The model is as follows, where b and d are the bias and drift, and σ_b and σ_d are the standard deviations of the zero-mean white noise of the bias and drift, respectively.
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
#[cfg_attr(
    feature = "python",
    pyo3(
        text_signature = "(tau, bias_sigma, omega_n_hz=None, damping_ratio=None, drift_sigma=None, bias=0.0, drift=0.0)"
    )
)]
pub struct GaussMarkov {
    /// Time constant of the Gauss-Markov process.
    #[serde(
        serialize_with = "duration_to_str",
        deserialize_with = "duration_from_str"
    )]
    pub tau: Duration,
    /// Standard deviation of the zero-mean white noise of the bias.
    pub bias_sigma: f64,
    /// Natural frequency of the drift Gauss-Markov process.
    pub omega_n_hz: Option<f64>,
    /// Damping ratio of the drift Gauss-Markov process.
    pub damping_ratio: Option<f64>,
    /// Standard deviation of the zero-mean white noise of the drift.
    pub drift_sigma: Option<f64>,
    /// Latest bias estimate
    #[serde(skip)]
    pub bias: f64,
    /// Latest drift estimate
    #[serde(skip)]
    pub drift: f64,
}

impl fmt::Display for GaussMarkov {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> fmt::Result {
        match self.omega_n_hz {
            None => {
                write!(
                    f,
                    "First order Gauss-Markov process with τ = {} and σ = {}: bias = {}",
                    self.tau, self.bias_sigma, self.bias
                )
            }
            Some(omega_n_hz) => {
                write!(
                    f,
                    "Coupled bias and drift Gauss-Markov process with τ = {}, ω_n = {} Hz, ζ = {}, σ_bias = {} σ_drift = {}: bias = {}, drift = {}",
                    self.tau, omega_n_hz, self.damping_ratio.unwrap(), self.bias_sigma, self.drift_sigma.unwrap(), self.bias, self.drift
                )
            }
        }
    }
}

impl GaussMarkov {
    /// Create a new Gauss-Markov process.
    ///
    /// # Arguments
    ///
    /// * `tau` - Time constant of the Gauss-Markov process.
    /// * `bias_sigma` - Standard deviation of the zero-mean white noise of the bias.
    /// * `omega_n_hz` - Natural frequency of the drift Gauss-Markov process.
    /// * `damping_ratio` - Damping ratio of the drift Gauss-Markov process.
    /// * `drift_sigma` - Standard deviation of the zero-mean white noise of the drift.
    /// * `bias` - Latest bias estimate.
    /// * `drift` - Latest drift estimate.
    ///
    /// # Panics
    ///
    /// If `tau` is zero.
    /// If `bias_sigma` is negative.
    /// If `omega_n_hz` is negative.
    /// If `damping_ratio` is negative.
    /// If `drift_sigma` is negative.
    pub fn new(
        tau: Duration,
        bias_sigma: f64,
        omega_n_hz: Option<f64>,
        damping_ratio: Option<f64>,
        drift_sigma: Option<f64>,
        init_bias: f64,
        init_drift: f64,
    ) -> Self {
        assert!(tau > Duration::ZERO);
        assert!(bias_sigma >= 0.0);
        if let Some(omega_n_hz) = omega_n_hz {
            assert!(omega_n_hz >= 0.0);
        }
        if let Some(damping_ratio) = damping_ratio {
            assert!(damping_ratio >= 0.0);
        }
        if let Some(drift_sigma) = drift_sigma {
            assert!(drift_sigma >= 0.0);
        }
        Self {
            tau,
            bias_sigma,
            omega_n_hz,
            damping_ratio,
            drift_sigma,
            bias: init_bias,
            drift: init_drift,
        }
    }

    /// Create a new first-order Gauss-Markov process.
    ///
    /// # Arguments
    ///
    /// * `tau` - Time constant of the Gauss-Markov process.
    /// * `bias_sigma` - Standard deviation of the zero-mean white noise of the bias.
    /// * `bias` - Latest bias estimate.
    ///
    pub fn first_order(tau: Duration, bias_sigma: f64, bias: f64) -> Self {
        Self::new(tau, bias_sigma, None, None, None, bias, 0.0)
    }

    /// Create a new bias and drift coupled Gauss-Markov process. The bias is a second order order GM process and the couple drift is a first order GM process.
    ///
    /// # Arguments
    ///
    /// * `tau` - Time constant of the Gauss-Markov process.
    /// * `bias_sigma` - Standard deviation of the zero-mean white noise of the bias.
    /// * `omega_n_hz` - Natural frequency of the drift Gauss-Markov process.
    /// * `damping_ratio` - Damping ratio of the drift Gauss-Markov process.
    /// * `drift_sigma` - Standard deviation of the zero-mean white noise of the drift.
    /// * `bias` - Latest bias estimate.
    /// * `drift` - Latest drift estimate.
    ///
    /// # Panics
    ///
    /// If `tau` is zero.
    /// If `bias_sigma` is negative.
    /// If `omega_n_hz` is negative.
    /// If `damping_ratio` is negative.
    /// If `drift_sigma` is negative.
    pub fn coupled_bias_drift(
        tau: Duration,
        bias_sigma: f64,
        omega_n_hz: f64,
        damping_ratio: f64,
        drift_sigma: f64,
        init_bias: f64,
        init_drift: f64,
    ) -> Self {
        Self::new(
            tau,
            bias_sigma,
            Some(omega_n_hz),
            Some(damping_ratio),
            Some(drift_sigma),
            init_bias,
            init_drift,
        )
    }

    pub fn next_bias<R: Rng>(&mut self, rng: &mut R) -> f64 {
        let omega_fact = match self.omega_n_hz {
            None => 0.0,
            Some(omega_n_hz) => -omega_n_hz.powi(2),
        };

        let zeta_fact = match self.damping_ratio {
            None => 0.0,
            Some(damping_ratio) => -2.0 * damping_ratio * self.omega_n_hz.unwrap(),
        };

        let a_mat = Matrix2::new(-1.0 / self.tau.to_seconds(), 1.0, omega_fact, zeta_fact);

        let b_mat = Matrix2::identity();

        let wn_bias = Normal::new(0.0, self.bias_sigma).unwrap();
        let wn_drift = Normal::new(0.0, self.drift_sigma.unwrap_or(0.0)).unwrap();

        let w_vec = Vector2::new(wn_bias.sample(rng), wn_drift.sample(rng));

        let rates = a_mat * Vector2::new(self.bias, self.drift) + b_mat * w_vec;

        // Update the bias and drift
        self.bias += rates[0];
        self.drift += rates[1];

        // Return the new bias
        self.bias
    }

    /// Zero noise Gauss-Markov process.
    pub const ZERO: Self = Self {
        tau: Duration::MAX,
        bias_sigma: 0.0,
        omega_n_hz: None,
        damping_ratio: None,
        drift_sigma: None,
        bias: 0.0,
        drift: 0.0,
    };

    /// Typical noise on the ranging data from a non-high-precision ground station.
    pub fn default_range_km() -> Self {
        Self {
            tau: 5.minutes(),
            bias_sigma: 50.0e-3, // 50 m
            omega_n_hz: None,
            damping_ratio: None,
            drift_sigma: Some(5e-3), // 5 m
            bias: 0.0,
            drift: 0.0,
        }
    }

    /// Typical noise on the Doppler data from a non-high-precision ground station.
    pub fn default_doppler_km_s() -> Self {
        Self {
            tau: 20.minutes(),
            bias_sigma: 50.0e-5, // 50 cm/s
            omega_n_hz: None,
            damping_ratio: None,
            drift_sigma: Some(7.5e-5), // 7.5 cm/s
            bias: 0.0,
            drift: 0.0,
        }
    }

    /// Example noise on the ranging data from a high-precision ground station.
    pub fn high_precision_range_km() -> Self {
        Self {
            tau: 12.hours(),
            bias_sigma: 5.0e-4, // 0.5 m
            omega_n_hz: None,
            damping_ratio: None,
            drift_sigma: Some(1e-4), // 0.1 m
            bias: -4e-3,
            drift: 0.0,
        }
    }

    /// Example noise on the Doppler data from a high-precision ground station.
    pub fn high_precision_doppler_km_s() -> Self {
        Self {
            tau: 12.hours(),
            bias_sigma: 5.0e-5, // 5 cm/s
            omega_n_hz: None,
            damping_ratio: None,
            drift_sigma: Some(1.5e-6), // 0.15 cm/s
            bias: -4.0e-5,
            drift: 0.0,
        }
    }
}

#[cfg_attr(feature = "python", pymethods)]
impl GaussMarkov {
    /// Simulate a Gauss Markov model and store the bias and drift in a parquet file.
    ///
    /// The unit is only used in the headers of the parquet file.
    ///
    /// This will simulate the model with 1000 different seeds, sampling the process 1140 times.
    /// This corresponds to one sample per minute for 24 hours.
    /// TODO: Add text_signature
    pub fn simulate(&self, path: String, unit: Option<String>) -> Result<(), NyxError> {
        struct GMState {
            run: u32,
            sample: u32,
            bias: f64,
            drift: f64,
        }

        let mut samples = vec![];

        for run in 0..1000 {
            let mut rng = Pcg64Mcg::from_entropy();

            let mut gm = self.clone();
            for sample in 0..1440 {
                gm.next_bias(&mut rng);
                samples.push(GMState {
                    run,
                    sample,
                    bias: gm.bias,
                    drift: gm.drift,
                });
            }
        }

        let unit = match unit {
            Some(unit) => format!("{unit}"),
            None => "(unitless)".to_string(),
        };

        // Build the parquet file
        let hdrs = vec![
            Field::new("Run", DataType::UInt32, false),
            Field::new("Sample", DataType::UInt32, false),
            Field::new(format!("Bias {unit}"), DataType::Float64, false),
            Field::new(format!("Drift {unit}"), DataType::Float64, false),
        ];

        let schema = Arc::new(Schema::new(hdrs));
        let mut record = Vec::new();

        record.push(Arc::new(UInt32Array::from(
            samples.iter().map(|s| s.run).collect::<Vec<u32>>(),
        )) as ArrayRef);
        record.push(Arc::new(UInt32Array::from(
            samples.iter().map(|s| s.sample).collect::<Vec<u32>>(),
        )) as ArrayRef);
        record.push(Arc::new(Float64Array::from(
            samples.iter().map(|s| s.drift).collect::<Vec<f64>>(),
        )) as ArrayRef);
        record.push(Arc::new(Float64Array::from(
            samples.iter().map(|s| s.bias).collect::<Vec<f64>>(),
        )) as ArrayRef);

        // Serialize all of the devices and add that to the parquet file too.
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
    fn py_new(
        tau: Duration,
        bias_sigma: f64,
        omega_n_hz: Option<f64>,
        damping_ratio: Option<f64>,
        drift_sigma: Option<f64>,
        bias: Option<f64>,
        drift: Option<f64>,
    ) -> Self {
        Self {
            tau,
            bias_sigma,
            omega_n_hz,
            damping_ratio,
            drift_sigma,
            bias: bias.unwrap_or(0.0),
            drift: drift.unwrap_or(0.0),
        }
    }
}

#[test]
fn fogm_test() {
    use hifitime::TimeUnits;
    use rstats::{triangmat::Vecops, Stats};

    let mut gm = GaussMarkov::first_order(24.hours(), 0.1, 0.0);

    let mut biases = Vec::with_capacity(1000);

    let mut rng = Pcg64Mcg::new(0);
    for _ in 0..1000 {
        biases.push(gm.next_bias(&mut rng));
    }

    // Result was inspected visually, not sure how to correctly test this.
    let min_max = biases.minmax();

    assert_eq!(biases.amean().unwrap(), -1.64625561108865);
    assert_eq!(min_max.min, -5.978770428734303);
    assert_eq!(min_max.max, 1.137123841255601);
}

#[test]
fn coupled_gm_test() {
    use hifitime::TimeUnits;
    use rand_pcg::Pcg64Mcg;
    use rstats::{triangmat::Vecops, Stats};

    let mut gm = GaussMarkov::coupled_bias_drift(24.hours(), 0.1, 0.1, 0.1, 0.1, 0.0, 0.0);

    let mut biases = Vec::with_capacity(1000);

    let mut rng = Pcg64Mcg::new(0);
    for _ in 0..1000 {
        biases.push(gm.next_bias(&mut rng));
    }

    // Result was inspected visually, not sure how to correctly test this.
    let min_max = biases.minmax();

    assert_eq!(biases.amean().unwrap(), -0.19487332357366507);
    assert_eq!(min_max.min, -17.878511214309935);
    assert_eq!(min_max.max, 18.13775350010581);
}

#[test]
fn zero_noise_test() {
    use rstats::{triangmat::Vecops, Stats};

    let mut gm = GaussMarkov::first_order(Duration::MAX, 0.0, 0.0);

    let mut biases = Vec::with_capacity(1000);

    let mut rng = Pcg64Mcg::new(0);
    for _ in 0..1000 {
        biases.push(gm.next_bias(&mut rng));
    }

    let min_max = biases.minmax();

    assert_eq!(biases.amean().unwrap(), 0.0);
    assert_eq!(min_max.min, 0.0);
    assert_eq!(min_max.max, 0.0);
}

#[test]
fn serde_test() {
    use serde_yaml;
    use std::str::FromStr;

    // Note that we set the initial bias to zero because it is not serialized.
    let gm = GaussMarkov::first_order(Duration::MAX, 0.1, 0.0);
    let serialized = serde_yaml::to_string(&gm).unwrap();
    println!("{serialized}");
    let gm_deser: GaussMarkov = serde_yaml::from_str(&serialized).unwrap();
    assert_eq!(gm_deser, gm);

    // Deserialize a coupled Gauss-Markov model
    let s = r#"
    tau: 1 hour 15 ns
    bias_sigma: 0.1
    omega_n_hz: 1.65
    damping_ratio: 0.1
    drift_sigma: 0.2
    "#;

    let gm_deser: GaussMarkov = serde_yaml::from_str(&s).unwrap();
    let gm = GaussMarkov::coupled_bias_drift(
        Duration::from_str("1 hour 15 ns").unwrap(),
        0.1,
        1.65,
        0.1,
        0.2,
        0.0,
        0.0,
    );

    assert_eq!(gm_deser, gm);
}
