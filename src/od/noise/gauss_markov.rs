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
use crate::io::{duration_from_str, duration_to_str, ConfigError, ConfigRepr, Configurable};
use crate::md::ui::Cosm;
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
/// The process is defined by the following stochastic differential equation:
///
/// \dot{b(t)} = -1/τ * b(t) + w(t)
///
/// Programmatically, it's calculated as follows:
///
/// b(t + Δt) = b(t) * exp(-Δt / τ) + q * (1 - exp(-Δt / τ)) * w(t)
///
/// Where w(t) ~ 𝓝(0, σ_{ss}) is a zero-mean white noise process with standard deviation σ_ss, the steady state sigma.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(text_signature = "(tau, sigma, state_state)"))]
pub struct GaussMarkov {
    /// The time constant, tau gives the correlation time, or the time over which the intensity of the time correlation will fade to 1/e of its prior value. (This is sometimes incorrectly referred to as the "half-life" of the process.)
    #[serde(
        serialize_with = "duration_to_str",
        deserialize_with = "duration_from_str"
    )]
    pub tau: Duration,
    /// Standard deviation (or covariance) of the zero-mean white noise of the initial bias, noted `σ_b`.
    pub bias_sigma: f64,
    /// The steady-state sigma is the zero-mean white noise as t → ∞, noted `σ_q`, and sometimes called the "constant".
    pub steady_state_sigma: f64,
    /// Latest bias, unset prior to the sample call and one will be generated from a zero mean normal distribution with standard deviation `bias_sigma`.
    #[serde(skip)]
    pub bias: Option<f64>,
    /// Epoch of the latest sample of the process.
    #[serde(skip)]
    pub epoch: Option<Epoch>,
}

impl fmt::Display for GaussMarkov {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "First order Gauss-Markov process with τ = {}, σ_b = {}, σ_q = {}",
            self.tau, self.bias_sigma, self.steady_state_sigma
        )
    }
}

impl GaussMarkov {
    /// Create a new first order Gauss-Markov process.
    /// # Arguments
    /// * `tau` - The time constant, tau gives the correlation time, or the time over which the intensity of the time correlation will fade to 1/e of its prior value.
    /// * `bias_sigma` - Standard deviation (or covariance) of the zero-mean white noise of the initial bias.
    /// * `steady_state_sigma` - The steady-state sigma is the zero-mean white noise as t → ∞, noted `q`, and sometimes called the "constant".
    pub fn new(
        tau: Duration,
        bias_sigma: f64,
        steady_state_sigma: f64,
    ) -> Result<Self, ConfigError> {
        if tau <= Duration::ZERO {
            return Err(ConfigError::InvalidConfig(format!(
                "tau must be positive but got {tau}"
            )));
        }

        Ok(Self {
            tau,
            bias_sigma,
            steady_state_sigma,
            bias: None,
            epoch: None,
        })
    }

    /// Return the next bias sample.
    pub fn next_bias<R: Rng>(&mut self, epoch: Epoch, rng: &mut R) -> f64 {
        // Compute the delta time in seconds between the previous epoch and the sample epoch.
        let dt_s = (match self.epoch {
            None => Duration::ZERO,
            Some(prev_epoch) => epoch - prev_epoch,
        })
        .to_seconds();
        self.epoch = Some(epoch);

        // If there is no bias, generate one using the standard deviation of the bias
        if self.bias.is_none() {
            self.bias = Some(rng.sample(Normal::new(0.0, self.bias_sigma).unwrap()));
        }

        // For some reason, there is a different expression for the decay and anti decay.
        let decay = (-dt_s / self.tau.to_seconds()).exp();
        let anti_decay = 1.0 - decay;

        // The steady state contribution. This is the bias that the process will converge to as t approaches infinity.
        let ss_contrib = rng.sample(Normal::new(0.0, self.steady_state_sigma).unwrap());

        let bias = self.bias.unwrap() * decay + ss_contrib * anti_decay;

        self.bias = Some(bias);

        // Return the new bias
        bias
    }

    /// Zero noise Gauss-Markov process.
    pub const ZERO: Self = Self {
        tau: Duration::MAX,
        bias_sigma: 0.0,
        steady_state_sigma: 0.0,
        bias: None,
        epoch: None,
    };

    /// Typical noise on the ranging data from a non-high-precision ground station.
    pub fn default_range_km() -> Self {
        Self {
            tau: 5.minutes(),
            bias_sigma: 50.0e-3,      // 50 m
            steady_state_sigma: 5e-3, // 5 m
            bias: None,
            epoch: None,
        }
    }

    /// Typical noise on the Doppler data from a non-high-precision ground station.
    pub fn default_doppler_km_s() -> Self {
        Self {
            tau: 20.minutes(),
            bias_sigma: 500.0e-6,      // 50 cm/s
            steady_state_sigma: 75e-6, // 7.5 cm/s
            bias: None,
            epoch: None,
        }
    }

    /// Example noise on the ranging data from a high-precision ground station, e.g. NASA Deep Space Network (DSN).
    pub fn high_precision_range_km() -> Self {
        Self {
            tau: 12.hours(),
            bias_sigma: 5.0e-3,         // 5 m
            steady_state_sigma: 0.1e-3, // 0.1 m
            bias: None,
            epoch: None,
        }
    }

    /// Example noise on the Doppler data from a high-precision ground station, e.g. NASA Deep Space Network (DSN).
    pub fn high_precision_doppler_km_s() -> Self {
        Self {
            tau: 12.hours(),
            bias_sigma: 50.0e-6,        // 5 cm/s
            steady_state_sigma: 1.5e-6, // 0.15 cm/s
            bias: None,
            epoch: None,
        }
    }

    /// Initializes a new Gauss Markov process for the provided kind of model.
    ///
    /// Available models are: `Range`, `Doppler`, `RangeHP`, `Doppler HP` (HP stands for high precision).
    pub fn from_default(kind: String) -> Result<Self, NyxError> {
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

#[cfg_attr(feature = "python", pymethods)]
impl GaussMarkov {
    /// Simulate a Gauss Markov model and store the bias in a parquet file.
    /// Python: call as `simulate(path, runs=25, unit=None)` where the path is the output Parquet file, runs is the number of runs, and unit is the unit of the bias, reflected only in the headers of the parquet file.
    ///
    /// The unit is only used in the headers of the parquet file.
    ///
    /// This will simulate the model with "runs" different seeds, sampling the process 500 times for a duration of 5 times the time constant.
    pub fn simulate(
        &self,
        path: String,
        runs: Option<u32>,
        unit: Option<String>,
    ) -> Result<(), NyxError> {
        struct GMState {
            run: u32,
            dt_s: f64,
            bias: f64,
        }

        let num_runs = runs.unwrap_or(25);

        let mut samples = Vec::with_capacity(num_runs as usize);
        let start = Epoch::now().unwrap();
        // Run until the steady state is clearly visible.
        let end = start + self.tau * 5;
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
    fn default(kind: String) -> Result<Self, NyxError> {
        Self::from_default(kind)
    }

    #[cfg(feature = "python")]
    #[staticmethod]
    fn load(path: &str) -> Result<Self, ConfigError> {
        <Self as ConfigRepr>::load(path)
    }

    #[cfg(feature = "python")]
    #[staticmethod]
    fn load_many(path: &str) -> Result<Vec<Self>, ConfigError> {
        <Self as ConfigRepr>::load_many(path)
    }

    #[cfg(feature = "python")]
    #[staticmethod]
    fn load_named(path: &str) -> Result<HashMap<String, Self>, ConfigError> {
        <Self as ConfigRepr>::load_named(path)
    }
}

impl ConfigRepr for GaussMarkov {}

impl Configurable for GaussMarkov {
    type IntermediateRepr = Self;

    fn from_config(
        cfg: Self::IntermediateRepr,
        _cosm: Arc<Cosm>,
    ) -> Result<Self, crate::io::ConfigError>
    where
        Self: Sized,
    {
        Ok(cfg)
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, crate::io::ConfigError> {
        let serded: Self::IntermediateRepr = (*self).into();
        Ok(serded)
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

    assert_eq!(biases.amean().unwrap(), 0.09422989888670787);
    assert_eq!(min_max.min, 0.09368618104727605);
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
    use std::env;
    use std::path::PathBuf;

    // Note that we set the initial bias to zero because it is not serialized.
    let gm = GaussMarkov::new(Duration::MAX, 0.1, 0.0).unwrap();
    let serialized = serde_yaml::to_string(&gm).unwrap();
    println!("{serialized}");
    let gm_deser: GaussMarkov = serde_yaml::from_str(&serialized).unwrap();
    assert_eq!(gm_deser, gm);

    let test_data: PathBuf = [
        env::var("CARGO_MANIFEST_DIR").unwrap(),
        "data".to_string(),
        "tests".to_string(),
        "config".to_string(),
        "high-prec-network.yaml".to_string(),
    ]
    .iter()
    .collect();

    let models = GaussMarkov::load_named(test_data).unwrap();
    assert_eq!(models.len(), 2);
    assert_eq!(
        models["range_noise_model"].tau,
        12.hours() + 159.milliseconds()
    );
    assert_eq!(models["range_noise_model"].bias_sigma, 5.0e-3);
    assert_eq!(models["range_noise_model"].steady_state_sigma, 0.1e-3);

    assert_eq!(models["doppler_noise_model"].tau, 11.hours() + 59.minutes());
    assert_eq!(models["doppler_noise_model"].bias_sigma, 50.0e-6);
    assert_eq!(models["doppler_noise_model"].steady_state_sigma, 1.5e-6);
}
