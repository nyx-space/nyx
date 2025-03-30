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

use crate::io::{ConfigError, ConfigRepr};
use hifitime::{Duration, Epoch, TimeUnits};

use rand::Rng;
use rand_distr::Normal;
use serde_derive::{Deserialize, Serialize};
use std::fmt;
use std::ops::{Mul, MulAssign};

use super::Stochastics;

/// A first order Gauss-Markov process for modeling biases as described in section 5.2.4 of the NASA Best Practices for Navigation Filters (D'Souza et al.).
///
/// The process is defined by the following stochastic differential equation:
///
/// \dot{b(t)} = -1/τ * b(t) + w(t)
///
/// Programmatically, it's calculated by sampling from b(t) ~ 𝓝(0, p_b(t)), where
///
/// p_b(t) = exp((-2 / τ) * (t - t_0)) * p_b(t_0) + s(t - t_0)
///
/// s(t - t_0) = ((q * τ) / 2) * (1 - exp((-2 / τ) * (t - t_0)))
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
pub struct GaussMarkov {
    /// The time constant, tau gives the correlation time, or the time over which the intensity of the time correlation will fade to 1/e of its prior value. (This is sometimes incorrectly referred to as the "half-life" of the process.)
    pub tau: Duration,
    pub process_noise: f64,
    /// An optional constant offset on top of the noise, defaults to zero.
    pub constant: Option<f64>,
    /// Epoch of the previous realization, used to compute the time delta for the process noise.
    #[serde(skip)]
    pub prev_epoch: Option<Epoch>,
    /// Sample of previous realization
    #[serde(skip)]
    pub init_sample: Option<f64>,
}

impl fmt::Display for GaussMarkov {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "First order Gauss-Markov process with τ = {}, σ = {}",
            self.tau, self.process_noise
        )
    }
}

impl GaussMarkov {
    /// Create a new first order Gauss-Markov process.
    /// # Arguments
    /// * `tau` - The time constant, tau gives the correlation time, or the time over which the intensity of the time correlation will fade to 1/e of its prior value.
    /// * `process_noise` - process noise of the system.
    pub fn new(tau: Duration, process_noise: f64) -> Result<Self, ConfigError> {
        if tau <= Duration::ZERO {
            return Err(ConfigError::InvalidConfig {
                msg: format!("tau must be positive but got {tau}"),
            });
        }

        Ok(Self {
            tau,
            process_noise,
            constant: None,
            init_sample: None,
            prev_epoch: None,
        })
    }

    /// Zero noise Gauss-Markov process.
    pub const ZERO: Self = Self {
        tau: Duration::MAX,
        process_noise: 0.0,
        constant: None,
        init_sample: None,
        prev_epoch: None,
    };

    /// Default Gauss Markov noise of the Deep Space Network, as per DESCANSO Chapter 3, Table 3-3.
    /// Used the range value of 60 cm over a 60 second average.
    pub fn default_range_km() -> Self {
        Self {
            tau: 1.minutes(),
            process_noise: 60.0e-5,
            constant: None,
            init_sample: None,
            prev_epoch: None,
        }
    }

    /// Default Gauss Markov noise of the Deep Space Network, as per DESCANSO Chapter 3, Table 3-3.
    /// Used the Doppler value of 0.03 mm/s over a 60 second average.
    pub fn default_doppler_km_s() -> Self {
        Self {
            tau: 1.minutes(),
            process_noise: 0.03e-6,
            constant: None,
            init_sample: None,
            prev_epoch: None,
        }
    }
}

impl Stochastics for GaussMarkov {
    fn covariance(&self, _epoch: Epoch) -> f64 {
        self.process_noise.powi(2)
    }

    /// Return the next bias sample.
    fn sample<R: Rng>(&mut self, epoch: Epoch, rng: &mut R) -> f64 {
        // Compute the delta time in seconds between the previous epoch and the sample epoch.
        let dt_s = (match self.prev_epoch {
            None => Duration::ZERO,
            Some(prev_epoch) => epoch - prev_epoch,
        })
        .to_seconds();
        self.prev_epoch = Some(epoch);

        // If there is no bias, generate one using the standard deviation of the bias
        if self.init_sample.is_none() {
            self.init_sample = Some(rng.sample(Normal::new(0.0, self.process_noise).unwrap()));
        }

        let decay = (-dt_s / self.tau.to_seconds()).exp();
        let anti_decay = 1.0 - decay;

        // The steady state contribution. This is the bias that the process will converge to as t approaches infinity.
        let steady_noise = 0.5 * self.process_noise * self.tau.to_seconds() * anti_decay;
        let ss_sample = rng.sample(Normal::new(0.0, steady_noise).unwrap());

        self.init_sample.unwrap() * decay + ss_sample + self.constant.unwrap_or(0.0)
    }
}

impl Mul<f64> for GaussMarkov {
    type Output = Self;

    /// Scale the Gauss Markov process by a constant, maintaining the same time constant.
    fn mul(mut self, rhs: f64) -> Self::Output {
        self.process_noise *= rhs;
        self.constant = None;
        self.init_sample = None;
        self.prev_epoch = None;
        self
    }
}

impl MulAssign<f64> for GaussMarkov {
    fn mul_assign(&mut self, rhs: f64) {
        *self = *self * rhs;
    }
}

impl ConfigRepr for GaussMarkov {}

#[cfg(test)]
mod ut_gm {

    use hifitime::{Duration, Epoch, TimeUnits};
    use rand_pcg::Pcg64Mcg;
    use rstats::{triangmat::Vecops, Stats};

    use crate::{
        io::ConfigRepr,
        od::noise::{GaussMarkov, Stochastics},
    };

    #[test]
    fn fogm_test() {
        let mut gm = GaussMarkov::new(24.hours(), 0.1).unwrap();

        let mut biases = Vec::with_capacity(1000);
        let epoch = Epoch::now().unwrap();

        let mut rng = Pcg64Mcg::new(0);
        for seconds in 0..1000 {
            biases.push(gm.sample(epoch + seconds.seconds(), &mut rng));
        }

        // Result was inspected visually with the test_gauss_markov.py Python script
        // I'm not sure how to correctly test this and open to ideas.
        let min_max = biases.minmax();

        assert_eq!(biases.amean().unwrap(), 0.09373233290645445);
        assert_eq!(min_max.max, 0.24067114622652647);
        assert_eq!(min_max.min, -0.045552031890295525);
    }

    #[test]
    fn zero_noise_test() {
        use rstats::{triangmat::Vecops, Stats};

        let mut gm = GaussMarkov::ZERO;

        let mut biases = Vec::with_capacity(1000);
        let epoch = Epoch::now().unwrap();

        let mut rng = Pcg64Mcg::new(0);
        for seconds in 0..1000 {
            biases.push(gm.sample(epoch + seconds.seconds(), &mut rng));
        }

        let min_max = biases.minmax();

        assert_eq!(biases.amean().unwrap(), 0.0);
        assert_eq!(min_max.min, 0.0);
        assert_eq!(min_max.max, 0.0);
    }

    #[test]
    fn serde_test() {
        use serde_yml;
        use std::env;
        use std::path::PathBuf;

        // Note that we set the initial bias to zero because it is not serialized.
        let gm = GaussMarkov::new(Duration::MAX, 0.1).unwrap();
        let serialized = serde_yml::to_string(&gm).unwrap();
        println!("{serialized}");
        let gm_deser: GaussMarkov = serde_yml::from_str(&serialized).unwrap();
        assert_eq!(gm_deser, gm);

        let test_data: PathBuf = [
            env::var("CARGO_MANIFEST_DIR").unwrap(),
            "data".to_string(),
            "03_tests".to_string(),
            "config".to_string(),
            "high-prec-network.yaml".to_string(),
        ]
        .iter()
        .collect();

        let models = <GaussMarkov as ConfigRepr>::load_named(test_data).unwrap();
        assert_eq!(models.len(), 2);
        assert_eq!(
            models["range_noise_model"].tau,
            12.hours() + 159.milliseconds()
        );
        assert_eq!(models["range_noise_model"].process_noise, 5.0e-3);

        assert_eq!(models["doppler_noise_model"].tau, 11.hours() + 59.minutes());
        assert_eq!(models["doppler_noise_model"].process_noise, 50.0e-6);
    }
}
