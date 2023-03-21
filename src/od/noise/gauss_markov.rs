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

use crate::linalg::{Matrix2, Vector2};
use hifitime::{Duration, Epoch};
use rand::Rng;
use rand_distr::{Distribution, Normal};
use std::fmt;

/// A Gauss-Markov process for modeling biases.
///
/// If the natural frequency of the process is zero, then the process is a first-order Gauss Markov process as described in section 5.2.4 of the NASA Best Practices for Navigation Filters (D'Souza et al.).
/// Otherwise, it is a bias and drift couple first- and second-order Gauss Markov process as described in section 5.3.3 of the NASA Best Practices for Navigation Filters (D'Souza et al.).
/// It is up to the caller to ensure that the units at initialization match the units used where the model is applied.
#[derive(Copy, Clone, Debug)]
pub struct GaussMarkov {
    /// Time constant of the Gauss-Markov process.
    pub tau: Duration,
    /// Standard deviation of the zero-mean white noise of the bias.
    pub bias_sigma: f64,
    /// Natural frequency of the drift Gauss-Markov process.
    pub omega_n_hz: Option<f64>,
    /// Damping ratio of the drift Gauss-Markov process.
    pub damping_ratio: Option<f64>,
    /// Standard deviation of the zero-mean white noise of the drift.
    pub drift_sigma: Option<f64>,
    /// Latest update epoch
    pub epoch: Epoch,
    /// Latest bias estimate
    pub bias: f64,
    /// Latest drift estimate
    pub drift: f64,
}

impl fmt::Display for GaussMarkov {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> fmt::Result {
        match self.omega_n_hz {
            None => {
                write!(
                    f,
                    "First order Gauss-Markov process with τ = {} and σ = {}: bias = {} @ {}",
                    self.tau, self.bias_sigma, self.bias, self.epoch
                )
            }
            Some(omega_n_hz) => {
                write!(
                    f,
                    "Coupled bias and drift Gauss-Markov process with τ = {}, ω_n = {} Hz, ζ = {}, σ_bias = {} σ_drift = {}: bias = {}, drift = {} @ {}",
                    self.tau, omega_n_hz, self.damping_ratio.unwrap(), self.bias_sigma, self.drift_sigma.unwrap(), self.bias, self.drift, self.epoch
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
    /// * `epoch` - Latest update epoch.
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
        epoch: Epoch,
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
            epoch,
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
    /// * `epoch` - Latest update epoch.
    /// * `bias` - Latest bias estimate.
    ///
    pub fn first_order(tau: Duration, bias_sigma: f64, epoch: Epoch, bias: f64) -> Self {
        Self::new(tau, bias_sigma, None, None, None, epoch, bias, 0.0)
    }

    pub fn next_bias<R: Rng>(&mut self, epoch: Epoch, rng: &mut R) -> f64 {
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
        self.epoch = epoch;

        // Return the new bias
        self.bias
    }
}

#[test]
fn fogm_test() {
    use hifitime::TimeUnits;
    use rand_pcg::Pcg64Mcg;
    use rstats::{triangmat::Vecops, Stats};

    let mut epoch = Epoch::now().unwrap();

    let mut gm = GaussMarkov::first_order(24.hours(), 0.1, epoch, 0.0);

    let mut biases = Vec::with_capacity(1000);

    let mut rng = Pcg64Mcg::new(0);
    for _ in 0..1000 {
        epoch = epoch + 1.hours();
        biases.push(gm.next_bias(epoch, &mut rng));
    }

    // Result was inspected visually, not sure how to correctly test this.
    let min_max = biases.minmax();

    assert_eq!(biases.amean().unwrap(), -1.64625561108865);
    assert_eq!(min_max.min, -5.978770428734303);
    assert_eq!(min_max.max, 1.137123841255601);
}
