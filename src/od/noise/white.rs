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

use std::ops::{Mul, MulAssign};

use anise::constants::SPEED_OF_LIGHT_KM_S;
use hifitime::{Duration, Epoch};
use rand::Rng;
use rand_distr::Normal;
use serde_derive::{Deserialize, Serialize};

use super::Stochastics;

/// White noise is an uncorrelated random variable.
#[derive(Copy, Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct WhiteNoise {
    /// Mean value of this white noise
    pub mean: f64,
    /// Process noise as a one-sigma of the Normal distribution.
    pub sigma: f64,
}

impl WhiteNoise {
    /// Initializes a new random walk stochastic noise model from the process noise and the integration time.
    /// This will compute the process noise per second automatically.
    pub fn new(process_noise: f64, integration_time: Duration) -> Self {
        Self {
            sigma: process_noise / integration_time.to_seconds(),
            ..Default::default()
        }
    }

    /// Initializes a new random walk stochastic noise model from the provided process noise, assuming that the noise level
    /// is fixed regardless of the integration time.
    pub fn constant_white_noise(process_noise: f64) -> Self {
        Self {
            sigma: process_noise,
            ..Default::default()
        }
    }

    // Initializes a new time-uncorrelated white noise process, using only the Pr/N0 value and the bandwidth.
    // This returns a white noise sigma in kilometers.
    //
    // # Equation
    // σ = c / (2 * B * √(Pr/N0))
    //
    // Where c is the speed of light, B is the bandwidth in Hz, and the Pr/N0 is the signal-to-noise ratio.
    pub fn from_pr_n0(pr_n0: f64, bandwidth_hz: f64) -> Self {
        Self {
            sigma: SPEED_OF_LIGHT_KM_S / (2.0 * bandwidth_hz * (pr_n0).sqrt()),
            mean: 0.0,
        }
    }
}

impl Stochastics for WhiteNoise {
    fn covariance(&self, _epoch: Epoch) -> f64 {
        self.sigma.powi(2)
    }

    fn sample<R: Rng>(&mut self, _epoch: Epoch, rng: &mut R) -> f64 {
        rng.sample(Normal::new(self.mean, self.sigma).unwrap())
    }
}

impl Mul<f64> for WhiteNoise {
    type Output = Self;

    /// Scale the white noise sigmas by a constant.
    fn mul(mut self, rhs: f64) -> Self::Output {
        self.sigma *= rhs;
        self
    }
}

impl MulAssign<f64> for WhiteNoise {
    fn mul_assign(&mut self, rhs: f64) {
        *self = *self * rhs;
    }
}

#[cfg(test)]
mod ut_wn {
    use hifitime::{Epoch, TimeUnits};
    use rand_pcg::Pcg64Mcg;

    use super::{Stochastics, WhiteNoise};

    #[test]
    fn white_noise_test() {
        let sigma = 10.0_f64;
        let mut wn = WhiteNoise { mean: 0.0, sigma };

        let mut larger_wn = WhiteNoise {
            mean: 0.0,
            sigma: sigma * 10.0,
        };

        let epoch = Epoch::now().unwrap();

        let mut rng = Pcg64Mcg::new(1000);
        let mut cnt_above_3sigma = 0;
        let mut cnt_below_3sigma = 0;
        let mut larger_cnt_above_3sigma = 0;
        let mut larger_cnt_below_3sigma = 0;
        for seconds in 0..1000_i64 {
            let bias = wn.sample(epoch + seconds.seconds(), &mut rng);

            if bias > 3.0 * sigma {
                cnt_above_3sigma += 1;
            } else if bias < -3.0 * sigma {
                cnt_below_3sigma += 1;
            }

            let larger_bias = larger_wn.sample(epoch + seconds.seconds(), &mut rng);
            if larger_bias > 30.0 * sigma {
                larger_cnt_above_3sigma += 1;
            } else if larger_bias < -30.0 * sigma {
                larger_cnt_below_3sigma += 1;
            }
        }

        assert!(dbg!(cnt_above_3sigma) <= 3);
        assert!(dbg!(cnt_below_3sigma) <= 3);

        assert!(dbg!(larger_cnt_above_3sigma) <= 3);
        assert!(dbg!(larger_cnt_below_3sigma) <= 3);
    }
}
