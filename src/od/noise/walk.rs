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

use hifitime::{Duration, Epoch};
use rand::Rng;
use rand_distr::Normal;
use serde_derive::{Deserialize, Serialize};

use super::Stochastics;

/// Random walk is the simplest non-constant systematic random process.
///
/// # Implementation
/// This implementation matches section 5.2.1 of the NASA Best Practices for Navigation Filters.
/// The variance is computed as p_b(t + Δt) = p_b(t) + qΔt , where q is the process noise.
#[derive(Copy, Clone, Debug, Default, PartialEq, Serialize, Deserialize)]
pub struct RandomWalk {
    /// Process noise in units of measurement units per second.
    /// E.g. if the measurement kind is range, this should be in kilometers per second.
    pub process_noise: f64,
    /// Epoch of the previous realization, used to compute the time delta for the process noise.
    #[serde(skip)]
    pub prev_epoch: Option<Epoch>,
    /// Epoch of the previous realization, used to compute the time delta for the process noise.
    #[serde(skip)]
    pub prev_sample: Option<f64>,
    /// Variance of the previous realization
    #[serde(skip)]
    pub init_variance: Option<f64>,
}

impl RandomWalk {
    /// Initializes a new random walk stochastic noise model from the process noise and the integration time.
    /// This will compute the process noise per second automatically.
    pub fn new(process_noise: f64, integration_time: Duration) -> Self {
        Self {
            process_noise: process_noise / integration_time.to_seconds(),
            ..Default::default()
        }
    }

    /// Initializes a new random walk stochastic noise model from the provided process noise, assuming that the noise level
    /// is fixed regardless of the integration time.
    pub fn constant_white_noise(process_noise: f64) -> Self {
        Self {
            process_noise,
            ..Default::default()
        }
    }
}

impl Stochastics for RandomWalk {
    fn variance(&self, epoch: Epoch) -> f64 {
        // TODO: Figure out how we can have this variance grow _a bit_.
        if let Some(prev_epoch) = self.prev_epoch {
            let delta_t = (epoch - prev_epoch).to_seconds();
            self.process_noise * delta_t
        } else {
            self.process_noise
        }
    }

    fn sample<R: Rng>(&mut self, epoch: Epoch, rng: &mut R) -> f64 {
        // Compute the delta time in seconds between the previous epoch and the sample epoch.
        let dt_s = (match self.prev_epoch {
            None => Duration::ZERO,
            Some(prev_epoch) => epoch - prev_epoch,
        })
        .to_seconds();
        self.prev_epoch = Some(epoch);

        if self.prev_sample.is_none() {
            // Perform the initial sample.
            self.prev_sample = Some(rng.sample(Normal::new(0.0, self.process_noise).unwrap()));
        } else {
            // Compute the next sample otherwise.
            self.prev_sample = Some(
                self.prev_sample.unwrap()
                    + rng.sample(Normal::new(0.0, self.process_noise * dt_s).unwrap()),
            );
        }

        self.prev_sample.unwrap()
    }
}

#[cfg(test)]
mod ut_walk {
    use hifitime::{Epoch, TimeUnits};
    use rand_pcg::Pcg64Mcg;

    use crate::od::noise::walk::{RandomWalk, Stochastics};

    #[test]
    fn white_noise_test() {
        let sigma = 10.0_f64;
        let mut walker = RandomWalk {
            process_noise: sigma.sqrt(),
            ..Default::default()
        };

        let mut larger_walker = RandomWalk {
            process_noise: sigma.sqrt() * 10.0,
            ..Default::default()
        };

        let epoch = Epoch::now().unwrap();

        let mut rng = Pcg64Mcg::new(1000);
        let mut cnt_above_3sigma = 0;
        let mut cnt_below_3sigma = 0;
        let mut larger_cnt_above_3sigma = 0;
        let mut larger_cnt_below_3sigma = 0;
        for seconds in 0..1000_i64 {
            let bias = walker.sample(epoch + seconds.seconds(), &mut rng);

            if bias > 3.0 * sigma {
                cnt_above_3sigma += 1;
            } else if bias < -3.0 * sigma {
                cnt_below_3sigma += 1;
            }

            let larger_bias = larger_walker.sample(epoch + seconds.seconds(), &mut rng);
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
