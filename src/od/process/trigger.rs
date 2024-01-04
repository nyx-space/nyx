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

use hifitime::Epoch;

use super::estimate::Estimate;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::time::Duration;
use crate::State;

/// An EkfTrigger on the number of measurements processed and a time between measurements.
pub struct EkfTrigger {
    pub num_msrs: usize,
    /// In seconds!
    pub disable_time: Duration,
    /// Set to the sigma number needed to switch to the EKF (cf. 68–95–99.7 rule). If number is negative, this is ignored.
    pub within_sigma: f64,
    prev_msr_dt: Option<Epoch>,
    cur_msrs: usize,
    inhibit: bool,
}

impl EkfTrigger {
    pub fn new(num_msrs: usize, disable_time: Duration) -> Self {
        Self {
            num_msrs,
            disable_time,
            within_sigma: -1.0,
            prev_msr_dt: None,
            cur_msrs: 0,
            inhibit: false,
        }
    }

    pub fn enable_ekf<T: State, E>(&mut self, est: &E) -> bool
    where
        E: Estimate<T>,
        DefaultAllocator: Allocator<f64, <T as State>::Size>
            + Allocator<f64, <T as State>::VecLength>
            + Allocator<f64, <T as State>::Size, <T as State>::Size>,
    {
        if self.inhibit {
            return false;
        }

        if !est.predicted() {
            // If this isn't a prediction, let's update the previous measurement time
            self.prev_msr_dt = Some(est.epoch());
        }
        self.cur_msrs += 1;
        self.cur_msrs >= self.num_msrs
            && ((self.within_sigma > 0.0 && est.within_sigma(self.within_sigma))
                || self.within_sigma <= 0.0)
    }

    pub fn disable_ekf(&mut self, epoch: Epoch) -> bool {
        if self.inhibit {
            return true;
        }
        // Return true if there is a prev msr dt, and the next measurement time is more than the disable time seconds away
        match self.prev_msr_dt {
            Some(prev_dt) => {
                if (epoch - prev_dt).abs() > self.disable_time {
                    self.cur_msrs = 0;
                    true
                } else {
                    false
                }
            }
            None => false,
        }
    }

    pub fn set_inhibit(&mut self, inhibit: bool) {
        self.inhibit = inhibit;
    }

    pub fn reset(&mut self) {
        self.cur_msrs = 0;
        self.prev_msr_dt = None;
    }
}
