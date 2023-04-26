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

use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;

pub use crate::od::estimate::*;
pub use crate::od::kalman::*;
pub use crate::od::ground_station::*;
pub use crate::od::residual::*;
pub use crate::od::snc::*;
pub use crate::od::*;

pub use crate::time::{Duration, Unit};
use crate::State;

use self::prelude::IterationConf;

/// A trait detailing when to switch to from a CKF to an EKF
pub trait KfTrigger {
    fn enable_ekf<T: State, E>(&mut self, est: &E) -> bool
    where
        E: Estimate<T>,
        DefaultAllocator: Allocator<f64, <T as State>::Size>
            + Allocator<f64, <T as State>::VecLength>
            + Allocator<f64, <T as State>::Size, <T as State>::Size>;

    /// Return true if the filter should not longer be as extended.
    /// By default, this returns false, i.e. when a filter has been switched to an EKF, it will
    /// remain as such.
    fn disable_ekf(&mut self, _epoch: Epoch) -> bool {
        false
    }

    /// If some iteration configuration is returned, the filter will iterate with it before enabling the EKF.
    fn interation_config(&self) -> Option<IterationConf> {
        None
    }
}

/// CkfTrigger will never switch a KF to an EKF
pub struct CkfTrigger;

impl KfTrigger for CkfTrigger {
    fn enable_ekf<T: State, E>(&mut self, _est: &E) -> bool
    where
        E: Estimate<T>,
        DefaultAllocator: Allocator<f64, <T as State>::Size>
            + Allocator<f64, <T as State>::VecLength>
            + Allocator<f64, <T as State>::Size, <T as State>::Size>,
    {
        false
    }
}

/// An EkfTrigger on the number of measurements processed and a time between measurements.
pub struct EkfTrigger {
    pub num_msrs: usize,
    /// In seconds!
    pub disable_time: Duration,
    /// Set to the sigma number needed to switch to the EKF (cf. 68–95–99.7 rule). If number is negative, this is ignored.
    pub within_sigma: f64,
    prev_msr_dt: Option<Epoch>,
    cur_msrs: usize,
}

impl EkfTrigger {
    pub fn new(num_msrs: usize, disable_time: Duration) -> Self {
        Self {
            num_msrs,
            disable_time,
            within_sigma: -1.0,
            prev_msr_dt: None,
            cur_msrs: 0,
        }
    }
}

impl KfTrigger for EkfTrigger {
    fn enable_ekf<T: State, E>(&mut self, est: &E) -> bool
    where
        E: Estimate<T>,
        DefaultAllocator: Allocator<f64, <T as State>::Size>
            + Allocator<f64, <T as State>::VecLength>
            + Allocator<f64, <T as State>::Size, <T as State>::Size>,
    {
        if !est.predicted() {
            // If this isn't a prediction, let's update the previous measurement time
            self.prev_msr_dt = Some(est.epoch());
        }
        self.cur_msrs += 1;
        self.cur_msrs >= self.num_msrs
            && ((self.within_sigma > 0.0 && est.within_sigma(self.within_sigma))
                || self.within_sigma <= 0.0)
    }

    fn disable_ekf(&mut self, epoch: Epoch) -> bool {
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
}
