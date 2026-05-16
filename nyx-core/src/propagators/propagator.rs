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

use std::sync::Arc;

use anise::almanac::Almanac;

use super::{IntegrationDetails, IntegratorMethod, IntegratorOptions, PropInstance};
use crate::dynamics::Dynamics;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, OVector};
use crate::time::Duration;
use crate::State;

/// A Propagator allows propagating a set of dynamics forward or backward in time.
/// It is an EventTracker, without any event tracking. It includes the options, the integrator
/// details of the previous step, and the set of coefficients used for the monomorphic instance.
#[derive(Clone, Debug)]
pub struct Propagator<D: Dynamics>
where
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>,
{
    pub dynamics: D, // Stores the dynamics used. *Must* use this to get the latest values
    pub opts: IntegratorOptions, // Stores the integration options (tolerance, min/max step, init step, etc.)
    pub method: IntegratorMethod,
}

/// The `Propagator` trait defines the functions of a propagator and of an event tracker.
impl<D: Dynamics> Propagator<D>
where
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>,
{
    /// Each propagator must be initialized with `new` which stores propagator information.
    pub fn new(dynamics: D, method: IntegratorMethod, opts: IntegratorOptions) -> Self {
        Self {
            dynamics,
            opts,
            method,
        }
    }

    /// Set the tolerance for the propagator
    pub fn set_tolerance(&mut self, tol: f64) {
        self.opts.tolerance = tol;
    }

    /// Set the maximum step size for the propagator and sets the initial step to that value if currently greater
    pub fn set_max_step(&mut self, step: Duration) {
        self.opts.set_max_step(step);
    }

    pub fn set_min_step(&mut self, step: Duration) {
        self.opts.set_min_step(step);
    }

    /// An RK89 propagator (the default) with custom propagator options.
    pub fn rk89(dynamics: D, opts: IntegratorOptions) -> Self {
        Self::new(dynamics, IntegratorMethod::RungeKutta89, opts)
    }

    /// A Dormand Prince 7-8 propagator with custom propagator options: it's about 20% faster than an RK98, and more stable in two body dynamics.
    /// WARNINGS: Dormand Prince may have issues with generating proper trajectories, leading to glitches in event finding.
    pub fn dp78(dynamics: D, opts: IntegratorOptions) -> Self {
        Self::new(dynamics, IntegratorMethod::DormandPrince78, opts)
    }

    pub fn with(&self, state: D::StateType, almanac: Arc<Almanac>) -> PropInstance<'_, D> {
        // Pre-allocate the k used in the propagator
        let mut k = Vec::with_capacity(self.method.stages() + 1);
        for _ in 0..self.method.stages() {
            k.push(OVector::<f64, <D::StateType as State>::VecLength>::zeros());
        }
        PropInstance {
            state,
            prop: self,
            details: IntegrationDetails {
                step: self.opts.init_step,
                error: 0.0,
                attempts: 1,
            },
            log_progress: true,
            almanac,
            step_size: self.opts.init_step,
            fixed_step: self.opts.fixed_step,
            k,
        }
    }

    /// Default propagator is an RK89 with the default PropOpts.
    pub fn default(dynamics: D) -> Self {
        Self::rk89(dynamics, IntegratorOptions::default())
    }

    /// A default Dormand Prince 78 propagator with the default PropOpts.
    /// Faster and more stable than an RK89 (`default`) but seems to cause issues for event finding.
    /// WARNINGS: Dormand Prince may have issues with generating proper trajectories, leading to glitches in event finding.
    pub fn default_dp78(dynamics: D) -> Self {
        Self::dp78(dynamics, IntegratorOptions::default())
    }
}
