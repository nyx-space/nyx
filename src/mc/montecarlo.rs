/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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
extern crate indicatif;
extern crate rand;
extern crate rand_distr;
use super::rayon::iter::ParallelBridge;
use super::rayon::prelude::ParallelIterator;
use super::Pcg64Mcg;
use crate::dynamics::Dynamics;
use crate::errors::NyxError;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, OVector};
use crate::md::trajectory::spline::INTERPOLATION_SAMPLES;
use crate::md::trajectory::{interpolate, InterpState, Traj, TrajError};
use crate::md::EventEvaluator;
use crate::propagators::{ErrorCtrl, Propagator};
use crate::time::{Duration, Epoch, TimeUnit};
use crate::{cosmic::Cosm, Orbit, State};
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rand::thread_rng;
use rand_distr::{Distribution, Normal};
use rayon::prelude::*;
use std::collections::BTreeMap;
use std::f64;
use std::fmt;
use std::sync::mpsc::{channel, Sender};
use std::sync::Arc;
use std::time::Instant as StdInstant;

// TODO: 1. Type with state generator; 2. Create a Results type for post processing and a serializer (no deserializer)

pub struct MonteCarlo<'a, D: Dynamics, E: ErrorCtrl>
where
    D::StateType: InterpState,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
    <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
{
    /// Propagator with dynamics to use in all runs
    pub prop: Propagator<'a, D, E>,
    /// Seed of the [64bit PCG random number generator](https://www.pcg-random.org/index.html)
    pub seed: u64,
    /// Set to True to disperse the initial state provided in the run
    pub disperse_template: bool,
    pub scenario: String,
}

impl<'a, D: Dynamics, E: ErrorCtrl> MonteCarlo<'a, D, E>
where
    D::StateType: InterpState,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
    <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
{
    // Just the template for the progress bar
    fn progress_bar(&self, num_runs: usize) -> ProgressBar {
        let pb = ProgressBar::new(num_runs.try_into().unwrap());
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:100.cyan/blue} {pos:>7}/{len:7} {msg}")
                .progress_chars("##-"),
        );
        pb.set_message(format!("{}", self));
        pb
    }

    /// Generate states and propagate each independently until a specific event is found `trigger` times.
    /// Returns the state found and the trajectory until `max_duration`
    pub fn run_until_nth_event<F: EventEvaluator<D::StateType>>(
        self,
        template: D::StateType,
        max_duration: Duration,
        event: &F,
        trigger: usize,
        num_runs: usize,
    ) {
        let cosm = Cosm::de438();
        let eme2k = cosm.frame("EME2000");
        let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
        let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

        // Around 1 km of error
        let sma_dist = Normal::new(0.0, 1.0).unwrap();

        let start = StdInstant::now();

        // Generate all 100 initial states
        let init_states: Vec<Orbit> = sma_dist
            .sample_iter(&mut thread_rng())
            .take(num_runs)
            .map(|delta_sma| state.with_sma(state.sma() + delta_sma))
            .collect();

        let pb = self.progress_bar(num_runs);
        let arc_prop = Arc::new(self.prop);

        init_states
            .par_iter()
            .progress_with(pb)
            .for_each_with(arc_prop, |setup, state| {
                let final_state =
                    setup
                        .with(template)
                        .until_nth_event(max_duration, event, trigger);
            });

        let clock_time = StdInstant::now() - start;
        println!(
            "Propagated {} states in {} seconds",
            init_states.len(),
            clock_time.as_secs_f64()
        );
    }

    /// Generate states and propagate each independently until a specific event is found `trigger` times.
    /// Returns the state found and the trajectory until `max_duration`
    pub fn run_until_epoch(self, template: D::StateType, end_time: Epoch, num_runs: usize) {
        let cosm = Cosm::de438();
        let eme2k = cosm.frame("EME2000");
        let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
        let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

        // Around 1 km of error
        let sma_dist = Normal::new(0.0, 1.0).unwrap();

        let start = StdInstant::now();

        // Generate all 100 initial states
        let init_states: Vec<Orbit> = sma_dist
            .sample_iter(&mut thread_rng())
            .take(num_runs)
            .map(|delta_sma| state.with_sma(state.sma() + delta_sma))
            .collect();

        let pb = self.progress_bar(num_runs);
        let arc_prop = Arc::new(self.prop);

        init_states
            .par_iter()
            .progress_with(pb)
            .for_each_with(arc_prop, |setup, state| {
                let final_state = setup.with(template).until_epoch(end_time);
            });

        let clock_time = StdInstant::now() - start;
        println!(
            "Propagated {} states in {} seconds",
            init_states.len(),
            clock_time.as_secs_f64()
        );
    }
}

impl<'a, D: Dynamics, E: ErrorCtrl> fmt::Display for MonteCarlo<'a, D, E>
where
    D::StateType: InterpState,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
    <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} - Nyx Monte Carlo - seed: {}",
            self.scenario, self.seed
        )
    }
}

#[test]
fn test_mc_setup() {
    use crate::cosmic::{Bodies, Cosm, Orbit};
    use crate::dynamics::{OrbitalDynamics, PointMasses};
    let cosm = Cosm::de438();
    let orbital_dyn = OrbitalDynamics::new(vec![PointMasses::new(
        &[Bodies::Sun, Bodies::Luna, Bodies::JupiterBarycenter],
        cosm.clone(),
    )]);

    // We need to wrap the propagator setup in an Arc to enable multithreading.
    let prop = Propagator::default(orbital_dyn);

    let my_mc = MonteCarlo {
        prop,
        seed: 0,
        disperse_template: false,
        scenario: "demo".to_string(),
    };

    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
    let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

    my_mc.run_until_epoch(state, dt + 1.0_f64 * TimeUnit::Day, 100);
}
