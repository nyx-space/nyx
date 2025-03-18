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

use super::Pcg64Mcg;
use crate::dynamics::Dynamics;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::mc::results::{PropResult, Results, Run};
use crate::mc::DispersedState;
use crate::md::trajectory::Interpolatable;
use crate::md::EventEvaluator;
use crate::propagators::Propagator;
#[cfg(not(target_arch = "wasm32"))]
use crate::time::Unit;
use crate::time::{Duration, Epoch};
use crate::State;
use anise::almanac::Almanac;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use log::info;
use rand::SeedableRng;
use rand_distr::Distribution;
use rayon::prelude::ParallelIterator;
use rayon::prelude::*;
use std::fmt;
use std::sync::mpsc::channel;
use std::sync::Arc;
#[cfg(not(target_arch = "wasm32"))]
use std::time::Instant as StdInstant;

/// A Monte Carlo framework, automatically running on all threads via a thread pool. This framework is targeted toward analysis of time-continuous variables.
/// One caveat of the design is that the trajectory is used for post processing, not each individual state. This may prevent some event switching from being shown in GNC simulations.
pub struct MonteCarlo<S: Interpolatable, Distr: Distribution<DispersedState<S>>>
where
    DefaultAllocator: Allocator<S::Size> + Allocator<S::Size, S::Size> + Allocator<S::VecLength>,
{
    /// Seed of the [64bit PCG random number generator](https://www.pcg-random.org/index.html)
    pub seed: Option<u128>,
    /// Generator of states for the Monte Carlo run
    pub random_state: Distr,
    /// Name of this run, will be reflected in the progress bar and in the output structure
    pub scenario: String,
    pub nominal_state: S,
}

impl<S: Interpolatable, Distr: Distribution<DispersedState<S>>> MonteCarlo<S, Distr>
where
    DefaultAllocator: Allocator<S::Size> + Allocator<S::Size, S::Size> + Allocator<S::VecLength>,
{
    pub fn new(
        nominal_state: S,
        random_variable: Distr,
        scenario: String,
        seed: Option<u128>,
    ) -> Self {
        Self {
            random_state: random_variable,
            seed,
            scenario,
            nominal_state,
        }
    }
    // Just the template for the progress bar
    fn progress_bar(&self, num_runs: usize) -> ProgressBar {
        let pb = ProgressBar::new(num_runs.try_into().unwrap());
        pb.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:100.cyan/blue} {pos:>7}/{len:7} {msg}")
                .unwrap()
                .progress_chars("##-"),
        );
        pb.set_message(format!("{self}"));
        pb
    }

    /// Generate states and propagate each independently until a specific event is found `trigger` times.
    #[allow(clippy::needless_lifetimes)]
    pub fn run_until_nth_event<D, F>(
        self,
        prop: Propagator<D>,
        almanac: Arc<Almanac>,
        max_duration: Duration,
        event: &F,
        trigger: usize,
        num_runs: usize,
    ) -> Results<S, PropResult<S>>
    where
        D: Dynamics<StateType = S>,

        F: EventEvaluator<S>,
        DefaultAllocator: Allocator<<D::StateType as State>::Size>
            + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
            + Allocator<<D::StateType as State>::VecLength>,
        <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    {
        self.resume_run_until_nth_event(prop, almanac, 0, max_duration, event, trigger, num_runs)
    }

    /// Generate states and propagate each independently until a specific event is found `trigger` times.
    #[must_use = "Monte Carlo result must be used"]
    #[allow(clippy::needless_lifetimes)]
    pub fn resume_run_until_nth_event<D, F>(
        &self,
        prop: Propagator<D>,
        almanac: Arc<Almanac>,
        skip: usize,
        max_duration: Duration,
        event: &F,
        trigger: usize,
        num_runs: usize,
    ) -> Results<S, PropResult<S>>
    where
        D: Dynamics<StateType = S>,

        F: EventEvaluator<S>,
        DefaultAllocator: Allocator<<D::StateType as State>::Size>
            + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
            + Allocator<<D::StateType as State>::VecLength>,
        <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    {
        // Generate the initial states
        let init_states = self.generate_states(skip, num_runs, self.seed);
        // Setup the progress bar
        let pb = self.progress_bar(num_runs);
        // Setup the thread friendly communication
        let (tx, rx) = channel();

        // Generate all states (must be done separately because the rng is not thread safe)
        #[cfg(not(target_arch = "wasm32"))]
        let start = StdInstant::now();

        init_states.par_iter().progress_with(pb).for_each_with(
            (prop, tx),
            |(prop, tx), (index, dispersed_state)| {
                let result = prop
                    .with(dispersed_state.state, almanac.clone())
                    .until_nth_event(max_duration, event, trigger);

                // Build a single run result
                let run = Run {
                    index: *index,
                    dispersed_state: dispersed_state.clone(),
                    result: result.map(|r| PropResult {
                        state: r.0,
                        traj: r.1,
                    }),
                };
                tx.send(run).unwrap();
            },
        );

        #[cfg(not(target_arch = "wasm32"))]
        {
            let clock_time = StdInstant::now() - start;
            info!(
                "Propagated {} states in {}",
                num_runs,
                clock_time.as_secs_f64() * Unit::Second
            );
        }

        // Collect all of the results and sort them by run index
        let mut runs = rx
            .iter()
            .collect::<Vec<Run<D::StateType, PropResult<D::StateType>>>>();
        runs.par_sort_by_key(|run| run.index);

        Results {
            runs,
            scenario: self.scenario.clone(),
        }
    }

    /// Generate states and propagate each independently until a specific event is found `trigger` times.
    #[must_use = "Monte Carlo result must be used"]
    #[allow(clippy::needless_lifetimes)]
    pub fn run_until_epoch<D>(
        self,
        prop: Propagator<D>,
        almanac: Arc<Almanac>,
        end_epoch: Epoch,
        num_runs: usize,
    ) -> Results<S, PropResult<S>>
    where
        D: Dynamics<StateType = S>,

        DefaultAllocator: Allocator<<D::StateType as State>::Size>
            + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
            + Allocator<<D::StateType as State>::VecLength>,
        <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    {
        self.resume_run_until_epoch(prop, almanac, 0, end_epoch, num_runs)
    }

    /// Resumes a Monte Carlo run by skipping the first `skip` items, generating states only after that, and propagate each independently until the specified epoch.
    #[must_use = "Monte Carlo result must be used"]
    #[allow(clippy::needless_lifetimes)]
    pub fn resume_run_until_epoch<D>(
        &self,
        prop: Propagator<D>,
        almanac: Arc<Almanac>,
        skip: usize,
        end_epoch: Epoch,
        num_runs: usize,
    ) -> Results<S, PropResult<S>>
    where
        D: Dynamics<StateType = S>,

        DefaultAllocator: Allocator<<D::StateType as State>::Size>
            + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
            + Allocator<<D::StateType as State>::VecLength>,
        <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    {
        // Generate the initial states
        let init_states = self.generate_states(skip, num_runs, self.seed);
        // Setup the progress bar
        let pb = self.progress_bar(num_runs);
        // Setup the thread friendly communication
        let (tx, rx) = channel();

        // And propagate on the thread pool
        #[cfg(not(target_arch = "wasm32"))]
        let start = StdInstant::now();
        init_states.par_iter().progress_with(pb).for_each_with(
            (prop, tx),
            |(arc_prop, tx), (index, dispersed_state)| {
                let result = arc_prop
                    .with(dispersed_state.state, almanac.clone())
                    .quiet()
                    .until_epoch_with_traj(end_epoch);

                // Build a single run result
                let run = Run {
                    index: *index,
                    dispersed_state: dispersed_state.clone(),
                    result: result.map(|r| PropResult {
                        state: r.0,
                        traj: r.1,
                    }),
                };

                tx.send(run).unwrap();
            },
        );

        #[cfg(not(target_arch = "wasm32"))]
        {
            let clock_time = StdInstant::now() - start;
            info!(
                "Propagated {} states in {}",
                num_runs,
                clock_time.as_secs_f64() * Unit::Second
            );
        }

        // Collect all of the results and sort them by run index
        let mut runs = rx.iter().collect::<Vec<Run<S, PropResult<S>>>>();
        runs.par_sort_by_key(|run| run.index);

        Results {
            runs,
            scenario: self.scenario.clone(),
        }
    }

    /// Set up the seed and generate the states. This is useful for checking the generated states before running a large scale Monte Carlo.
    #[must_use = "Generated states for a Monte Carlo run must be used"]
    pub fn generate_states(
        &self,
        skip: usize,
        num_runs: usize,
        seed: Option<u128>,
    ) -> Vec<(usize, DispersedState<S>)> {
        // Setup the RNG
        let rng = match seed {
            Some(seed) => Pcg64Mcg::new(seed),
            None => Pcg64Mcg::from_os_rng(),
        };

        // Generate the states, forcing the borrow as specified in the `sample_iter` docs.
        (&self.random_state)
            .sample_iter(rng)
            .skip(skip)
            .take(num_runs)
            .enumerate()
            .collect::<Vec<(usize, DispersedState<S>)>>()
    }
}

impl<S: Interpolatable, Distr: Distribution<DispersedState<S>>> fmt::Display
    for MonteCarlo<S, Distr>
where
    DefaultAllocator: Allocator<S::Size> + Allocator<S::Size, S::Size> + Allocator<S::VecLength>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{} - Nyx Monte Carlo - seed: {:?}",
            self.scenario, self.seed
        )
    }
}

impl<S: Interpolatable, Distr: Distribution<DispersedState<S>>> fmt::LowerHex
    for MonteCarlo<S, Distr>
where
    DefaultAllocator: Allocator<S::Size> + Allocator<S::Size, S::Size> + Allocator<S::VecLength>,
{
    /// Returns a filename friendly name
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "mc-data-{}-seed-{:?}",
            self.scenario.replace(' ', "-"),
            self.seed
        )
    }
}
