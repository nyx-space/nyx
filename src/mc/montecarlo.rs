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
use super::rand_distr::Distribution;
use super::rayon::prelude::ParallelIterator;
use super::{Generator, Pcg64Mcg};
use crate::dynamics::Dynamics;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::mc::results::{McResults, Run};
use crate::mc::DispersedState;
use crate::md::trajectory::InterpState;
use crate::md::EventEvaluator;
use crate::propagators::{ErrorCtrl, Propagator};
use crate::time::{Duration, Epoch, TimeUnit};
use crate::State;
use indicatif::{ParallelProgressIterator, ProgressBar, ProgressStyle};
use rayon::prelude::*;
use std::f64;
use std::fmt;
use std::sync::mpsc::channel;
use std::sync::Arc;
use std::time::Instant as StdInstant;

/// A Monte Carlo framework, automatically running on all threads via a thread pool. This framework is targeted toward analysis of time-continuous variables.
/// One caveat of the design is that the trajectory is used for post processing, not each individual state. This may prevent some event switching from being shown in GNC simulations.
pub struct MonteCarlo<'a, D: Dynamics, E: ErrorCtrl, Distr: Distribution<f64> + Copy>
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
    /// Generator of states for the Monte Carlo run
    pub generator: Generator<D::StateType, Distr>,
    /// Name of this run, will be reflected in the progress bar and in the output structure
    pub scenario: String,
}

impl<'a, D: Dynamics, E: ErrorCtrl, Distr: Distribution<f64> + Copy> MonteCarlo<'a, D, E, Distr>
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
    pub fn run_until_nth_event<F: EventEvaluator<D::StateType>>(
        self,
        max_duration: Duration,
        event: &F,
        trigger: usize,
        num_runs: usize,
    ) -> McResults<D::StateType> {
        self.resume_run_until_nth_event(0, max_duration, event, trigger, num_runs)
    }

    /// Generate states and propagate each independently until a specific event is found `trigger` times.
    #[must_use = "Monte Carlo result must be used"]
    pub fn resume_run_until_nth_event<F: EventEvaluator<D::StateType>>(
        self,
        skip: usize,
        max_duration: Duration,
        event: &F,
        trigger: usize,
        num_runs: usize,
    ) -> McResults<D::StateType> {
        // Setup the RNG
        let rng = Pcg64Mcg::new(self.seed.into());
        // Setup the progress bar
        let pb = self.progress_bar(num_runs);
        // Wrap the propagator in an Arc
        let arc_prop = Arc::new(self.prop);

        // Setup the channels
        let (tx, rx) = channel();

        // Generate all states (must be done separately because the rng is not thread safe)
        let start = StdInstant::now();

        let init_states = self
            .generator
            .sample_iter(rng)
            .skip(skip)
            .take(num_runs)
            .enumerate()
            .collect::<Vec<(usize, DispersedState<D::StateType>)>>();

        // And propagate
        init_states.par_iter().progress_with(pb).for_each_with(
            (arc_prop, tx),
            |(arc_prop, sender), (index, dispersed_state)| {
                let result = arc_prop.with(dispersed_state.state).until_nth_event(
                    max_duration,
                    event,
                    trigger,
                );

                // Build a single run result
                let run = Run {
                    index: *index,
                    dispersed_state: dispersed_state.clone(),
                    result,
                };
                sender.send(run).unwrap();
            },
        );

        let clock_time = StdInstant::now() - start;
        println!(
            "Propagated {} states in {}",
            num_runs,
            clock_time.as_secs_f64() * TimeUnit::Second
        );

        // Collect all of the results and sort them by run index
        let mut runs = rx.iter().collect::<Vec<Run<D::StateType>>>();
        runs.par_sort_by_key(|run| run.index);

        McResults {
            runs,
            scenario: self.scenario,
        }
    }

    /// Generate states and propagate each independently until a specific event is found `trigger` times.
    #[must_use = "Monte Carlo result must be used"]
    pub fn run_until_epoch(self, end_epoch: Epoch, num_runs: usize) -> McResults<D::StateType> {
        self.resume_run_until_epoch(0, end_epoch, num_runs)
    }

    /// Resumes a Monte Carlo run by skipping the first `skip` items, generating states only after that, and propagate each independently until the specified epoch.
    #[must_use = "Monte Carlo result must be used"]
    pub fn resume_run_until_epoch(
        self,
        skip: usize,
        end_epoch: Epoch,
        num_runs: usize,
    ) -> McResults<D::StateType> {
        // Setup the RNG
        let rng = Pcg64Mcg::new(self.seed.into());
        // Setup the progress bar
        let pb = self.progress_bar(num_runs);
        // Wrap the propagator in an Arc
        let arc_prop = Arc::new(self.prop);
        // Setup the channels
        let (tx, rx) = channel();

        // Generate all states (must be done separately because the rng is not thread safe)
        let start = StdInstant::now();

        let init_states = self
            .generator
            .sample_iter(rng)
            .skip(skip)
            .take(num_runs)
            .enumerate()
            .collect::<Vec<(usize, DispersedState<D::StateType>)>>();

        // And propagate on the thread pool
        init_states.par_iter().progress_with(pb).for_each_with(
            (arc_prop, tx),
            |(arc_prop, sender), (index, dispersed_state)| {
                let result = arc_prop
                    .with(dispersed_state.state)
                    .until_epoch_with_traj(end_epoch);

                // Build a single run result
                let run = Run {
                    index: *index,
                    dispersed_state: dispersed_state.clone(),
                    result,
                };
                sender.send(run).unwrap();
            },
        );

        let clock_time = StdInstant::now() - start;
        println!(
            "Propagated {} states in {}",
            num_runs,
            clock_time.as_secs_f64() * TimeUnit::Second
        );

        // Collect all of the results and sort them by run index
        let mut runs = rx.iter().collect::<Vec<Run<D::StateType>>>();
        runs.par_sort_by_key(|run| run.index);

        McResults {
            runs,
            scenario: self.scenario,
        }
    }
}

impl<'a, D: Dynamics, E: ErrorCtrl, Distr: Distribution<f64> + Copy> fmt::Display
    for MonteCarlo<'a, D, E, Distr>
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

impl<'a, D: Dynamics, E: ErrorCtrl, Distr: Distribution<f64> + Copy> fmt::LowerHex
    for MonteCarlo<'a, D, E, Distr>
where
    D::StateType: InterpState,
    DefaultAllocator: Allocator<f64, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<usize, <D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<f64, <D::StateType as State>::VecLength>,
    <DefaultAllocator as Allocator<f64, <D::StateType as State>::VecLength>>::Buffer: Send,
{
    /// Returns a filename friendly name
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "mc-data-{}-seed-{}",
            self.scenario.replace(" ", "-"),
            self.seed
        )
    }
}
