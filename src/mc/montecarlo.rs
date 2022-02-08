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
use crate::mc::results::McResults;
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

// TODO: 2. Create a Results type for post processing and a serializer to fit in any serialization (no deserializer)

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
    /// Returns the state found and the trajectory until `max_duration`
    #[must_use = "Monte Carlo result must be used"]
    pub fn run_until_nth_event<F: EventEvaluator<D::StateType>>(
        self,
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
            .take(num_runs)
            .enumerate()
            .collect::<Vec<(usize, D::StateType)>>();

        // And propagate
        init_states.par_iter().progress_with(pb).for_each_with(
            (arc_prop, tx),
            |(arc_prop, sender), (idx, state)| {
                let rslt = arc_prop
                    .with(*state)
                    .until_nth_event(max_duration, event, trigger);
                sender.send((*idx, rslt)).unwrap();
            },
        );

        let clock_time = StdInstant::now() - start;
        println!(
            "Propagated {} states in {}",
            num_runs,
            clock_time.as_secs_f64() * TimeUnit::Second
        );

        McResults {
            data: rx.iter().collect(),
        }
    }

    /// Generate states and propagate each independently until a specific event is found `trigger` times.
    /// Returns the state found and the trajectory until `max_duration`
    #[must_use = "Monte Carlo result must be used"]
    pub fn run_until_epoch(self, end_time: Epoch, num_runs: usize) -> McResults<D::StateType> {
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
            .take(num_runs)
            .enumerate()
            .collect::<Vec<(usize, D::StateType)>>();

        // And propagate on the thread pool
        init_states.par_iter().progress_with(pb).for_each_with(
            (arc_prop, tx),
            |(arc_prop, sender), (idx, state)| {
                let rslt = arc_prop.with(*state).until_epoch_with_traj(end_time);
                sender.send((*idx, rslt)).unwrap();
            },
        );

        let clock_time = StdInstant::now() - start;
        println!(
            "Propagated {} states in {}",
            num_runs,
            clock_time.as_secs_f64() * TimeUnit::Second
        );

        McResults {
            data: rx.iter().collect(),
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

#[test]
fn test_monte_carlo_epoch() {
    use super::GaussianGenerator;
    use crate::cosmic::{Bodies, Cosm, Orbit};
    use crate::dynamics::{OrbitalDynamics, PointMasses};
    use crate::mc::results::Stats;
    use crate::md::StateParameter;
    use crate::time::TimeUnitHelper;
    let cosm = Cosm::de438();
    // Build the state generator

    let eme2k = cosm.frame("EME2000");
    let dt = Epoch::from_gregorian_utc_at_midnight(2021, 1, 31);
    let state = Orbit::keplerian(8_191.93, 1e-6, 12.85, 306.614, 314.19, 99.887_7, dt, eme2k);

    // 5% error on SMA and 5% on Eccentricity
    let generator = GaussianGenerator::from_1σs_prct(
        state,
        &[
            (StateParameter::SMA, 0.05),
            (StateParameter::Eccentricity, 0.05),
        ],
    )
    .unwrap();

    // Set up the dynamics
    let orbital_dyn = OrbitalDynamics::new(vec![PointMasses::new(
        &[Bodies::Sun, Bodies::Luna, Bodies::JupiterBarycenter],
        cosm.clone(),
    )]);

    let prop = Propagator::default(orbital_dyn);

    let my_mc = MonteCarlo {
        prop,
        generator,
        seed: 0,
        scenario: "demo".to_string(),
    };

    let rslts = my_mc.run_until_epoch(dt + 1.0_f64 * TimeUnit::Day, 100);

    let average_sma = rslts
        .report_every(StateParameter::SMA, 5_i32.minutes(), None)
        .amean()
        .unwrap();

    println!("Average SMA = {} km", average_sma);
}
