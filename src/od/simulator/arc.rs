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

use anise::almanac::Almanac;
use hifitime::{Duration, Epoch, TimeSeries, TimeUnits};
use log::{info, warn};
use num::integer::gcd;
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;

use crate::dynamics::NyxError;
use crate::io::ConfigError;
use crate::md::trajectory::Interpolatable;
use crate::od::msr::TrackingDataArc;
use crate::od::prelude::Strand;
use crate::od::simulator::Cadence;
use crate::od::GroundStation;
use crate::Spacecraft;
use crate::State;
use crate::{linalg::allocator::Allocator, od::TrackingDevice};
use crate::{linalg::DefaultAllocator, md::prelude::Traj};
use std::collections::BTreeMap;
use std::fmt::Display;
use std::marker::PhantomData;
use std::sync::Arc;

use super::{Handoff, TrkConfig};

#[derive(Clone)]
pub struct TrackingArcSim<MsrIn, D>
where
    D: TrackingDevice<MsrIn>,
    MsrIn: State,
    MsrIn: Interpolatable,
    DefaultAllocator: Allocator<<MsrIn as State>::Size>
        + Allocator<<MsrIn as State>::Size, <MsrIn as State>::Size>
        + Allocator<<MsrIn as State>::VecLength>,
{
    /// Map of devices from their names.
    pub devices: BTreeMap<String, D>,
    /// Receiver trajectory
    pub trajectory: Traj<MsrIn>,
    /// Configuration of each device
    pub configs: BTreeMap<String, TrkConfig>,
    /// Random number generator used for this tracking arc, ensures repeatability
    rng: Pcg64Mcg,
    /// Greatest common denominator time series that allows this arc to meet all of the conditions.
    time_series: TimeSeries,
    _msr_in: PhantomData<MsrIn>,
}

impl<MsrIn, D> TrackingArcSim<MsrIn, D>
where
    D: TrackingDevice<MsrIn>,
    MsrIn: State,
    MsrIn: Interpolatable,
    DefaultAllocator: Allocator<<MsrIn as State>::Size>
        + Allocator<<MsrIn as State>::Size, <MsrIn as State>::Size>
        + Allocator<<MsrIn as State>::VecLength>,
{
    /// Build a new tracking arc simulator using the provided seeded random number generator.
    pub fn with_rng(
        devices: BTreeMap<String, D>,
        trajectory: Traj<MsrIn>,
        configs: BTreeMap<String, TrkConfig>,
        rng: Pcg64Mcg,
    ) -> Result<Self, ConfigError> {
        // Check that each device has an associated configurations.
        // We don't care if there are more configurations than chosen devices.
        let mut sampling_rates_ns = Vec::with_capacity(devices.len());
        for name in devices.keys() {
            if let Some(cfg) = configs.get(name) {
                if let Err(e) = cfg.sanity_check() {
                    warn!("Ignoring device {name}: {e}");
                    continue;
                }
                sampling_rates_ns.push(cfg.sampling.truncated_nanoseconds());
            } else {
                warn!("Ignoring device {name}: no associated tracking configuration",);
                continue;
            }
        }

        if sampling_rates_ns.is_empty() {
            return Err(ConfigError::InvalidConfig {
                msg: "None of the devices are properly configured".to_string(),
            });
        }

        let common_sampling_rate_ns = sampling_rates_ns
            .iter()
            .fold(sampling_rates_ns[0], |a, &b| gcd(a, b));

        // The overall time series is the one going from the start to the end of the trajectory with the smallest time step
        // of all the tracking configurations.
        let time_series = TimeSeries::inclusive(
            trajectory.first().epoch(),
            trajectory.last().epoch(),
            Duration::from_truncated_nanoseconds(common_sampling_rate_ns),
        );

        let me = Self {
            devices,
            trajectory,
            configs,
            rng,
            time_series,
            _msr_in: PhantomData,
        };

        info!("{me}");

        Ok(me)
    }

    /// Build a new tracking arc simulator using the provided seed to initialize the random number generator.
    pub fn with_seed(
        devices: BTreeMap<String, D>,
        trajectory: Traj<MsrIn>,
        configs: BTreeMap<String, TrkConfig>,
        seed: u64,
    ) -> Result<Self, ConfigError> {
        let rng = Pcg64Mcg::new(seed as u128);

        Self::with_rng(devices, trajectory, configs, rng)
    }

    /// Build a new tracking arc simulator using the system entropy to seed the random number generator.
    pub fn new(
        devices: BTreeMap<String, D>,
        trajectory: Traj<MsrIn>,
        configs: BTreeMap<String, TrkConfig>,
    ) -> Result<Self, ConfigError> {
        let rng = Pcg64Mcg::from_os_rng();

        Self::with_rng(devices, trajectory, configs, rng)
    }

    /// Generates measurements for the tracking arc using the defined strands
    ///
    /// # Warning
    /// This function will return an error if any of the devices defines as a scheduler.
    /// You must create the schedule first using `build_schedule` first.
    ///
    /// # Notes
    /// Although mutable, this function may be called several times to generate different measurements.
    ///
    /// # Algorithm
    /// For each tracking device, and for each strand within that device, sample the trajectory at the sample
    /// rate of the tracking device, adding a measurement whenever the spacecraft is visible.
    /// Build the measurements as a vector, ordered chronologically.
    ///
    pub fn generate_measurements(
        &mut self,
        almanac: Arc<Almanac>,
    ) -> Result<TrackingDataArc, NyxError> {
        let mut measurements = BTreeMap::new();

        for (name, device) in self.devices.iter_mut() {
            if let Some(cfg) = self.configs.get(name) {
                if cfg.scheduler.is_some() {
                    if cfg.strands.is_none() {
                        return Err(NyxError::CustomError {
                            msg: format!(
                                "schedule for {name} must be built before generating measurements"
                            ),
                        });
                    } else {
                        warn!("scheduler for {name} is ignored, using the defined tracking strands instead")
                    }
                }

                let init_msr_count = measurements.len();
                let tick = Epoch::now().unwrap();

                match cfg.strands.as_ref() {
                    Some(strands) => {
                        // Strands are defined at this point
                        'strands: for (ii, strand) in strands.iter().enumerate() {
                            // Build the time series for this strand, sampling at the correct rate
                            for epoch in
                                TimeSeries::inclusive(strand.start, strand.end, cfg.sampling)
                            {
                                match device.measure(
                                    epoch,
                                    &self.trajectory,
                                    Some(&mut self.rng),
                                    almanac.clone(),
                                ) {
                                    Ok(msr_opt) => {
                                        if let Some(msr) = msr_opt {
                                            measurements.insert(epoch, msr);
                                        }
                                    }
                                    Err(e) => {
                                        if epoch != strand.end {
                                            warn!(
                                            "Skipping the remaining strand #{ii} ending on {}: {e}",
                                            strand.end
                                        );
                                        }
                                        continue 'strands;
                                    }
                                }
                            }
                        }

                        info!(
                            "Simulated {} measurements for {name} for {} tracking strands in {}",
                            measurements.len() - init_msr_count,
                            strands.len(),
                            (Epoch::now().unwrap() - tick).round(1.0_f64.milliseconds())
                        );
                    }
                    None => {
                        warn!("No tracking strands defined for {name}, skipping");
                    }
                }
            }
        }

        // Build the tracking arc.
        let trk_data = TrackingDataArc {
            measurements,
            source: None,
            moduli: None,
            force_reject: false,
        };

        Ok(trk_data)
    }
}

impl<MsrIn, D> Display for TrackingArcSim<MsrIn, D>
where
    D: TrackingDevice<MsrIn>,
    MsrIn: Interpolatable,
    DefaultAllocator: Allocator<<MsrIn as State>::Size>
        + Allocator<<MsrIn as State>::Size, <MsrIn as State>::Size>
        + Allocator<<MsrIn as State>::VecLength>,
{
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "Tracking Arc Simulator on {} with devices {:?} over {}",
            self.trajectory,
            self.devices.keys(),
            self.time_series
        )
    }
}

// Literally the same as above, but can't make it generic =(
impl TrackingArcSim<Spacecraft, GroundStation> {
    /// Builds the schedule provided the config. Requires the tracker to be a ground station.
    ///
    /// # Algorithm
    ///
    /// 1. For each tracking device:
    /// 2. Find when the vehicle trajectory has an elevation greater or equal to zero, and use that as the first start of the first tracking arc for this station
    /// 3. Find when the vehicle trajectory has an elevation less than zero (i.e. disappears below the horizon), after that initial epoch
    /// 4. Repeat 2, 3 until the end of the trajectory
    /// 5. Build each of these as "tracking strands" for this tracking device.
    /// 6. Organize all of the built tracking strands chronologically.
    /// 7. Iterate through all of the strands:
    ///    7.a. if that tracker is marked as `Greedy` and it ends after the start of the next strand, change the start date of the next strand.
    ///    7.b. if that tracker is marked as `Eager` and it ends after the start of the next strand, change the end date of the current strand.
    pub fn generate_schedule(
        &self,
        almanac: Arc<Almanac>,
    ) -> Result<BTreeMap<String, TrkConfig>, NyxError> {
        // Consider using find_all via the heuristic
        let mut built_cfg = self.configs.clone();

        for (name, device) in self.devices.iter() {
            if let Some(cfg) = self.configs.get(name) {
                if let Some(scheduler) = cfg.scheduler {
                    info!("Building schedule for {name}");
                    if scheduler.handoff == Handoff::Overlap {
                        warn!("Overlapping measurements on {name} is no longer supported on identical epochs.");
                    }
                    built_cfg.get_mut(name).unwrap().scheduler = None;
                    built_cfg.get_mut(name).unwrap().strands = Some(Vec::new());

                    // Convert the trajectory into the ground station frame.
                    let traj = self.trajectory.to_frame(device.frame, almanac.clone())?;

                    match traj.find_arcs(&device, None, almanac.clone()) {
                        Err(_) => info!("No measurements from {name}"),
                        Ok(elevation_arcs) => {
                            for arc in elevation_arcs {
                                let strand_start = arc.rise.state.epoch();
                                let strand_end = arc.fall.state.epoch();

                                if strand_end - strand_start
                                    < cfg.sampling * i64::from(scheduler.min_samples)
                                {
                                    info!(
                                    "Too few samples from {name} opportunity from {strand_start} to {strand_end}, discarding strand",
                                );
                                    continue;
                                }

                                let mut strand_range = Strand {
                                    start: strand_start,
                                    end: strand_end,
                                };

                                // If there is an alignment, apply it
                                if let Some(alignment) = scheduler.sample_alignment {
                                    strand_range.start = strand_range.start.round(alignment);
                                    strand_range.end = strand_range.end.round(alignment);
                                }

                                if let Cadence::Intermittent { on, off } = scheduler.cadence {
                                    // Check that the next start time is within the allocated time
                                    if let Some(prev_strand) =
                                        built_cfg[name].strands.as_ref().unwrap().last()
                                    {
                                        if prev_strand.end + off > strand_range.start {
                                            // We're turning on the tracking sooner than the schedule allows, so let's fix that.
                                            strand_range.start = prev_strand.end + off;
                                            // Check that we didn't eat into the whole tracking opportunity
                                            if strand_range.start > strand_end {
                                                // Lost this whole opportunity.
                                                info!("Discarding {name} opportunity from {strand_start} to {strand_end} due to cadence {:?}", scheduler.cadence);
                                                continue;
                                            }
                                        }
                                    }
                                    // Check that we aren't tracking for longer than configured
                                    if strand_range.end - strand_range.start > on {
                                        strand_range.end = strand_range.start + on;
                                    }
                                }

                                // We've found when the spacecraft is below the horizon, so this is a new strand.
                                built_cfg
                                    .get_mut(name)
                                    .unwrap()
                                    .strands
                                    .as_mut()
                                    .unwrap()
                                    .push(strand_range);
                            }

                            info!(
                                "Built {} tracking strands for {name}",
                                built_cfg[name].strands.as_ref().unwrap().len()
                            );
                        }
                    }
                }
            }
        }
        // Build all of the strands, remembering which tracker they come from.
        let mut cfg_as_vec = Vec::new();
        for (name, cfg) in &built_cfg {
            if let Some(strands) = &cfg.strands {
                for (ii, strand) in strands.iter().enumerate() {
                    cfg_as_vec.push((name.clone(), ii, *strand));
                }
            }
        }
        // Iterate through the strands by chronological order. Cannot use maps because we change types.
        cfg_as_vec.sort_by_key(|(_, _, strand)| strand.start);
        for (ii, (this_name, this_pos, this_strand)) in
            cfg_as_vec.iter().take(cfg_as_vec.len() - 1).enumerate()
        {
            // Grab the config
            if let Some(config) = self.configs[this_name].scheduler.as_ref() {
                // Grab the next strand, chronologically
                if let Some((next_name, next_pos, next_strand)) = cfg_as_vec.get(ii + 1) {
                    if config.handoff == Handoff::Greedy && this_strand.end >= next_strand.start {
                        // Modify the built configurations to change the start time of the next strand because the current one is greedy.
                        let next_config = built_cfg.get_mut(next_name).unwrap();
                        let new_start = this_strand.end + next_config.sampling;
                        next_config.strands.as_mut().unwrap()[*next_pos].start = new_start;
                        info!(
                            "{this_name} configured as {:?}, so {next_name} now starts on {new_start}",
                            config.handoff
                        );
                    } else if config.handoff == Handoff::Eager
                        && this_strand.end >= next_strand.start
                    {
                        let this_config = built_cfg.get_mut(this_name).unwrap();
                        let new_end = next_strand.start - this_config.sampling;
                        this_config.strands.as_mut().unwrap()[*this_pos].end = new_end;
                        info!(
                            "{this_name} now hands off to {next_name} on {new_end} because it's configured as {:?}",
                            config.handoff
                        );
                    }
                } else {
                    // Reached the end
                    break;
                }
            }
        }

        Ok(built_cfg)
    }

    /// Sets the schedule to that built in `build_schedule`
    pub fn build_schedule(&mut self, almanac: Arc<Almanac>) -> Result<(), NyxError> {
        self.configs = self.generate_schedule(almanac)?;

        Ok(())
    }
}
