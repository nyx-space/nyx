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

use hifitime::{Duration, Epoch, TimeSeries};
use num::integer::gcd;
use rand::SeedableRng;
use rand_pcg::Pcg64Mcg;

pub use crate::dynamics::{Dynamics, NyxError};
use crate::io::ConfigError;
use crate::md::trajectory::InterpState;
use crate::od::msr::arc::TrackingArc;
use crate::od::simulator::{Availability, Schedule};
use crate::od::{EstimateFrom, SimMeasurement};
pub use crate::{cosmic::Cosm, State, TimeTagged};
use crate::{linalg::allocator::Allocator, od::TrackingDeviceSim};
use crate::{linalg::DefaultAllocator, md::ui::Traj};
use std::collections::{HashMap, HashSet};
use std::fmt::Display;
use std::marker::PhantomData;
use std::sync::Arc;

use super::TrkConfig;

#[derive(Clone)]
pub struct TrackingArcSim<S, MsrIn, Msr, D>
where
    D: TrackingDeviceSim<MsrIn, Msr>,
    MsrIn: State,
    Msr: SimMeasurement<State = S>,
    S: EstimateFrom<MsrIn>,
    MsrIn: InterpState,
    DefaultAllocator: Allocator<f64, <Msr::State as State>::Size>
        + Allocator<f64, <MsrIn as State>::Size>
        + Allocator<f64, <MsrIn as State>::Size, <MsrIn as State>::Size>
        + Allocator<f64, <MsrIn as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize, <MsrIn as State>::Size>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, <S as State>::Size>
        + Allocator<f64, <S as State>::Size, <S as State>::Size>
        + Allocator<f64, <S as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize, <S as State>::Size>,
{
    /// Map of devices from their names.
    pub devices: HashMap<String, D>,
    /// Receiver trajectory
    pub trajectory: Traj<MsrIn>,
    /// Configuration of each device
    pub configs: HashMap<String, TrkConfig>,
    /// Set to true to allow for overlapping measurements (enabled by default),
    /// i.e. if N ground stations are available at a given epoch, generate N measurements not just one.
    pub allow_overlap: bool,
    /// Random number generator used for this tracking arc, ensures repeatability
    rng: Pcg64Mcg,
    /// Greatest common denominator time series that allows this arc to meet all of the conditions.
    time_series: TimeSeries,
    _msr_in: PhantomData<MsrIn>,
    _msr: PhantomData<Msr>,
}

impl<S, MsrIn, Msr, D> TrackingArcSim<S, MsrIn, Msr, D>
where
    D: TrackingDeviceSim<MsrIn, Msr>,
    MsrIn: State,
    Msr: SimMeasurement<State = S>,
    S: EstimateFrom<MsrIn>,
    MsrIn: InterpState,
    DefaultAllocator: Allocator<f64, <Msr::State as State>::Size>
        + Allocator<f64, <MsrIn as State>::Size>
        + Allocator<f64, <MsrIn as State>::Size, <MsrIn as State>::Size>
        + Allocator<f64, <MsrIn as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize, <MsrIn as State>::Size>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, <S as State>::Size>
        + Allocator<f64, <S as State>::Size, <S as State>::Size>
        + Allocator<f64, <S as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize, <S as State>::Size>,
{
    pub fn with_rng(
        devices: Vec<D>,
        trajectory: Traj<MsrIn>,
        configs: HashMap<String, TrkConfig>,
        rng: Pcg64Mcg,
    ) -> Result<Self, ConfigError> {
        // Check that each device has an associated configurations.
        // We don't care if there are more configurations than chosen devices.
        let mut devices_map = HashMap::new();
        let mut sampling_rates_ns = Vec::with_capacity(devices.len());
        for device in devices {
            if let Some(cfg) = configs.get(&device.name()) {
                sampling_rates_ns.push(cfg.sampling.truncated_nanoseconds());
            } else {
                return Err(ConfigError::InvalidConfig(format!(
                    "device {} has no associated configuration",
                    device.name()
                )));
            }
            devices_map.insert(device.name(), device);
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
            devices: devices_map,
            trajectory,
            configs,
            allow_overlap: false,
            rng,
            time_series,
            _msr_in: PhantomData,
            _msr: PhantomData,
        };

        info!("{me}");

        Ok(me)
    }

    pub fn with_seed(
        devices: Vec<D>,
        trajectory: Traj<MsrIn>,
        configs: HashMap<String, TrkConfig>,
        seed: u64,
    ) -> Result<Self, ConfigError> {
        let rng = Pcg64Mcg::new(seed as u128);

        Self::with_rng(devices, trajectory, configs, rng)
    }

    pub fn new(
        devices: Vec<D>,
        trajectory: Traj<MsrIn>,
        configs: HashMap<String, TrkConfig>,
    ) -> Result<Self, ConfigError> {
        let rng = Pcg64Mcg::from_entropy();

        Self::with_rng(devices, trajectory, configs, rng)
    }

    /// Disallows overlapping measurements
    pub fn disallow_overlap(&mut self) {
        self.allow_overlap = false;
    }

    /// Generates measurements from the simulated tracking arc.
    ///
    /// Notes:
    /// Although mutable, this function may be called several times to generate different measurements.
    pub fn generate_measurements(&mut self, cosm: Arc<Cosm>) -> Result<TrackingArc<Msr>, NyxError> {
        // Stores the first measurement and last measurement of a given sub-arc for each device.
        #[derive(Copy, Clone, Debug)]
        struct SchedData {
            start: Option<Epoch>,
            prev: Option<Epoch>,
            end: Option<Epoch>,
        }
        let mut sched: HashMap<String, SchedData> = HashMap::new();

        let mut start_trace_msg = HashSet::new();
        let mut end_trace_msg = HashSet::new();
        let mut sched_trace_msg = HashSet::new();

        let start = Epoch::now().unwrap();
        let mut measurements = Vec::new();
        // Clone the time series so we don't consume it.
        let ts = self.time_series.clone();
        'ts: for epoch in ts {
            'devices: for (name, device) in self.devices.iter_mut() {
                let cfg = &self.configs[name];
                // Check the start condition
                if let Availability::Epoch(start_epoch) = cfg.start {
                    if start_epoch > epoch {
                        if !start_trace_msg.contains(name) {
                            debug!(
                                "{name} tracking starts in {} (start = {start_epoch})",
                                start_epoch - epoch
                            );
                            start_trace_msg.insert(name.clone());
                        }
                        continue;
                    }
                }
                // Check the end condition
                if let Availability::Epoch(end_epoch) = cfg.end {
                    if end_epoch < epoch {
                        if !end_trace_msg.contains(name) {
                            debug!(
                                "{name} tracking ended {} ago (end = {end_epoch})",
                                epoch - end_epoch
                            );
                            end_trace_msg.insert(name.clone());
                        }
                        continue;
                    }
                }

                // Check the schedule
                if let Some(device_sched) = sched.get(name) {
                    // Check that we aren't generating too many measurements
                    if let Some(prev_epoch) = device_sched.prev {
                        if (epoch - prev_epoch) < cfg.sampling {
                            continue;
                        }
                    }
                    match cfg.schedule {
                        Schedule::Intermittent { on, off } => {
                            // Check that we haven't been on for too long
                            if let Some(start_epoch) = device_sched.start {
                                if (epoch - start_epoch) > on {
                                    if !sched_trace_msg.contains(name) {
                                        debug!(
                                            "{name} is now turned off after being on for {}",
                                            epoch - start_epoch
                                        );
                                        sched_trace_msg.insert(name.clone());
                                    }
                                    continue;
                                }
                            }
                            // Check that we haven't been off for enough time
                            if let Some(end_epoch) = device_sched.end {
                                if (epoch - end_epoch) <= off {
                                    if !sched_trace_msg.contains(name) {
                                        debug!(
                                            "{name} will be available again in {}",
                                            epoch - end_epoch
                                        );
                                        sched_trace_msg.insert(name.clone());
                                    }
                                    continue;
                                }
                            }
                        }
                        Schedule::Continuous => {
                            // No filtering, pass through
                        }
                    }
                }

                // Check the exclusion epochs
                if let Some(excl_list) = &cfg.exclusion_epochs {
                    for excl in excl_list {
                        if excl.contains(epoch) {
                            // We are in an exclusion epoch, move to next device.
                            continue 'devices;
                        }
                    }
                }

                // Check the inclusion epochs
                if let Some(incl_list) = &cfg.inclusion_epochs {
                    for incl in incl_list {
                        if !incl.contains(epoch) {
                            // Current epoch is not included in the inclusion epochs list, move to next device.
                            continue 'devices;
                        }
                    }
                }

                // Availability OK, so let's remove this device from the trace messages if needed.
                start_trace_msg.remove(name);
                sched_trace_msg.remove(name);
                end_trace_msg.remove(name);

                if let Some(msr) =
                    device.measure(epoch, &self.trajectory, Some(&mut self.rng), cosm.clone())
                {
                    measurements.push((name.clone(), msr));
                    // We have a new measurement, let's update the schedule.
                    if let Some(device_sched) = sched.get_mut(name) {
                        if device_sched.start.is_none() {
                            // Set the start time of this pass
                            device_sched.start = Some(epoch);
                            debug!("{name} is now tracking {epoch}");
                        }
                        // In any case, set the end to none and set the prev to now.
                        device_sched.prev = Some(epoch);
                        device_sched.end = None;
                    } else {
                        debug!("{name} is now tracking {epoch}");
                        // Oh, great, first measurement for this device!
                        sched.insert(
                            name.to_string(),
                            SchedData {
                                start: Some(epoch),
                                prev: Some(epoch),
                                end: None,
                            },
                        );
                    }

                    if !self.allow_overlap {
                        // No need to check for the availability of other around station for this epoch
                        // so let's move to the next epoch in the time series.
                        continue 'ts;
                    }
                } else if let Some(device_sched) = sched.get_mut(name) {
                    // No measurement was generated, so let's update the schedule if there is one.
                    if device_sched.end.is_none() {
                        device_sched.start = None;
                        device_sched.end = Some(epoch);
                    }
                }
            }
        }

        info!(
            "Generated {} measurements in {}",
            measurements.len(),
            Epoch::now().unwrap() - start
        );

        let mut devices = Vec::new();
        for device in self.devices.values() {
            let repr = device.to_config().unwrap();
            devices.push(repr);
        }

        // Build the tracking arc.
        let trk = TrackingArc {
            device_cfg: serde_yaml::to_string(&devices).unwrap(),
            measurements,
        };

        Ok(trk)
    }
}

impl<S, MsrIn, Msr, D> Display for TrackingArcSim<S, MsrIn, Msr, D>
where
    D: TrackingDeviceSim<MsrIn, Msr>,
    MsrIn: State,
    Msr: SimMeasurement<State = S>,
    S: EstimateFrom<MsrIn>,
    MsrIn: InterpState,
    DefaultAllocator: Allocator<f64, <Msr::State as State>::Size>
        + Allocator<f64, <MsrIn as State>::Size>
        + Allocator<f64, <MsrIn as State>::Size, <MsrIn as State>::Size>
        + Allocator<f64, <MsrIn as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize, <MsrIn as State>::Size>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, <S as State>::Size>
        + Allocator<f64, <S as State>::Size, <S as State>::Size>
        + Allocator<f64, <S as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize, <S as State>::Size>,
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
