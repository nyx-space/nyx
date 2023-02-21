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

use hifitime::{Epoch, TimeSeries};
use rand::rngs::StdRng;
use rand::SeedableRng;

pub use crate::dynamics::{Dynamics, NyxError};
use crate::md::trajectory::InterpState;
use crate::od::Measurement;
pub use crate::{cosmic::Cosm, State, TimeTagged};
use crate::{linalg::allocator::Allocator, od::TrackingDataSim};
use crate::{linalg::DefaultAllocator, md::ui::Traj};
use std::marker::PhantomData;
use std::sync::Arc;

pub struct TrackingArcSim<MsrIn, Msr, D>
where
    D: TrackingDataSim<MsrIn, Msr>,
    MsrIn: State,
    Msr: Measurement<State = MsrIn>,
    <Msr as Measurement>::State: InterpState,
    DefaultAllocator: Allocator<f64, <Msr::State as State>::Size>
        + Allocator<f64, <Msr::State as State>::Size, <Msr::State as State>::Size>
        + Allocator<f64, <Msr::State as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <Msr::State as State>::Size>,
{
    /// List of devices.
    pub devices: Vec<D>,
    /// Receiver trajectory
    pub trajectory: Traj<Msr::State>,
    /// Random number generator used for this tracking arc, ensures repeatability
    rng: StdRng,
    /// Greatest common denominator time series that allows this arc to meet all of the conditions.
    gcd_time_series: TimeSeries,
    _msr_in: PhantomData<MsrIn>,
    _msr: PhantomData<Msr>,
}

impl<MsrIn, Msr, D> TrackingArcSim<MsrIn, Msr, D>
where
    D: TrackingDataSim<MsrIn, Msr>,
    MsrIn: State,
    Msr: Measurement<State = MsrIn>,
    <Msr as Measurement>::State: InterpState,
    DefaultAllocator: Allocator<f64, <Msr::State as State>::Size>
        + Allocator<f64, <Msr::State as State>::Size, <Msr::State as State>::Size>
        + Allocator<f64, <Msr::State as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <Msr::State as State>::Size>,
{
    pub fn with_rng(devices: Vec<D>, trajectory: Traj<Msr::State>, rng: StdRng) -> Self {
        // TODO: Support other schedules than Continuous.

        // Compute the smallest sampling of all the devices
        let step = devices
            .iter()
            .map(|d| d.config())
            .map(|c| c.sampling)
            .min()
            .unwrap();

        let ts = TimeSeries::inclusive(trajectory.first().epoch(), trajectory.last().epoch(), step);

        Self {
            devices,
            trajectory,
            rng,
            gcd_time_series: ts,
            _msr_in: PhantomData,
            _msr: PhantomData,
        }
    }

    pub fn with_seed(devices: Vec<D>, trajectory: Traj<Msr::State>, seed: u64) -> Self {
        let rng = StdRng::seed_from_u64(seed);

        Self::with_rng(devices, trajectory, rng)
    }

    pub fn new(devices: Vec<D>, trajectory: Traj<Msr::State>) -> Self {
        let rng = StdRng::from_entropy();

        Self::with_rng(devices, trajectory, rng)
    }

    /// Generates measurements from the simulated tracking arc.
    ///
    /// Notes:
    /// Although mutable, this function may be called several times to generate different measurements.
    pub fn generate_measurements(&mut self, cosm: Arc<Cosm>) -> Result<Vec<Msr>, NyxError> {
        let start = Epoch::now().unwrap();
        let mut measurements = Vec::new();
        // Clone the time series so we don't consume it.
        let ts = self.gcd_time_series.clone();
        for epoch in ts {
            // Get the state
            let state = self.trajectory.at(epoch)?;
            for device in self.devices.iter_mut() {
                if let Some(msr) = device.measure(&state, &mut self.rng, cosm.clone()) {
                    measurements.push(msr);
                }
            }
        }

        info!(
            "Generated {} measurements in {}",
            measurements.len(),
            Epoch::now().unwrap() - start
        );

        Ok(measurements)
    }
}
