/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

use super::solution::kalman::KalmanVariant;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName};
use crate::md::trajectory::Interpolatable;
pub use crate::od::snc::*;
pub use crate::od::*;
use crate::propagators::Propagator;
use anise::prelude::Almanac;
use hifitime::Unit;
use msr::sensitivity::TrackerSensitivity;
use std::collections::BTreeMap;
use std::marker::PhantomData;
use std::ops::Add;

use super::{KalmanODProcess, ResidRejectCrit};

impl<
        D: Dynamics,
        MsrSize: DimName,
        Accel: DimName,
        Trk: TrackerSensitivity<D::StateType, D::StateType>,
    > KalmanODProcess<D, MsrSize, Accel, Trk>
where
    D::StateType:
        Interpolatable + Add<OVector<f64, <D::StateType as State>::Size>, Output = D::StateType>,
    <DefaultAllocator as Allocator<<D::StateType as State>::VecLength>>::Buffer<f64>: Send,
    <DefaultAllocator as Allocator<<D::StateType as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>>::Buffer<f64>: Copy,
    DefaultAllocator: Allocator<<D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <D::StateType as State>::Size>
        + Allocator<<D::StateType as State>::Size, MsrSize>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<D::StateType as State>::Size, <D::StateType as State>::Size>
        + Allocator<Accel>
        + Allocator<Accel, Accel>
        + Allocator<<D::StateType as State>::Size, Accel>
        + Allocator<Accel, <D::StateType as State>::Size>
        + Allocator<nalgebra::Const<1>, MsrSize>,
{
    /// Initialize a new Kalman sequential filter for the orbit determination process,
    /// setting the max step of the STM to one minute, and the measurement epoch precision to 1 microsecond.
    pub fn new(
        prop: Propagator<D>,
        kf_variant: KalmanVariant,
        resid_crit: Option<ResidRejectCrit>,
        devices: BTreeMap<String, Trk>,
        almanac: Arc<Almanac>,
    ) -> Self {
        Self {
            prop,
            kf_variant,
            devices,
            resid_crit,
            process_noise: vec![],
            max_step: Unit::Minute * 1,
            epoch_precision: Unit::Microsecond * 1,
            almanac,
            _msr_size: PhantomData::<MsrSize>,
        }
    }

    /// Set (or replaces) the existing process noise configuration.
    pub fn from_process_noise(
        prop: Propagator<D>,
        kf_variant: KalmanVariant,
        devices: BTreeMap<String, Trk>,
        resid_crit: Option<ResidRejectCrit>,
        process_noise: ProcessNoise<Accel>,
        almanac: Arc<Almanac>,
    ) -> Self {
        let mut me = Self::new(
            prop,
            kf_variant,
            resid_crit,
            devices,
            almanac
        );
        me.process_noise.push(process_noise);
        me
    }

    /// Set (or replaces) the existing process noise configuration.
    pub fn with_process_noise(mut self, process_noise: ProcessNoise<Accel>) -> Self {
        self.process_noise.clear();
        self.process_noise.push(process_noise);
        self
    }

    /// Pushes the provided process noise to the list the existing process noise configurations.
    pub fn and_with_process_noise(mut self, process_noise: ProcessNoise<Accel>) -> Self {
        self.process_noise.push(process_noise);
        self
    }
}
