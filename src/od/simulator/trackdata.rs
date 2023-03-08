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

use std::sync::Arc;

use rand::Rng;

use crate::linalg::DefaultAllocator;
use crate::od::SimMeasurement;
use crate::State;
use crate::{io::Configurable, linalg::allocator::Allocator};

use super::Cosm;

/// Tracking device simulator.
pub trait TrackingDeviceSim<MsrIn, Msr>: Configurable
where
    MsrIn: State,
    Msr: SimMeasurement,
    DefaultAllocator: Allocator<f64, <Msr::State as State>::Size>
        + Allocator<f64, <Msr::State as State>::Size, <Msr::State as State>::Size>
        + Allocator<f64, <Msr::State as State>::VecLength>
        + Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, Msr::MeasurementSize, <Msr::State as State>::Size>
        + Allocator<f64, MsrIn::Size>
        + Allocator<f64, MsrIn::Size, MsrIn::Size>
        + Allocator<f64, MsrIn::VecLength>,
{
    /// Returns the name of this tracking data simulator
    fn name(&self) -> String;

    /// Observes the input state and returns a measurement from itself to the input state, and returns None of the object is not visible.
    fn measure<R: Rng>(&mut self, input: &MsrIn, rng: &mut R, cosm: Arc<Cosm>) -> Option<Msr>;
}
