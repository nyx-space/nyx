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

use hifitime::Epoch;
use rand_pcg::Pcg64Mcg;

use crate::linalg::DefaultAllocator;
use crate::md::trajectory::InterpState;
use crate::md::ui::Traj;
use crate::od::Measurement;
use crate::State;
use crate::{io::Configurable, linalg::allocator::Allocator};

use super::Cosm;

/// Tracking device simulator.
pub trait TrackingDeviceSim<MsrIn, Msr>: Configurable
where
    MsrIn: InterpState,
    Msr: Measurement,
    DefaultAllocator: Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, MsrIn::Size>
        + Allocator<f64, MsrIn::Size, MsrIn::Size>
        + Allocator<f64, MsrIn::VecLength>,
{
    /// Returns the name of this tracking data simulator
    fn name(&self) -> String;

    /// Observes the input trajectory at the provided epoch, and returns a measurement from itself to the input state, and returns None of the object is not visible.
    /// This trait function takes in a trajectory and epoch so it can properly simulate integration times for the measurements.
    /// If the random number generator is provided, it shall be used to add noise to the measurement.
    ///
    /// # Choice of the random number generator
    /// The Pcg64Mcg is chosen because it is fast, space efficient, and has a good statistical distribution.
    fn measure(
        &mut self,
        epoch: Epoch,
        traj: &Traj<MsrIn>,
        rng: Option<&mut Pcg64Mcg>,
        cosm: Arc<Cosm>,
    ) -> Option<(Msr, MsrIn)>;
}
