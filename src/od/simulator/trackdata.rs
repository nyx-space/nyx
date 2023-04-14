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
use crate::md::trajectory::Interpolatable;
use crate::md::ui::{Frame, Traj};
use crate::od::Measurement;
use crate::{io::Configurable, linalg::allocator::Allocator};
use crate::{NyxError, Orbit};

use super::Cosm;

/// Tracking device simulator.
pub trait TrackingDeviceSim<MsrIn, Msr>: Configurable
where
    MsrIn: Interpolatable,
    Msr: Measurement,
    DefaultAllocator: Allocator<f64, Msr::MeasurementSize>
        + Allocator<f64, MsrIn::Size>
        + Allocator<f64, MsrIn::Size, MsrIn::Size>
        + Allocator<f64, MsrIn::VecLength>,
{
    /// Returns the name of this tracking data simulator
    fn name(&self) -> String;

    /// Performs a measurement of the input trajectory at the provided epoch (with integration times if relevant), and returns a measurement from itself to the input state. Returns None of the object is not visible.
    /// This trait function takes in a trajectory and epoch so it can properly simulate integration times for the measurements.
    /// If the random number generator is provided, it shall be used to add noise to the measurement.
    ///
    /// # Choice of the random number generator
    /// The Pcg64Mcg is chosen because it is fast, space efficient, and has a good statistical distribution.
    ///
    /// # Errors
    ///     + A specific measurement is requested but the noise on that measurement type is not configured.
    fn measure(
        &mut self,
        epoch: Epoch,
        traj: &Traj<MsrIn>,
        rng: Option<&mut Pcg64Mcg>,
        cosm: Arc<Cosm>,
    ) -> Result<Option<Msr>, NyxError>;

    /// Returns the device location at the given epoch and in the given frame.
    fn location(&self, epoch: Epoch, frame: Frame, cosm: &Cosm) -> Orbit;

    // Perform an instantaneous measurement (without integration times, i.e. one-way). Returns None if the object is not visible, else returns the measurement.
    fn measure_instantaneous(
        &mut self,
        rx: MsrIn,
        rng: Option<&mut Pcg64Mcg>,
        cosm: Arc<Cosm>,
    ) -> Result<Option<Msr>, NyxError>;
}
