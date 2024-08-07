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

use std::sync::Arc;

use anise::almanac::Almanac;
use anise::errors::AlmanacResult;
use hifitime::Epoch;
use rand_pcg::Pcg64Mcg;

use crate::io::ConfigRepr;
use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, OMatrix};
use crate::md::prelude::{Frame, Traj};
use crate::md::trajectory::Interpolatable;
use crate::od::{Measurement, ODError};
use crate::Orbit;

/// Tracking device simulator.
pub trait TrackingDeviceSim<MsrIn, Msr>: ConfigRepr
where
    MsrIn: Interpolatable,
    Msr: Measurement,
    DefaultAllocator: Allocator<Msr::MeasurementSize>
        + Allocator<MsrIn::Size>
        + Allocator<MsrIn::Size, MsrIn::Size>
        + Allocator<Msr::MeasurementSize, Msr::MeasurementSize>
        + Allocator<MsrIn::VecLength>,
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
    /// + A specific measurement is requested but the noise on that measurement type is not configured.
    ///
    fn measure(
        &mut self,
        epoch: Epoch,
        traj: &Traj<MsrIn>,
        rng: Option<&mut Pcg64Mcg>,
        almanac: Arc<Almanac>,
    ) -> Result<Option<Msr>, ODError>;

    /// Returns the device location at the given epoch and in the given frame.
    fn location(&self, epoch: Epoch, frame: Frame, almanac: Arc<Almanac>) -> AlmanacResult<Orbit>;

    // Perform an instantaneous measurement (without integration times, i.e. one-way). Returns None if the object is not visible, else returns the measurement.
    fn measure_instantaneous(
        &mut self,
        rx: MsrIn,
        rng: Option<&mut Pcg64Mcg>,
        almanac: Arc<Almanac>,
    ) -> Result<Option<Msr>, ODError>;

    // Return the noise statistics of this tracking device at the requested epoch.
    fn measurement_noise(
        &mut self,
        epoch: Epoch,
    ) -> Result<OMatrix<f64, Msr::MeasurementSize, Msr::MeasurementSize>, ODError>;
}
