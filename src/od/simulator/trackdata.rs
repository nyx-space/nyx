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
use indexmap::IndexSet;
use nalgebra::{DimName, OMatrix};
use rand_pcg::Pcg64Mcg;

use crate::io::ConfigRepr;
use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::md::prelude::{Frame, Traj};
use crate::md::trajectory::Interpolatable;
use crate::od::msr::measurement::Measurement as NewMeasurement;
use crate::od::msr::MeasurementType;
use crate::od::ODError;
use crate::Orbit;

/// Tracking device simulator.
pub trait TrackingDevice<MsrIn>: ConfigRepr
where
    MsrIn: Interpolatable,
    DefaultAllocator:
        Allocator<MsrIn::Size> + Allocator<MsrIn::Size, MsrIn::Size> + Allocator<MsrIn::VecLength>,
{
    /// Returns the name of this tracking data simulator
    fn name(&self) -> String;

    /// Returns the _enabled_ measurement types for thie device.
    fn measurement_types(&self) -> &IndexSet<MeasurementType>;

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
    ) -> Result<Option<NewMeasurement>, ODError>;

    /// Returns the device location at the given epoch and in the given frame.
    fn location(&self, epoch: Epoch, frame: Frame, almanac: Arc<Almanac>) -> AlmanacResult<Orbit>;

    // Perform an instantaneous measurement (without integration times, i.e. one-way). Returns None if the object is not visible, else returns the measurement.
    fn measure_instantaneous(
        &mut self,
        rx: MsrIn,
        rng: Option<&mut Pcg64Mcg>,
        almanac: Arc<Almanac>,
    ) -> Result<Option<NewMeasurement>, ODError>;

    // Return the noise statistics of this tracking device for the provided measurement type at the requested epoch.
    fn measurement_covar(&self, msr_type: MeasurementType, epoch: Epoch) -> Result<f64, ODError>;

    fn measurement_covar_matrix<M: DimName>(
        &self,
        msr_types: &IndexSet<MeasurementType>,
        epoch: Epoch,
    ) -> Result<OMatrix<f64, M, M>, ODError>
    where
        DefaultAllocator: Allocator<M, M>,
    {
        // Rebuild the R matrix of the measurement noise.
        let mut r_mat = OMatrix::<f64, M, M>::zeros();

        for (i, msr_type) in msr_types.iter().enumerate() {
            if self.measurement_types().contains(msr_type) {
                r_mat[(i, i)] = self.measurement_covar(*msr_type, epoch)?;
            }
        }

        Ok(r_mat)
    }
}
