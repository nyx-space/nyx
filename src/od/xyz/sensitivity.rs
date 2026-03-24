use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::od::ODError;
use crate::{Spacecraft, State};
use anise::prelude::Almanac;
use indexmap::IndexSet;
use nalgebra::{DimName, OMatrix, U1};
use std::marker::PhantomData;
use std::sync::Arc;

use super::XyzDevice;
use crate::od::msr::measurement::Measurement;
use crate::od::msr::sensitivity::{ScalarSensitivity, ScalarSensitivityT, TrackerSensitivity};
use crate::od::msr::MeasurementType;

impl TrackerSensitivity<Spacecraft, Spacecraft> for XyzDevice
where
    DefaultAllocator: Allocator<<Spacecraft as State>::Size>
        + Allocator<<Spacecraft as State>::VecLength>
        + Allocator<<Spacecraft as State>::Size, <Spacecraft as State>::Size>,
{
    fn h_tilde<M: DimName>(
        &self,
        msr: &Measurement,
        msr_types: &IndexSet<MeasurementType>,
        rx: &Spacecraft,
        almanac: Arc<Almanac>,
    ) -> Result<OMatrix<f64, M, <Spacecraft as State>::Size>, ODError>
    where
        DefaultAllocator: Allocator<M> + Allocator<M, <Spacecraft as State>::Size>,
    {
        // Rebuild each row of the scalar sensitivities.
        let mut mat = OMatrix::<f64, M, <Spacecraft as State>::Size>::zeros();
        for (ith_row, msr_type) in msr_types.iter().enumerate() {
            if !msr.data.contains_key(msr_type) {
                // Skip computation, this row is zero anyway.
                continue;
            }
            let scalar_h =
                <ScalarSensitivity<Spacecraft, Spacecraft, XyzDevice> as ScalarSensitivityT<
                    Spacecraft,
                    Spacecraft,
                    XyzDevice,
                >>::new(*msr_type, msr, rx, self, almanac.clone())?;

            mat.set_row(ith_row, &scalar_h.sensitivity_row);
        }
        Ok(mat)
    }
}

impl ScalarSensitivityT<Spacecraft, Spacecraft, XyzDevice>
    for ScalarSensitivity<Spacecraft, Spacecraft, XyzDevice>
{
    fn new(
        msr_type: MeasurementType,
        _msr: &Measurement,
        _rx: &Spacecraft,
        _tx: &XyzDevice,
        _almanac: Arc<Almanac>,
    ) -> Result<Self, ODError> {
        let idx = match msr_type {
            MeasurementType::X => 0,
            MeasurementType::Y => 1,
            MeasurementType::Z => 2,
            _ => {
                return Err(ODError::MeasurementSimError {
                    details: format!("{msr_type:?} is not supported by XyzDevice"),
                })
            }
        };

        let mut sensitivity_row = OMatrix::<f64, U1, <Spacecraft as State>::Size>::zeros();
        sensitivity_row[(0, idx)] = 1.0;

        Ok(Self {
            sensitivity_row,
            _rx: PhantomData::<_>,
            _tx: PhantomData::<_>,
        })
    }
}
