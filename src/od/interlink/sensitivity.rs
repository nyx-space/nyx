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

use crate::linalg::allocator::Allocator;
use crate::linalg::DefaultAllocator;
use crate::od::interlink::InterlinkTxSpacecraft;
use crate::od::msr::{Measurement, MeasurementType};
use crate::od::prelude::sensitivity::{ScalarSensitivityT, TrackerSensitivity};
use crate::od::{ODAlmanacSnafu, ODError, TrackingDevice};
use crate::{Spacecraft, State};
use anise::prelude::Almanac;
use indexmap::IndexSet;
use nalgebra::{DimName, OMatrix, U1};
use snafu::ResultExt;
use std::marker::PhantomData;
use std::sync::Arc;

struct ScalarSensitivity<SolveState: State, Rx, Tx>
where
    DefaultAllocator: Allocator<SolveState::Size>
        + Allocator<SolveState::VecLength>
        + Allocator<SolveState::Size, SolveState::Size>
        + Allocator<U1, SolveState::Size>,
{
    sensitivity_row: OMatrix<f64, U1, SolveState::Size>,
    _rx: PhantomData<Rx>,
    _tx: PhantomData<Tx>,
}

impl TrackerSensitivity<Spacecraft, Spacecraft> for InterlinkTxSpacecraft
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
        let mut mat = OMatrix::<f64, M, <Spacecraft as State>::Size>::identity();
        for (ith_row, msr_type) in msr_types.iter().enumerate() {
            if !msr.data.contains_key(msr_type) {
                // Skip computation, this row is zero anyway.
                continue;
            }
            let scalar_h =
                <ScalarSensitivity<Spacecraft, Spacecraft, InterlinkTxSpacecraft> as ScalarSensitivityT<
                    Spacecraft,
                    Spacecraft,
                    InterlinkTxSpacecraft,
                >>::new(*msr_type, msr, rx, self, almanac.clone())?;

            mat.set_row(ith_row, &scalar_h.sensitivity_row);
        }
        Ok(mat)
    }
}

impl ScalarSensitivityT<Spacecraft, Spacecraft, InterlinkTxSpacecraft>
    for ScalarSensitivity<Spacecraft, Spacecraft, InterlinkTxSpacecraft>
{
    fn new(
        msr_type: MeasurementType,
        msr: &Measurement,
        rx: &Spacecraft,
        tx: &InterlinkTxSpacecraft,
        almanac: Arc<Almanac>,
    ) -> Result<Self, ODError> {
        let receiver = rx.orbit;

        // Compute the device location in the receiver frame because we compute the sensitivity in that frame.
        // This frame is required because the scalar measurements are frame independent, but the sensitivity
        // must be in the estimation frame.
        let transmitter = tx
            .location(rx.orbit.epoch, rx.orbit.frame, almanac.clone())
            .context(ODAlmanacSnafu {
                action: "computing transmitter location when computing sensitivity matrix",
            })?;

        let delta_r = receiver.radius_km - transmitter.radius_km;
        let delta_v = receiver.velocity_km_s - transmitter.velocity_km_s;

        match msr_type {
            MeasurementType::Doppler => {
                // If we have a simultaneous measurement of the range, use that, otherwise we compute the expected range.
                let ρ_km = msr.data.get(&MeasurementType::Range).unwrap();
                let ρ_dot_km_s = msr.data.get(&MeasurementType::Doppler).unwrap();
                let m11 = delta_r.x / ρ_km;
                let m12 = delta_r.y / ρ_km;
                let m13 = delta_r.z / ρ_km;
                let m21 = delta_v.x / ρ_km - ρ_dot_km_s * delta_r.x / ρ_km.powi(2);
                let m22 = delta_v.y / ρ_km - ρ_dot_km_s * delta_r.y / ρ_km.powi(2);
                let m23 = delta_v.z / ρ_km - ρ_dot_km_s * delta_r.z / ρ_km.powi(2);

                let sensitivity_row =
                    OMatrix::<f64, U1, <Spacecraft as State>::Size>::from_row_slice(&[
                        m21, m22, m23, m11, m12, m13, 0.0, 0.0, 0.0,
                    ]);

                Ok(Self {
                    sensitivity_row,
                    _rx: PhantomData::<_>,
                    _tx: PhantomData::<_>,
                })
            }
            MeasurementType::Range => {
                let ρ_km = msr.data.get(&MeasurementType::Range).unwrap();
                let m11 = delta_r.x / ρ_km;
                let m12 = delta_r.y / ρ_km;
                let m13 = delta_r.z / ρ_km;

                let sensitivity_row =
                    OMatrix::<f64, U1, <Spacecraft as State>::Size>::from_row_slice(&[
                        m11, m12, m13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    ]);

                Ok(Self {
                    sensitivity_row,
                    _rx: PhantomData::<_>,
                    _tx: PhantomData::<_>,
                })
            }
            MeasurementType::Azimuth
            | MeasurementType::Elevation
            | MeasurementType::ReceiveFrequency
            | MeasurementType::TransmitFrequency => Err(ODError::MeasurementSimError {
                details: format!("{msr_type:?} is not supported for interlink"),
            }),
        }
    }
}
