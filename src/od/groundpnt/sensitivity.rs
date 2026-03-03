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
use crate::od::groundpnt::GroundAsset;
use crate::od::interlink::InterlinkTxSpacecraft;
use crate::od::msr::{sensitivity::ScalarSensitivity, Measurement, MeasurementType};
use crate::od::prelude::sensitivity::{ScalarSensitivityT, TrackerSensitivity};
use crate::od::{ODAlmanacSnafu, ODError};
use crate::{Spacecraft, State};
use anise::ephemerides::EphemerisPhysicsSnafu;
use anise::errors::EphemerisSnafu;
use anise::prelude::Almanac;
use indexmap::IndexSet;
use nalgebra::{DimName, OMatrix, U1};
use snafu::ResultExt;
use std::marker::PhantomData;
use std::sync::Arc;

impl TrackerSensitivity<GroundAsset, GroundAsset> for InterlinkTxSpacecraft
where
    DefaultAllocator: Allocator<<Spacecraft as State>::Size>
        + Allocator<<Spacecraft as State>::VecLength>
        + Allocator<<Spacecraft as State>::Size, <Spacecraft as State>::Size>,
{
    fn h_tilde<M: DimName>(
        &self,
        msr: &Measurement,
        msr_types: &IndexSet<MeasurementType>,
        rx: &GroundAsset,
        almanac: Arc<Almanac>,
    ) -> Result<OMatrix<f64, M, <GroundAsset as State>::Size>, ODError>
    where
        DefaultAllocator: Allocator<M> + Allocator<M, <GroundAsset as State>::Size>,
    {
        // Rebuild each row of the scalar sensitivities.
        let mut mat = OMatrix::<f64, M, <GroundAsset as State>::Size>::identity();
        for (ith_row, msr_type) in msr_types.iter().enumerate() {
            if !msr.data.contains_key(msr_type) {
                // Skip computation, this row is zero anyway.
                continue;
            }
            let scalar_h =
                <ScalarSensitivity<GroundAsset, GroundAsset, InterlinkTxSpacecraft> as ScalarSensitivityT<
                    GroundAsset,
                    GroundAsset,
                    InterlinkTxSpacecraft,
                >>::new(*msr_type, msr, rx, self, almanac.clone())?;

            mat.set_row(ith_row, &scalar_h.sensitivity_row);
        }
        Ok(mat)
    }
}

impl ScalarSensitivityT<GroundAsset, GroundAsset, InterlinkTxSpacecraft>
    for ScalarSensitivity<GroundAsset, GroundAsset, InterlinkTxSpacecraft>
{
    /// First, we ensure that the transmitter vehicle is expressed in the same frame as the
    /// ground asset. Then we compute the AER as seen from the ground asset.
    fn new(
        msr_type: MeasurementType,
        msr: &Measurement,
        rx: &GroundAsset,
        tx: &InterlinkTxSpacecraft,
        almanac: Arc<Almanac>,
    ) -> Result<Self, ODError> {
        let rx_orbit = rx.orbit();

        // Compute the SEZ DCM
        // SEZ DCM is topo to fixed
        let sez_dcm = rx_orbit
            .dcm_from_topocentric_to_body_fixed()
            .context(EphemerisPhysicsSnafu { action: "" })
            .context(EphemerisSnafu {
                action: "computing SEZ DCM for sensitivity",
            })
            .context(ODAlmanacSnafu { action: "" })?;

        let rx_sez = (sez_dcm.transpose() * rx_orbit)
            .context(EphemerisPhysicsSnafu { action: "" })
            .context(EphemerisSnafu {
                action: "transforming ground asset to SEZ",
            })
            .context(ODAlmanacSnafu { action: "" })?;

        // Convert the transmitter/PNT vehicle into the body fixed transmitter frame.
        let tx_in_rx_frame = almanac
            .transform_to(tx.traj.at(rx.epoch).unwrap().orbit, rx.frame, None)
            .context(ODAlmanacSnafu {
                action: "computing transmitter location when computing sensitivity matrix",
            })?;

        // Convert into SEZ frame
        let tx_sez = (sez_dcm.transpose() * tx_in_rx_frame)
            .context(EphemerisPhysicsSnafu { action: "" })
            .context(EphemerisSnafu {
                action: "transforming received to SEZ",
            })
            .context(ODAlmanacSnafu { action: "" })?;

        // Compute the range ρ in the SEZ frame
        let delta_r_km = tx_sez.radius_km - rx_sez.radius_km;
        let ρ_km_sez = delta_r_km.norm();
        // Compute the velocity difference - BUT note that rx_in_tx_frame is already the relative velocity of rx wrt tx!
        let delta_v_km_s = tx_in_rx_frame.velocity_km_s;

        match msr_type {
            MeasurementType::Doppler => {
                // If we have a simultaneous measurement of the range, use that, otherwise we compute the expected range.
                let ρ_km = match msr.data.get(&MeasurementType::Range) {
                    Some(range_km) => *range_km,
                    None => ρ_km_sez,
                };

                let ρ_dot_km_s = msr.data.get(&MeasurementType::Doppler).unwrap();
                let m11 = delta_r_km.x / ρ_km;
                let m12 = delta_r_km.y / ρ_km;
                let m13 = delta_r_km.z / ρ_km;
                let m21 = delta_v_km_s.x / ρ_km - ρ_dot_km_s * delta_r_km.x / ρ_km.powi(2);
                let m22 = delta_v_km_s.y / ρ_km - ρ_dot_km_s * delta_r_km.y / ρ_km.powi(2);
                let m23 = delta_v_km_s.z / ρ_km - ρ_dot_km_s * delta_r_km.z / ρ_km.powi(2);

                let sensitivity_row =
                    OMatrix::<f64, U1, <GroundAsset as State>::Size>::from_row_slice(&[
                        m21, m22, m23, m11, m12, m13,
                    ]);

                Ok(Self {
                    sensitivity_row,
                    _rx: PhantomData::<_>,
                    _tx: PhantomData::<_>,
                })
            }
            MeasurementType::Range => {
                let ρ_km = msr.data.get(&MeasurementType::Range).unwrap();
                let m11 = delta_r_km.x / ρ_km;
                let m12 = delta_r_km.y / ρ_km;
                let m13 = delta_r_km.z / ρ_km;

                let sensitivity_row =
                    OMatrix::<f64, U1, <GroundAsset as State>::Size>::from_row_slice(&[
                        m11, m12, m13, 0.0, 0.0, 0.0,
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
