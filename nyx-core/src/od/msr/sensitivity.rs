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

use crate::io::ConfigError;
use crate::linalg::DefaultAllocator;
use crate::linalg::allocator::Allocator;
use crate::md::prelude::Interpolatable;
use crate::od::{GroundStation, ODAlmanacSnafu, ODError, TrackingDevice};
use crate::{Spacecraft, State};
use anise::constants::SPEED_OF_LIGHT_KM_S;
use anise::errors::OrientationSnafu;
use anise::prelude::Almanac;
use indexmap::IndexSet;
use nalgebra::{DimName, OMatrix, U1};
use snafu::ResultExt;
use std::marker::PhantomData;

use super::measurement::Measurement;
use super::{MeasurementType, TrackingDataArc};

pub trait ScalarSensitivityT<SolveState: State, Rx, Tx>
where
    Self: Sized,
    DefaultAllocator: Allocator<SolveState::Size>
        + Allocator<SolveState::VecLength>
        + Allocator<SolveState::Size, SolveState::Size>,
{
    fn new(
        msr_type: MeasurementType,
        msr: &Measurement,
        rx: &Rx,
        tx: &Tx,
        almanac: &Almanac,
    ) -> Result<Self, ODError>;
}

/// Trait required to build a triplet of a solve-for state, a receiver, and a transmitter.
pub trait TrackerSensitivity<SolveState: Interpolatable, Rx>: TrackingDevice<SolveState>
where
    Self: Sized,
    DefaultAllocator: Allocator<SolveState::Size>
        + Allocator<SolveState::VecLength>
        + Allocator<SolveState::Size, SolveState::Size>,
{
    /// Returns the sensitivity matrix of size MxS where M is the number of simultaneous measurements
    /// and S is the size of the state being solved for.
    fn h_tilde<M: DimName>(
        &self,
        msr: &Measurement,
        msr_types: &IndexSet<MeasurementType>, // Consider switching to array
        rx: &Rx,
        almanac: &Almanac,
    ) -> Result<OMatrix<f64, M, SolveState::Size>, ODError>
    where
        DefaultAllocator: Allocator<M> + Allocator<M, SolveState::Size>;

    /// Returns whether this tracker is expected to be compatible with the tracking data arc
    fn is_compatible(&self, _tracker: &str, _arc: &TrackingDataArc) -> Result<(), ODError> {
        Ok(())
    }
}

pub struct ScalarSensitivity<SolveState: State, Rx, Tx>
where
    DefaultAllocator: Allocator<SolveState::Size>
        + Allocator<SolveState::VecLength>
        + Allocator<SolveState::Size, SolveState::Size>
        + Allocator<U1, SolveState::Size>,
{
    pub sensitivity_row: OMatrix<f64, U1, SolveState::Size>,
    pub _rx: PhantomData<Rx>,
    pub _tx: PhantomData<Tx>,
}

impl TrackerSensitivity<Spacecraft, Spacecraft> for GroundStation
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
        almanac: &Almanac,
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
                <ScalarSensitivity<Spacecraft, Spacecraft, GroundStation> as ScalarSensitivityT<
                    Spacecraft,
                    Spacecraft,
                    GroundStation,
                >>::new(*msr_type, msr, rx, self, almanac)?;

            mat.set_row(ith_row, &scalar_h.sensitivity_row);
        }
        Ok(mat)
    }

    fn is_compatible(&self, tracker: &str, arc: &TrackingDataArc) -> Result<(), ODError> {
        // Ensure that the arc doppler config matches this ground station doppler config.
        for msr in &arc.measurements {
            if msr.doppler_config != self.doppler_config {
                return Err(ODError::ODConfigError {
                    source: ConfigError::InvalidConfig {
                        msg: format!(
                            "Tracker `{tracker}` Doppler config does not match measurement config @ {}\nGround Station:\n{:?}\nMeasurement:\n{:?}",
                            msr.epoch, self.doppler_config, msr.doppler_config
                        ),
                    },
                });
            }
        }
        Ok(())
    }
}

impl ScalarSensitivityT<Spacecraft, Spacecraft, GroundStation>
    for ScalarSensitivity<Spacecraft, Spacecraft, GroundStation>
{
    fn new(
        msr_type: MeasurementType,
        _msr: &Measurement,
        rx: &Spacecraft,
        tx: &GroundStation,
        almanac: &Almanac,
    ) -> Result<Self, ODError> {
        let receiver = rx.orbit;

        // Compute the device location in the receiver frame because we compute the sensitivity in that frame.
        // This frame is required because the scalar measurements are frame independent, but the sensitivity
        // must be in the estimation frame.
        let transmitter = tx
            .location(receiver.epoch, receiver.frame, almanac)
            .context(ODAlmanacSnafu {
                action: "computing transmitter location when computing sensitivity matrix",
            })?;

        // Relative geometry in estimation frame
        let delta_r = receiver.radius_km - transmitter.radius_km;
        let delta_v = receiver.velocity_km_s - transmitter.velocity_km_s;

        let rho_km = delta_r.norm();
        // let rho_km = *msr.data.get(&MeasurementType::Range).unwrap();
        if rho_km < 1e-6 {
            return Err(ODError::MeasurementSimError {
                details: "Zero separation between ground station and spacecraft".to_string(),
            });
        }

        // Line-of-sight unit vector pointing from station to spacecraft: d(rho)/d(r)
        let u_los = delta_r / rho_km;

        let sensitivity_row = match msr_type {
            MeasurementType::Doppler => {
                // Nominal line-of-sight range-rate from trajectory geometry
                let rho_dot_km_s = u_los.dot(&delta_v);
                // let rho_dot_km_s = msr.data.get(&MeasurementType::Doppler).unwrap();
                let m11 = delta_r.x / rho_km;
                let m12 = delta_r.y / rho_km;
                let m13 = delta_r.z / rho_km;
                let m21 = delta_v.x / rho_km - rho_dot_km_s * delta_r.x / rho_km.powi(2);
                let m22 = delta_v.y / rho_km - rho_dot_km_s * delta_r.y / rho_km.powi(2);
                let m23 = delta_v.z / rho_km - rho_dot_km_s * delta_r.z / rho_km.powi(2);

                OMatrix::<f64, U1, <Spacecraft as State>::Size>::from_row_slice(&[
                    m21, m22, m23, m11, m12, m13, 0.0, 0.0, 0.0,
                ])
            }
            MeasurementType::Range => {
                // Velocity sensitivity due to retarded bounce epoch: d(rho)/d(v) = -tau * u_los
                // This is required because we're computing the sensitivity at the reception epoch
                // and not at the bounce epoch.
                let tau_s = if tx.light_time_correction {
                    rho_km / SPEED_OF_LIGHT_KM_S
                } else {
                    0.0
                };

                OMatrix::<f64, U1, <Spacecraft as State>::Size>::from_row_slice(&[
                    u_los.x,
                    u_los.y,
                    u_los.z,
                    -tau_s * u_los.x,
                    -tau_s * u_los.y,
                    -tau_s * u_los.z,
                    0.0,
                    0.0,
                    0.0,
                ])
            }

            MeasurementType::Azimuth | MeasurementType::Elevation => {
                // Transform line-of-sight vector into station topocentric SEZ frame
                let rx_in_tx = almanac
                    .transform_to(receiver, tx.location.frame.into(), None)
                    .context(ODAlmanacSnafu {
                        action: "transforming receiver to station frame for topocentric angles",
                    })?;

                let d_sez = rx_in_tx.radius_km; // Topocentric vector
                let s = d_sez.x;
                let e = d_sez.y;
                let z = d_sez.z - tx.location.height_km;
                let rho2_horiz = s * s + e * e;
                let rho2_total = rho2_horiz + z * z;

                if rho2_horiz < 1e-12 {
                    return Err(ODError::MeasurementSimError {
                        details: "Singularity at zenith for topocentric angle sensitivities".into(),
                    });
                }

                // Partials in local SEZ frame
                let d_sez_partials = if msr_type == MeasurementType::Azimuth {
                    // Azimuth beta = atan2(E, -S) => d(beta)/d(S) = E / (S^2 + E^2), d(beta)/d(E) = -S / (S^2 + E^2)
                    nalgebra::Vector3::new(e / rho2_horiz, -s / rho2_horiz, 0.0)
                } else {
                    // Elevation el = asin(Z / rho)
                    let sqrt_horiz = rho2_horiz.sqrt();
                    nalgebra::Vector3::new(
                        -(s * z) / (rho2_total * sqrt_horiz),
                        -(e * z) / (rho2_total * sqrt_horiz),
                        sqrt_horiz / rho2_total,
                    )
                };

                // Rotation matrix from SEZ back to estimation frame: R_sez_to_inertial
                let rot_mat = almanac
                    .rotate(tx.location.frame.into(), receiver.frame, receiver.epoch)
                    .context(OrientationSnafu {
                        action: "computing SEZ-to-inertial rotation",
                    })
                    .context(ODAlmanacSnafu {
                        action: "computing SEZ-to-inertial rotation for angles",
                    })?;

                let d_inertial = rot_mat * d_sez_partials;

                OMatrix::<f64, U1, <Spacecraft as State>::Size>::from_row_slice(&[
                    d_inertial.x,
                    d_inertial.y,
                    d_inertial.z,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                ])
            }

            MeasurementType::ReceiveFrequency
            | MeasurementType::TransmitFrequency
            | MeasurementType::TransmitFrequencyRate => {
                return Err(ODError::MeasurementSimError {
                    details: format!("{msr_type:?} is only supported in CCSDS TDM parsing"),
                });
            }

            MeasurementType::X | MeasurementType::Y | MeasurementType::Z => {
                return Err(ODError::MeasurementSimError {
                    details: format!("{msr_type:?} is not supported for ground stations"),
                });
            }
        };

        Ok(Self {
            sensitivity_row,
            _rx: PhantomData,
            _tx: PhantomData,
        })
    }
}
