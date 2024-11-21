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
use crate::linalg::{DefaultAllocator, OVector};
use crate::od::{GroundStation, TrackingDeviceSim};
use crate::{Spacecraft, State};
use anise::prelude::Almanac;
use std::marker::PhantomData;
use std::sync::Arc;

use super::measurement::Measurement;
use super::MeasurementType;

trait Sensitivity<SolveState: State, Rx, Tx>
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
        almanac: Arc<Almanac>,
    ) -> Result<Self, String>;
}

struct ScalarSensitivity<SolveState: State, Rx, Tx>
where
    DefaultAllocator: Allocator<SolveState::Size>
        + Allocator<SolveState::VecLength>
        + Allocator<SolveState::Size, SolveState::Size>,
{
    sensitivity_row: OVector<f64, SolveState::Size>,
    msr_type: MeasurementType,
    _rx: PhantomData<Rx>,
    _tx: PhantomData<Tx>,
}

impl Sensitivity<Spacecraft, Spacecraft, GroundStation>
    for ScalarSensitivity<Spacecraft, Spacecraft, GroundStation>
{
    fn new(
        msr_type: MeasurementType,
        msr: &Measurement,
        rx: &Spacecraft,
        tx: &GroundStation,
        almanac: Arc<Almanac>,
    ) -> Result<Self, String> {
        let receiver = rx.orbit;
        // Compute the device location
        let transmitter = tx
            .location(rx.orbit.epoch, rx.orbit.frame, almanac.clone())
            .unwrap();

        let delta_r = receiver.radius_km - transmitter.radius_km;
        let delta_v = receiver.velocity_km_s - transmitter.velocity_km_s;

        match msr_type {
            MeasurementType::Doppler => {
                // If we have a simultaneous measurement of the range, use that, otherwise we compute the expected range.
                let ρ_km = match msr.data.get(&MeasurementType::Range) {
                    Some(range_km) => *range_km,
                    None => {
                        tx.azimuth_elevation_of(receiver, None, &almanac)
                            .unwrap()
                            .range_km
                    }
                };

                let ρ_dot_km_s = msr.data.get(&MeasurementType::Doppler).unwrap();
                let m11 = delta_r.x / ρ_km;
                let m12 = delta_r.y / ρ_km;
                let m13 = delta_r.z / ρ_km;
                let m21 = delta_v.x / ρ_km - ρ_dot_km_s * delta_r.x / ρ_km.powi(2);
                let m22 = delta_v.y / ρ_km - ρ_dot_km_s * delta_r.y / ρ_km.powi(2);
                let m23 = delta_v.z / ρ_km - ρ_dot_km_s * delta_r.z / ρ_km.powi(2);

                let sensitivity_row =
                    OVector::<f64, <Spacecraft as State>::Size>::from_row_slice(&[
                        m21, m22, m23, m11, m12, m13, 0.0, 0.0, 0.0,
                    ]);

                Ok(Self {
                    sensitivity_row,
                    msr_type,
                    _rx: PhantomData::<_>,
                    _tx: PhantomData::<_>,
                })
            }
            MeasurementType::Range => {
                let ρ_km = msr.data.get(&MeasurementType::Doppler).unwrap();
                let m11 = delta_r.x / ρ_km;
                let m12 = delta_r.y / ρ_km;
                let m13 = delta_r.z / ρ_km;

                let sensitivity_row =
                    OVector::<f64, <Spacecraft as State>::Size>::from_row_slice(&[
                        m11, m12, m13, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    ]);

                Ok(Self {
                    sensitivity_row,
                    msr_type,
                    _rx: PhantomData::<_>,
                    _tx: PhantomData::<_>,
                })
            }
        }
    }
}
