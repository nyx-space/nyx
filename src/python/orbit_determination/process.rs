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

use hifitime::Duration;
use nalgebra::Matrix2;
use pyo3::prelude::*;

use crate::{
    io::tracking_data::DynamicTrackingArc,
    md::ui::{Cosm, Propagator, SpacecraftDynamics},
    od::{
        kalman::KF,
        process::{EkfTrigger, ODProcess},
    },
    NyxError, Spacecraft,
};

use super::{estimate::OrbitEstimate, GroundStation};

/// Propagates the provided spacecraft with the provided dynamics until the provided stopping condition (duration, epoch, or event [and optionally the count]).
#[pyfunction]
#[pyo3(
    text_signature = "(spacecraft, dynamics, duration=None, epoch=None, event=None, event_count=None, min_step=None, max_step=None, fixed_step=None, tolerance=None)"
)]
pub(crate) fn process_tracking_arc(
    dynamics: SpacecraftDynamics,
    spacecraft: Spacecraft, // TODO: Consolidate spacecraft and estimate better
    initial_estimate: OrbitEstimate,
    measurement_noise: Vec<f64>, // TODO: Switch to a serializable matrix with Numpy?
    arc: &DynamicTrackingArc,
    ekf_num_meas: usize,
    ekf_disable_time: Duration,
) -> Result<Vec<OrbitEstimate>, NyxError> {
    // TODO: Return a navigation trajectory or use a class that mimics the better ODProcess
    let msr_noise = Matrix2::from_iterator(measurement_noise);

    let init_sc = spacecraft.with_orbit(initial_estimate.0.nominal_state.with_stm());

    // Build KF without SNC
    let kf = KF::no_snc(initial_estimate.0, msr_noise);

    let prop = Propagator::default(dynamics);
    let prop_est = prop.with(init_sc);

    let mut odp = ODProcess::ekf(
        prop_est,
        kf,
        EkfTrigger::new(ekf_num_meas, ekf_disable_time),
        Cosm::de438(),
    );

    let concrete_arc = arc.to_tracking_arc()?;

    odp.process_tracking_arc::<GroundStation>(&concrete_arc)
        .unwrap();

    // Now build a vector of orbit estimates.
    let mut rslt = Vec::with_capacity(odp.estimates.len());
    for est in odp.estimates {
        rslt.push(OrbitEstimate(est));
    }

    Ok(rslt)
}
