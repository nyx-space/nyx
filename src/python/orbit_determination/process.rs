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

use hifitime::{Duration, Epoch};
use nalgebra::Matrix2;
use pyo3::prelude::*;

use crate::{
    io::tracking_data::DynamicTrackingArc,
    io::ExportCfg,
    md::prelude::{Cosm, Propagator, SpacecraftDynamics},
    od::{
        filter::kalman::KF,
        process::{EkfTrigger, FltResid, IterationConf, ODProcess},
    },
    NyxError, Spacecraft,
};

use super::{estimate::OrbitEstimate, GroundStation};

/// Runs an orbit determination process given the dynamics, the initial spacecraft object and its orbit estimate, the measurement noise for the filter, the tracking arc data.
/// You must also provide an export path and optionally and export configuration to export the results to a Parquet file.
#[pyfunction]
#[pyo3(
    text_signature = "(dynamics, spacecraft, initial_estimate, measurement_noise, arc, export_path, export_cfg, ekf_num_meas=None, ekf_disable_time=None, resid_crit=None, predict_until=None, predict_for=None, predict_step=None, fixed_step=False)"
)]
pub(crate) fn process_tracking_arc(
    dynamics: SpacecraftDynamics,
    spacecraft: Spacecraft,
    initial_estimate: OrbitEstimate,
    measurement_noise: Vec<f64>,
    arc: &DynamicTrackingArc,
    export_path: String,
    export_cfg: Option<ExportCfg>,
    ekf_num_meas: Option<usize>,
    ekf_disable_time: Option<Duration>,
    resid_crit: Option<FltResid>,
    predict_until: Option<Epoch>,
    predict_for: Option<Duration>,
    predict_step: Option<Duration>,
    fixed_step: Option<bool>,
    iter_conf: Option<IterationConf>,
) -> Result<String, NyxError> {
    let msr_noise = Matrix2::from_iterator(measurement_noise);

    let init_sc = spacecraft.with_orbit(initial_estimate.0.nominal_state.with_stm());

    // Build KF without SNC
    let kf = KF::no_snc(initial_estimate.0, msr_noise);

    let prop = Propagator::default(dynamics);
    let prop_est = prop.with(init_sc);

    if (ekf_disable_time.is_some() && ekf_num_meas.is_none())
        || (ekf_disable_time.is_none() && ekf_num_meas.is_some())
    {
        return Err(NyxError::CustomError(format!(
            "For an EKF trigger, you must provide both a disable time and a num measurements."
        )));
    }

    let trigger = match ekf_num_meas {
        Some(ekf_num_meas) => Some(EkfTrigger::new(ekf_num_meas, ekf_disable_time.unwrap())),
        None => None,
    };

    let mut odp = ODProcess::new(prop_est, kf, trigger, resid_crit, Cosm::de438());

    let concrete_arc = arc.to_tracking_arc()?;

    odp.process_arc::<GroundStation>(&concrete_arc)?;

    if let Some(iter_conf) = iter_conf {
        odp.iterate_arc::<GroundStation>(&concrete_arc, iter_conf)?;
    }

    if let Some(epoch) = predict_until {
        let max_step =
            predict_step.ok_or_else(|| NyxError::CustomError("predict_step unset".to_string()))?;
        odp.predict_until(max_step, fixed_step.unwrap_or_else(|| false), epoch)?;
    } else if let Some(duration) = predict_for {
        let max_step =
            predict_step.ok_or_else(|| NyxError::CustomError("predict_step unset".to_string()))?;
        odp.predict_for(max_step, fixed_step.unwrap_or_else(|| false), duration)?;
    }

    let maybe = odp.to_parquet(
        export_path,
        export_cfg.unwrap_or_else(|| ExportCfg::default()),
    );

    match maybe {
        Ok(path) => Ok(format!("{}", path.to_str().unwrap())),
        Err(e) => Err(NyxError::CustomError(e.to_string())),
    }
}

/// Runs an orbit determination prediction-only process given the dynamics, the initial spacecraft object, its orbit estimate, the desired step size, and either a prediction epoch or a prediction duration
/// You must also provide an export path and optionally and export configuration to export the results to a Parquet file.
#[pyfunction]
#[pyo3(
    text_signature = "(dynamics, spacecraft, initial_estimate, step, export_path, export_cfg, predict_until=None, predict_for=None, fixed_step=False)"
)]
pub(crate) fn predictor(
    dynamics: SpacecraftDynamics,
    spacecraft: Spacecraft,
    initial_estimate: OrbitEstimate,
    step: Duration,
    export_path: String,
    export_cfg: Option<ExportCfg>,
    predict_until: Option<Epoch>,
    predict_for: Option<Duration>,
    fixed_step: Option<bool>,
) -> Result<String, NyxError> {
    // TODO: Return a navigation trajectory or use a class that mimics the better ODProcess -- https://github.com/nyx-space/nyx/issues/134
    let msr_noise = Matrix2::from_iterator(vec![1e-10, 0.0, 0.0, 1e-10]);

    let init_sc = spacecraft.with_orbit(initial_estimate.0.nominal_state.with_stm());

    // Build KF without SNC
    let kf = KF::no_snc(initial_estimate.0, msr_noise);

    let prop = Propagator::default(dynamics);
    let prop_est = prop.with(init_sc);

    if (predict_until.is_some() && predict_for.is_some())
        || (predict_until.is_none() && predict_for.is_none())
    {
        return Err(NyxError::CustomError(format!(
            "You must provide either the predict_until argument or the predict_for argument"
        )));
    }

    let mut odp = ODProcess::ckf(prop_est, kf, None, Cosm::de438());

    if let Some(epoch) = predict_until {
        odp.predict_until(step, fixed_step.unwrap_or_else(|| false), epoch)
            .unwrap();
    } else if let Some(duration) = predict_for {
        odp.predict_for(step, fixed_step.unwrap_or_else(|| false), duration)
            .unwrap();
    }

    let maybe = odp.to_parquet(
        export_path,
        export_cfg.unwrap_or_else(|| ExportCfg::default()),
    );

    match maybe {
        Ok(path) => Ok(format!("{}", path.to_str().unwrap())),
        Err(e) => Err(NyxError::CustomError(e.to_string())),
    }
}
