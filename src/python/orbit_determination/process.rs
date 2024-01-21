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
use snafu::ResultExt;

use crate::{
    io::tracking_data::DynamicTrackingArc,
    io::ExportCfg,
    md::prelude::{Cosm, Propagator, SpacecraftDynamics},
    od::{
        filter::kalman::KF,
        process::{EkfTrigger, FltResid, IterationConf, ODIOSnafu, ODProcess},
        snc::SNC3,
        ODError,
    },
    Spacecraft,
};

use super::{estimate::OrbitEstimate, ConfigError, GroundStation};

/// Runs an orbit determination process and returns the path to those results.
#[pyfunction]
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
    iter_conf: Option<IterationConf>,
    snc_disable_time: Option<Duration>,
    snc_diagonals: Option<Vec<f64>>,
) -> Result<String, ODError> {
    let msr_noise = Matrix2::from_iterator(measurement_noise);

    let init_sc = spacecraft.with_orbit(initial_estimate.0.nominal_state.with_stm());

    // Build KF without SNC
    let kf = if (snc_disable_time.is_some() && snc_diagonals.as_ref().is_none())
        || (snc_disable_time.is_none() && snc_diagonals.as_ref().is_some())
        || (snc_diagonals.as_ref().is_some() && snc_diagonals.as_ref().unwrap().len() != 3)
    {
        return Err(ODError::ODConfigError {
            source: ConfigError::InvalidConfig {
                msg: "SNC requires a disable time and the snc_diagonals (3 items required)."
                    .to_string(),
            },
        });
    } else if snc_disable_time.is_some() && snc_diagonals.is_some() {
        let snc = SNC3::from_diagonal(snc_disable_time.unwrap(), &snc_diagonals.unwrap());
        KF::new(initial_estimate.0, snc, msr_noise)
    } else {
        KF::no_snc(initial_estimate.0, msr_noise)
    };

    let prop = Propagator::default(dynamics);
    let prop_est = prop.with(init_sc);

    if (ekf_disable_time.is_some() && ekf_num_meas.is_none())
        || (ekf_disable_time.is_none() && ekf_num_meas.is_some())
    {
        return Err(ODError::ODConfigError {
            source: ConfigError::InvalidConfig {
                msg: "For an EKF trigger, you must provide both a disable time and a num measurements."
                    .to_string(),
            },
        });
    }

    let trigger = match ekf_num_meas {
        Some(ekf_num_meas) => Some(EkfTrigger::new(ekf_num_meas, ekf_disable_time.unwrap())),
        None => None,
    };

    let mut odp = ODProcess::new(prop_est, kf, trigger, resid_crit, Cosm::de438());

    let concrete_arc = arc.to_tracking_arc().with_context(|_| ODIOSnafu)?;

    odp.process_arc::<GroundStation>(&concrete_arc)?;

    if let Some(iter_conf) = iter_conf {
        odp.iterate_arc::<GroundStation>(&concrete_arc, iter_conf)?;
    }

    if let Some(epoch) = predict_until {
        let max_step = predict_step.ok_or_else(|| ODError::ODConfigError {
            source: ConfigError::InvalidConfig {
                msg: "predict_step unset.".to_string(),
            },
        })?;
        odp.predict_until(max_step, epoch)?;
    } else if let Some(duration) = predict_for {
        let max_step = predict_step.ok_or_else(|| ODError::ODConfigError {
            source: ConfigError::InvalidConfig {
                msg: "predict_step unset.".to_string(),
            },
        })?;
        odp.predict_for(max_step, duration)?;
    }

    let path = odp.to_parquet(
        export_path,
        export_cfg.unwrap_or_else(|| ExportCfg::default()),
    )?;

    Ok(format!("{}", path.to_str().unwrap()))
}

/// Runs an orbit determination prediction-only process and returns the path to those results.
#[pyfunction]
pub(crate) fn predictor(
    dynamics: SpacecraftDynamics,
    spacecraft: Spacecraft,
    initial_estimate: OrbitEstimate,
    step: Duration,
    export_path: String,
    export_cfg: Option<ExportCfg>,
    predict_until: Option<Epoch>,
    predict_for: Option<Duration>,
) -> Result<String, ODError> {
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
        return Err(ODError::ODConfigError {
            source: ConfigError::InvalidConfig {
                msg: "predict_step and predict_for unset.".to_string(),
            },
        });
    }

    let mut odp = ODProcess::ckf(prop_est, kf, None, Cosm::de438());

    if let Some(epoch) = predict_until {
        odp.predict_until(step, epoch)?;
    } else if let Some(duration) = predict_for {
        odp.predict_for(step, duration)?;
    }

    let path = odp.to_parquet(
        export_path,
        export_cfg.unwrap_or_else(|| ExportCfg::default()),
    )?;

    Ok(format!("{}", path.to_str().unwrap()))
}
