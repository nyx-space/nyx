use crate::md::prelude::Traj;
use crate::od::msr::measurement::Measurement;
use crate::od::msr::MeasurementType;
use crate::od::{ODError, ODTrajSnafu, TrackingDevice};
use crate::time::Epoch;
use crate::Spacecraft;
use anise::errors::AlmanacResult;
use anise::frames::Frame;
use anise::prelude::{Almanac, Orbit};
use indexmap::IndexSet;
use rand_pcg::Pcg64Mcg;
use snafu::ResultExt;
use std::sync::Arc;

use super::XyzDevice;

impl TrackingDevice<Spacecraft> for XyzDevice {
    fn measurement_types(&self) -> &IndexSet<MeasurementType> {
        &self.measurement_types
    }

    fn measure(
        &mut self,
        epoch: Epoch,
        traj: &Traj<Spacecraft>,
        rng: Option<&mut Pcg64Mcg>,
        almanac: Arc<Almanac>,
    ) -> Result<Option<Measurement>, ODError> {
        let rx = traj.at(epoch).context(ODTrajSnafu {
            details: "fetching state for instantaneous measurement".to_string(),
        })?;
        self.measure_instantaneous(rx, rng, almanac)
    }

    fn name(&self) -> String {
        self.name.clone()
    }

    fn location(&self, _epoch: Epoch, _frame: Frame, _almanac: Arc<Almanac>) -> AlmanacResult<Orbit> {
        // XyzDevice does not have a location
        unimplemented!("XyzDevice does not have a location")
    }

    fn measure_instantaneous(
        &mut self,
        rx: Spacecraft,
        rng: Option<&mut Pcg64Mcg>,
        _almanac: Arc<Almanac>,
    ) -> Result<Option<Measurement>, ODError> {
        let mut msr = Measurement::new(self.name.clone(), rx.orbit.epoch);

        let mut noise_x = 0.0;
        let mut noise_y = 0.0;
        let mut noise_z = 0.0;

        if let Some(rng_gen) = rng {
            use rand_distr::{Distribution, Normal};

            if self.measurement_types.contains(&MeasurementType::X) {
                let cov = self.measurement_covar(MeasurementType::X, rx.orbit.epoch)?;
                let normal = Normal::new(0.0, cov.sqrt()).unwrap();
                noise_x = normal.sample(rng_gen);
            }
            if self.measurement_types.contains(&MeasurementType::Y) {
                let cov = self.measurement_covar(MeasurementType::Y, rx.orbit.epoch)?;
                let normal = Normal::new(0.0, cov.sqrt()).unwrap();
                noise_y = normal.sample(rng_gen);
            }
            if self.measurement_types.contains(&MeasurementType::Z) {
                let cov = self.measurement_covar(MeasurementType::Z, rx.orbit.epoch)?;
                let normal = Normal::new(0.0, cov.sqrt()).unwrap();
                noise_z = normal.sample(rng_gen);
            }
        }

        if self.measurement_types.contains(&MeasurementType::X) {
            msr.push(MeasurementType::X, rx.orbit.radius_km.x + noise_x + self.measurement_bias(MeasurementType::X, rx.orbit.epoch)?);
        }
        if self.measurement_types.contains(&MeasurementType::Y) {
            msr.push(MeasurementType::Y, rx.orbit.radius_km.y + noise_y + self.measurement_bias(MeasurementType::Y, rx.orbit.epoch)?);
        }
        if self.measurement_types.contains(&MeasurementType::Z) {
            msr.push(MeasurementType::Z, rx.orbit.radius_km.z + noise_z + self.measurement_bias(MeasurementType::Z, rx.orbit.epoch)?);
        }

        Ok(Some(msr))
    }

    fn measurement_covar(&self, msr_type: MeasurementType, epoch: Epoch) -> Result<f64, ODError> {
        if let Some(stochastics) = &self.stochastic_noises {
            Ok(stochastics
                .get(&msr_type)
                .ok_or(ODError::NoiseNotConfigured {
                    kind: format!("{msr_type:?}"),
                })?
                .covariance(epoch))
        } else {
            Ok(0.0)
        }
    }

    fn measurement_bias(&self, msr_type: MeasurementType, _epoch: Epoch) -> Result<f64, ODError> {
        if let Some(stochastics) = &self.stochastic_noises {
            if let Some(gm) = stochastics
                .get(&msr_type)
                .ok_or(ODError::NoiseNotConfigured {
                    kind: format!("{msr_type:?}"),
                })?
                .bias
            {
                Ok(gm.constant.unwrap_or(0.0))
            } else {
                Ok(0.0)
            }
        } else {
            Ok(0.0)
        }
    }
}
