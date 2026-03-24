use crate::md::prelude::Traj;
use crate::od::msr::measurement::Measurement;
use crate::od::msr::MeasurementType;
use crate::od::NoiseNotConfiguredSnafu;
use crate::od::{ODError, ODTrajSnafu, TrackingDevice};
use crate::time::Epoch;
use crate::Spacecraft;
use anise::errors::AlmanacResult;
use anise::frames::Frame;
use anise::prelude::{Almanac, Orbit};
use indexmap::{IndexMap, IndexSet};
use rand_pcg::Pcg64Mcg;
use snafu::{ensure, ResultExt};
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

    fn location(
        &self,
        _epoch: Epoch,
        _frame: Frame,
        _almanac: Arc<Almanac>,
    ) -> AlmanacResult<Orbit> {
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
        let mut noises = IndexMap::with_capacity(self.measurement_types.len());

        if let Some(rng) = rng {
            ensure!(
                self.stochastic_noises.is_some(),
                NoiseNotConfiguredSnafu {
                    kind: "ground station stochastics".to_string(),
                }
            );

            let stochastics = self.stochastic_noises.as_mut().unwrap();

            for msr_type in &self.measurement_types {
                noises.insert(
                    *msr_type,
                    stochastics
                        .get_mut(msr_type)
                        .ok_or(ODError::NoiseNotConfigured {
                            kind: format!("{msr_type:?}"),
                        })?
                        .sample(rx.orbit.epoch, rng),
                );
            }
        }

        for (ii, msr_type) in self.measurement_types.iter().copied().enumerate() {
            msr.push(
                msr_type,
                rx.orbit.radius_km[ii]
                    + noises.get(&msr_type).unwrap_or(&0.0)
                    + self.measurement_bias(msr_type, rx.orbit.epoch)?,
            );
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
