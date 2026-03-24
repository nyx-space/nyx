pub mod sensitivity;
pub mod trk_device;

use crate::io::ConfigRepr;
use crate::od::msr::MeasurementType;
use crate::od::noise::StochasticNoise;
use indexmap::IndexMap;
use indexmap::IndexSet;
use serde::{Deserialize, Serialize};
use std::fmt::Display;

/// Position device can be used to post-filter position measurements from GNSS/GPS devices.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct PositionDevice {
    pub name: String,
    pub stochastic_noises: Option<IndexMap<MeasurementType, StochasticNoise>>,
    pub measurement_types: IndexSet<MeasurementType>,
}

impl PositionDevice {
    pub fn new(name: String) -> Self {
        Self {
            name,
            stochastic_noises: None,
            measurement_types: IndexSet::new(),
        }
    }

    pub fn with_noise(mut self, msr_type: MeasurementType, noise: StochasticNoise) -> Self {
        if self.stochastic_noises.is_none() {
            self.stochastic_noises = Some(IndexMap::new());
        }
        self.stochastic_noises
            .as_mut()
            .unwrap()
            .insert(msr_type, noise);
        self.measurement_types.insert(msr_type);
        self
    }
}

impl ConfigRepr for PositionDevice {}

impl Display for PositionDevice {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "PositionDevice({})", self.name)
    }
}
