pub mod sensitivity;
pub mod trk_device;

#[cfg(feature = "python")]
pub mod python;

use crate::io::ConfigRepr;
use crate::od::msr::MeasurementType;
use crate::od::noise::StochasticNoise;
use anise::prelude::Frame;
use indexmap::IndexMap;
use indexmap::IndexSet;
use serde::{Deserialize, Serialize};
use std::fmt::Display;

#[cfg(feature = "python")]
use pyo3::prelude::*;

/// Position device can be used to post-filter position measurements from GNSS/GPS devices.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass(from_py_object))]
pub struct PositionDevice {
    pub name: String,
    pub frame: Option<Frame>,
    pub stochastic_noises: Option<IndexMap<MeasurementType, StochasticNoise>>,
    pub measurement_types: IndexSet<MeasurementType>,
}

impl PositionDevice {
    pub fn new(name: String, frame: Option<Frame>) -> Self {
        Self {
            name,
            frame,
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
