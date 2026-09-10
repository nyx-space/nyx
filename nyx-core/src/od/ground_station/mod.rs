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

use anise::astro::{Aberration, AzElRange, Location};
use anise::errors::{AlmanacError, AlmanacResult};
use anise::frames::FrameUid;
use anise::prelude::{Almanac, Frame, Orbit};
use indexmap::{IndexMap, IndexSet};
use snafu::ensure;

use super::msr::MeasurementType;
use super::noise::{GaussMarkov, StochasticNoise};
use super::{ODAlmanacSnafu, ODError, ODTrajSnafu, TrackingDevice};
use crate::od::NoiseNotConfiguredSnafu;
use crate::time::Epoch;
use rand_pcg::Pcg64Mcg;
use serde::{Deserialize, Serialize};
use std::fmt::{self, Debug};

mod asn1;
pub mod builtin;
mod doppler_config;
pub mod trk_device;

pub use doppler_config::DopplerConfig;

#[cfg(feature = "python")]
use pyo3::exceptions::PyValueError;
#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
use pyo3::types::{PyBytes, PyType};
#[cfg(feature = "python")]
mod python;

#[cfg(feature = "python")]
use der::{Decode, Encode};
/// GroundStation defines a one-way or two-way ranging and doppler station. Set the doppler config for two-way.
///
/// :type name: str
/// :type location: Location
/// :type stochastic_noises: dict[MeasurementType, StochasticNoise]
/// :type doppler_config: DopplerConfig | None
/// :type light_time_correction: bool | None
/// :type timestamp_noise_s: StochasticNoise | None
/// :type obstruction_body: FrameUid | None
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
#[cfg_attr(feature = "python", pyclass(from_py_object))]
pub struct GroundStation {
    pub name: String,
    pub location: Location,
    pub measurement_types: IndexSet<MeasurementType>,
    /// Doppler tracking loop settings (required if tracking Doppler)
    #[serde(
        default = "default_doppler_config",
        skip_serializing_if = "Option::is_none"
    )]
    pub doppler_config: Option<DopplerConfig>,
    /// If light-time correction is enabled, then Range and Doppler are assumed coherent Two-Way; Az/El is OneWay.
    pub light_time_correction: bool,
    /// Noise on the timestamp of the measurement
    pub timestamp_noise_s: Option<StochasticNoise>,
    pub stochastic_noises: Option<IndexMap<MeasurementType, StochasticNoise>>,
    /// Body that obstructs the line of sight (e.g. Moon if tracking a lunar spacecraft from Earth)
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub obstruction_body: Option<FrameUid>,
}

#[cfg_attr(feature = "python", pymethods)]
impl GroundStation {
    /// Computes the azimuth and elevation of the provided object seen from this ground station, both in degrees.
    /// This is a shortcut to almanac.azimuth_elevation_range_sez.
    ///
    /// :type rx: Orbit
    /// :type obstructing_body: Frame | None
    /// :type almanac: Almanac
    /// :rtype: AzElRange
    pub fn azimuth_elevation_of(
        &self,
        rx: Orbit,
        obstructing_body: Option<Frame>,
        almanac: &Almanac,
    ) -> AlmanacResult<AzElRange> {
        let ab_corr = if self.light_time_correction {
            Aberration::LT
        } else {
            Aberration::NONE
        };
        almanac.azimuth_elevation_range_sez(
            rx,
            self.to_orbit(rx.epoch, almanac)?,
            obstructing_body,
            ab_corr,
        )
    }

    /// Return this ground station as an orbit in its current frame
    ///
    /// :type epoch: Epoch
    /// :type almanac: Almanac
    /// :rtype: Orbit
    pub fn to_orbit(&self, epoch: Epoch, almanac: &Almanac) -> AlmanacResult<Orbit> {
        Orbit::try_latlongalt(
            self.location.latitude_deg,
            self.location.longitude_deg,
            self.location.height_km,
            epoch,
            almanac.frame_info(self.location.frame).map_err(|source| {
                AlmanacError::GenericError {
                    err: source.to_string(),
                }
            })?,
        )
        .map_err(|source| AlmanacError::AlmanacPhysics {
            action: "building ground station location",
            source: Box::new(source),
        })
    }
}

impl GroundStation {
    /// Initializes a point on the surface of a celestial object.
    /// This is meant for analysis, not for spacecraft navigation.
    pub fn from_point(
        name: String,
        latitude_deg: f64,
        longitude_deg: f64,
        height_km: f64,
        frame: Frame,
    ) -> Self {
        Self {
            name,
            location: Location {
                latitude_deg,
                longitude_deg,
                height_km,
                frame: frame.into(),
                terrain_mask: vec![],
                terrain_mask_ignored: true,
            },
            measurement_types: IndexSet::new(),
            doppler_config: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            stochastic_noises: None,
            obstruction_body: None,
        }
    }

    /// Returns a copy of this ground station with the new measurement type added (or replaced)
    pub fn with_msr_type(mut self, msr_type: MeasurementType, noise: StochasticNoise) -> Self {
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

    /// Returns a copy of this ground station without the provided measurement type (if defined, else no error)
    pub fn without_msr_type(mut self, msr_type: MeasurementType) -> Self {
        if let Some(noises) = self.stochastic_noises.as_mut() {
            noises.swap_remove(&msr_type);
        }

        self.measurement_types.swap_remove(&msr_type);

        self
    }

    pub fn with_doppler_config(mut self, doppler_config: Option<DopplerConfig>) -> Self {
        self.doppler_config = doppler_config;

        self
    }

    pub fn with_obstruction_body(mut self, obstruction_body: Option<FrameUid>) -> Self {
        self.obstruction_body = obstruction_body;

        self
    }

    /// Returns a copy of this ground station with the measurement type noises' constant bias set to the provided value.
    pub fn with_msr_bias_constant(
        mut self,
        msr_type: MeasurementType,
        bias_constant: f64,
    ) -> Result<Self, ODError> {
        if self.stochastic_noises.is_none() {
            self.stochastic_noises = Some(IndexMap::new());
        }

        let stochastics = self.stochastic_noises.as_mut().unwrap();

        let this_noise = stochastics
            .get_mut(&msr_type)
            .ok_or(ODError::NoiseNotConfigured {
                kind: format!("{msr_type:?}"),
            })
            .unwrap();

        if this_noise.bias.is_none() {
            this_noise.bias = Some(GaussMarkov::ZERO);
        }

        this_noise.bias.unwrap().constant = Some(bias_constant);

        Ok(self)
    }

    /// Returns the noises for all measurement types configured for this ground station at the provided epoch, timestamp noise is the first entry.
    fn noises(&mut self, epoch: Epoch, rng: Option<&mut Pcg64Mcg>) -> Result<Vec<f64>, ODError> {
        let mut noises = vec![0.0; self.measurement_types.len() + 1];

        if let Some(rng) = rng {
            ensure!(
                self.stochastic_noises.is_some(),
                NoiseNotConfiguredSnafu {
                    kind: "ground station stochastics".to_string(),
                }
            );
            // Add the timestamp noise first

            if let Some(mut timestamp_noise) = self.timestamp_noise_s {
                noises[0] = timestamp_noise.sample(epoch, rng);
            }

            let stochastics = self.stochastic_noises.as_mut().unwrap();

            for (ii, msr_type) in self.measurement_types.iter().enumerate() {
                noises[ii + 1] = stochastics
                    .get_mut(msr_type)
                    .ok_or(ODError::NoiseNotConfigured {
                        kind: format!("{msr_type:?}"),
                    })?
                    .sample(epoch, rng);
            }
        }

        Ok(noises)
    }

    pub(crate) fn available_data(&self) -> u8 {
        let mut bits: u8 = 0;

        if self.doppler_config.is_some() {
            bits |= 1 << 0;
        }
        if self.timestamp_noise_s.is_some() {
            bits |= 1 << 1;
        }
        if self.stochastic_noises.is_some() {
            bits |= 1 << 2;
        }
        if self.obstruction_body.is_some() {
            bits |= 1 << 3;
        }
        bits
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl GroundStation {
    /// Decodes an ASN.1 DER encoded byte array into a GroundStation object.
    ///
    /// :type data: bytes
    /// :rtype: GroundStation
    #[classmethod]
    pub fn from_asn1(_cls: &Bound<'_, PyType>, data: &[u8]) -> PyResult<Self> {
        match Self::from_der(data) {
            Ok(obj) => Ok(obj),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 decoding error: {e}"))),
        }
    }

    /// Encodes this GroundStation object into an ASN.1 DER encoded byte array.
    ///
    /// :rtype: bytes
    pub fn to_asn1<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let mut buf = Vec::new();
        match self.encode_to_vec(&mut buf) {
            Ok(_) => Ok(PyBytes::new(py, &buf)),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 encoding error: {e}"))),
        }
    }
}

impl Default for GroundStation {
    fn default() -> Self {
        let mut measurement_types = IndexSet::new();
        measurement_types.insert(MeasurementType::Range);
        measurement_types.insert(MeasurementType::Doppler);
        Self {
            name: "UNDEFINED".to_string(),
            measurement_types,
            location: Location::default(),
            doppler_config: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            stochastic_noises: None,
            obstruction_body: None,
        }
    }
}

impl fmt::Display for GroundStation {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} ({})", self.name, self.location)
    }
}

fn default_doppler_config() -> Option<DopplerConfig> {
    Some(DopplerConfig::default())
}

#[cfg(test)]
mod gs_ut {

    use anise::astro::{Location, TerrainMask};
    use anise::constants::frames::IAU_EARTH_FRAME;
    use indexmap::{IndexMap, IndexSet};

    use crate::io::ConfigRepr;
    use crate::od::prelude::*;

    #[ignore = "github cache is a pain, works locally"]
    #[test]
    fn test_load_single() {
        use std::env;
        use std::path::PathBuf;

        use hifitime::TimeUnits;

        let test_data: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "../data",
            "03_tests",
            "config",
            "one_ground_station.yaml",
        ]
        .iter()
        .collect();

        assert!(test_data.exists(), "Could not find the test data");

        let gs = GroundStation::load(test_data).unwrap();

        dbg!(&gs);

        let mut measurement_types = IndexSet::new();
        measurement_types.insert(MeasurementType::Range);
        measurement_types.insert(MeasurementType::Doppler);

        let mut stochastics = IndexMap::new();
        stochastics.insert(
            MeasurementType::Range,
            StochasticNoise {
                bias: Some(GaussMarkov::new(1.days(), 5e-3).unwrap()),
                ..Default::default()
            },
        );
        stochastics.insert(
            MeasurementType::Doppler,
            StochasticNoise {
                bias: Some(GaussMarkov::new(1.days(), 5e-5).unwrap()),
                ..Default::default()
            },
        );

        let expected_gs = GroundStation {
            name: "Demo ground station".to_string(),
            location: Location {
                latitude_deg: 2.3522,
                longitude_deg: 48.8566,
                height_km: 0.4,
                frame: IAU_EARTH_FRAME.into(),
                terrain_mask: TerrainMask::from_flat_terrain(5.0),
                terrain_mask_ignored: false,
            },
            measurement_types,
            stochastic_noises: Some(stochastics),

            light_time_correction: false,
            timestamp_noise_s: None,
            doppler_config: Some(DopplerConfig::default()),
            obstruction_body: None,
        };

        println!("{}", serde_yml::to_string(&expected_gs).unwrap());

        assert_eq!(expected_gs, gs);
    }

    #[test]
    fn test_load_many() {
        use hifitime::TimeUnits;
        use std::env;
        use std::path::PathBuf;

        let test_file: PathBuf = [
            env!("CARGO_MANIFEST_DIR"),
            "../data",
            "03_tests",
            "config",
            "many_ground_stations.yaml",
        ]
        .iter()
        .collect();

        let stations = GroundStation::load_many(test_file).unwrap();

        dbg!(&stations);

        let mut measurement_types = IndexSet::new();
        measurement_types.insert(MeasurementType::Range);
        measurement_types.insert(MeasurementType::Doppler);

        let mut stochastics = IndexMap::new();
        stochastics.insert(
            MeasurementType::Range,
            StochasticNoise {
                bias: Some(GaussMarkov::new(1.days(), 5e-3).unwrap()),
                ..Default::default()
            },
        );
        stochastics.insert(
            MeasurementType::Doppler,
            StochasticNoise {
                bias: Some(GaussMarkov::new(1.days(), 5e-5).unwrap()),
                ..Default::default()
            },
        );

        let expected = vec![
            GroundStation {
                name: "Demo ground station".to_string(),
                location: Location {
                    latitude_deg: 2.3522,
                    longitude_deg: 48.8566,
                    height_km: 0.4,
                    frame: IAU_EARTH_FRAME.into(),
                    terrain_mask: TerrainMask::from_flat_terrain(5.0),
                    terrain_mask_ignored: false,
                },
                measurement_types: measurement_types.clone(),
                stochastic_noises: Some(stochastics.clone()),
                light_time_correction: false,
                timestamp_noise_s: None,
                doppler_config: None,
                obstruction_body: None,
            },
            GroundStation {
                name: "Canberra".to_string(),
                location: Location {
                    latitude_deg: -35.398333,
                    longitude_deg: 148.981944,
                    height_km: 0.691750,
                    frame: IAU_EARTH_FRAME.into(),
                    terrain_mask: TerrainMask::from_flat_terrain(5.0),
                    terrain_mask_ignored: false,
                },
                measurement_types,
                stochastic_noises: Some(stochastics),
                light_time_correction: false,
                timestamp_noise_s: None,
                doppler_config: None,
                obstruction_body: None,
            },
        ];

        assert_eq!(expected, stations);

        // Serialize back
        let reser = serde_yml::to_string(&expected).unwrap();
        dbg!(reser);
    }
}
