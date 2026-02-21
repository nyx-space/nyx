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

use anise::astro::{Aberration, AzElRange, Location, PhysicsResult};
use anise::errors::AlmanacResult;
use anise::prelude::{Almanac, Frame, Orbit};
use der::{Decode, Encode, Reader};
use indexmap::{IndexMap, IndexSet};
use snafu::ensure;

use super::msr::MeasurementType;
use super::noise::{GaussMarkov, StochasticNoise};
use super::{ODAlmanacSnafu, ODError, ODTrajSnafu, TrackingDevice};
use crate::io::ConfigRepr;
use crate::od::NoiseNotConfiguredSnafu;
use crate::time::Epoch;
use hifitime::Duration;
use rand_pcg::Pcg64Mcg;
use serde_derive::{Deserialize, Serialize};
use std::fmt::{self, Debug};

pub mod builtin;
pub mod trk_device;

/// GroundStation defines a two-way ranging and doppler station.
#[derive(Debug, Clone, Serialize, Deserialize, PartialEq)]
pub struct GroundStation {
    pub name: String,
    pub location: Location,
    pub measurement_types: IndexSet<MeasurementType>,
    /// Duration needed to generate a measurement (if unset, it is assumed to be instantaneous)
    pub integration_time: Option<Duration>,
    /// Whether to correct for light travel time
    pub light_time_correction: bool,
    /// Noise on the timestamp of the measurement
    pub timestamp_noise_s: Option<StochasticNoise>,
    pub stochastic_noises: Option<IndexMap<MeasurementType, StochasticNoise>>,
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
            integration_time: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            stochastic_noises: None,
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

    pub fn with_integration_time(mut self, integration_time: Option<Duration>) -> Self {
        self.integration_time = integration_time;

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

    /// Computes the azimuth and elevation of the provided object seen from this ground station, both in degrees.
    /// This is a shortcut to almanac.azimuth_elevation_range_sez.
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
            self.to_orbit(rx.epoch, almanac).unwrap(),
            obstructing_body,
            ab_corr,
        )
    }

    /// Return this ground station as an orbit in its current frame
    pub fn to_orbit(&self, epoch: Epoch, almanac: &Almanac) -> PhysicsResult<Orbit> {
        Orbit::try_latlongalt(
            self.location.latitude_deg,
            self.location.longitude_deg,
            self.location.height_km,
            epoch,
            almanac.frame_info(self.location.frame).unwrap(),
        )
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

    fn available_data(&self) -> u8 {
        let mut bits: u8 = 0;

        if self.integration_time.is_some() {
            bits |= 1 << 0;
        }
        if self.timestamp_noise_s.is_some() {
            bits |= 1 << 1;
        }
        if self.stochastic_noises.is_some() {
            bits |= 1 << 2;
        }
        bits
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
            integration_time: None,
            light_time_correction: false,
            timestamp_noise_s: None,
            stochastic_noises: None,
        }
    }
}

impl ConfigRepr for GroundStation {}

impl<'a> Decode<'a> for GroundStation {
    fn decode<R: Reader<'a>>(decoder: &mut R) -> der::Result<Self> {
        let name: String = decoder.decode()?;
        let location = decoder.decode()?;
        // Measurement types are stored as a sequence of measurement types
        let msr_types_vec: Vec<MeasurementType> = decoder.decode()?;
        let measurement_types = IndexSet::from_iter(msr_types_vec);

        let light_time_correction = decoder.decode()?;

        // The flags tell us what happens next
        let flags: u8 = decoder.decode()?;

        let integration_time = if flags & (1 << 0) != 0 {
            Some(Duration::from_total_nanoseconds(decoder.decode()?))
        } else {
            None
        };

        let timestamp_noise_s = if flags & (1 << 1) != 0 {
            Some(decoder.decode()?)
        } else {
            None
        };

        let stochastic_noises = if flags & (1 << 2) != 0 {
            // Stochastic noises are stored as a sequence of (MeasurementType, StochasticNoise) tuples (SEQUENCE of SEQUENCE)
            // We define a helper struct for decoding
            #[derive(der::Sequence)]
            struct MsrNoisePair {
                msr_type: MeasurementType,
                noise: StochasticNoise,
            }

            let stochastics_vec: Vec<MsrNoisePair> = decoder.decode()?;
            let mut map = IndexMap::new();
            for pair in stochastics_vec {
                map.insert(pair.msr_type, pair.noise);
            }
            Some(map)
        } else {
            None
        };

        Ok(GroundStation {
            name,
            location,
            measurement_types,
            integration_time,
            light_time_correction,
            timestamp_noise_s,
            stochastic_noises,
        })
    }
}

impl Encode for GroundStation {
    fn encoded_len(&self) -> der::Result<der::Length> {
        let msr_types_vec: Vec<MeasurementType> = self.measurement_types.iter().copied().collect();

        let integration_time_ns = self.integration_time.map(|d| d.total_nanoseconds());

        // Helper for encoding map
        #[derive(der::Sequence)]
        struct MsrNoisePair {
            msr_type: MeasurementType,
            noise: StochasticNoise,
        }

        let stochastics_vec = self.stochastic_noises.as_ref().map(|map| {
            map.iter()
                .map(|(k, v)| MsrNoisePair {
                    msr_type: *k,
                    noise: *v,
                })
                .collect::<Vec<MsrNoisePair>>()
        });

        self.name.encoded_len()?
            + self.location.encoded_len()?
            + msr_types_vec.encoded_len()?
            + self.light_time_correction.encoded_len()?
            + self.available_data().encoded_len()?
            + integration_time_ns.encoded_len()?
            + self.timestamp_noise_s.encoded_len()?
            + stochastics_vec.encoded_len()?
    }

    fn encode(&self, encoder: &mut impl der::Writer) -> der::Result<()> {
        self.name.encode(encoder)?;
        self.location.encode(encoder)?;

        let msr_types_vec: Vec<MeasurementType> = self.measurement_types.iter().cloned().collect();
        msr_types_vec.encode(encoder)?;

        self.light_time_correction.encode(encoder)?;
        self.available_data().encode(encoder)?;

        let integration_time_ns = self.integration_time.map(|d| d.total_nanoseconds());
        integration_time_ns.encode(encoder)?;

        self.timestamp_noise_s.encode(encoder)?;

        // Helper for encoding map
        #[derive(der::Sequence)]
        struct MsrNoisePair {
            msr_type: MeasurementType,
            noise: StochasticNoise,
        }

        let stochastics_vec = self.stochastic_noises.as_ref().map(|map| {
            map.iter()
                .map(|(k, v)| MsrNoisePair {
                    msr_type: *k,
                    noise: *v,
                })
                .collect::<Vec<MsrNoisePair>>()
        });
        stochastics_vec.encode(encoder)?;

        Ok(())
    }
}

impl fmt::Display for GroundStation {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "{} ({})", self.name, self.location,)
    }
}

#[cfg(test)]
mod gs_ut {

    use anise::astro::{Location, TerrainMask};
    use anise::constants::frames::IAU_EARTH_FRAME;
    use indexmap::{IndexMap, IndexSet};

    use crate::io::ConfigRepr;
    use crate::od::prelude::*;

    #[test]
    fn test_load_single() {
        use std::env;
        use std::path::PathBuf;

        use hifitime::TimeUnits;

        let test_data: PathBuf = [
            env::var("CARGO_MANIFEST_DIR").unwrap(),
            "data".to_string(),
            "03_tests".to_string(),
            "config".to_string(),
            "one_ground_station.yaml".to_string(),
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
            integration_time: Some(60 * Unit::Second),
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
            env::var("CARGO_MANIFEST_DIR").unwrap(),
            "data".to_string(),
            "03_tests".to_string(),
            "config".to_string(),
            "many_ground_stations.yaml".to_string(),
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
                integration_time: None,
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
                integration_time: None,
            },
        ];

        assert_eq!(expected, stations);

        // Serialize back
        let reser = serde_yml::to_string(&expected).unwrap();
        dbg!(reser);
    }
}
