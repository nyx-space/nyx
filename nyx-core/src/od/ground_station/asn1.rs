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

use super::GroundStation;
use crate::io::ConfigRepr;
use crate::od::msr::MeasurementType;
use crate::od::noise::StochasticNoise;
use der::{Decode, Encode, Reader};
use indexmap::{IndexMap, IndexSet};

impl ConfigRepr for GroundStation {}

#[derive(der::Sequence)]
struct MsrNoisePair {
    msr_type: MeasurementType,
    noise: StochasticNoise,
}

impl<'a> Decode<'a> for GroundStation {
    fn decode<R: Reader<'a>>(decoder: &mut R) -> der::Result<Self> {
        let name: String = decoder.decode()?;
        let location = decoder.decode()?;
        // Measurement types are stored as a sequence of measurement types
        let msr_types_vec: Vec<MeasurementType> = decoder.decode()?;
        let measurement_types = IndexSet::from_iter(msr_types_vec);

        let light_time_correction = decoder.decode()?;
        let relativistic_corrections = decoder.decode()?;

        // The flags tell us what happens next
        let flags: u8 = decoder.decode()?;

        let doppler_config = if flags & (1 << 0) != 0 {
            Some(decoder.decode()?)
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

            let stochastics_vec: Vec<MsrNoisePair> = decoder.decode()?;
            let mut map = IndexMap::new();
            for pair in stochastics_vec {
                map.insert(pair.msr_type, pair.noise);
            }
            Some(map)
        } else {
            None
        };

        let obstruction_body = if flags & (1 << 3) != 0 {
            Some(decoder.decode()?)
        } else {
            None
        };

        Ok(GroundStation {
            name,
            location,
            measurement_types,
            doppler_config,
            light_time_correction,
            timestamp_noise_s,
            stochastic_noises,
            obstructing_body: obstruction_body,
            relativistic_corrections,
        })
    }
}

impl Encode for GroundStation {
    fn encoded_len(&self) -> der::Result<der::Length> {
        let msr_types_vec: Vec<MeasurementType> = self.measurement_types.iter().copied().collect();

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
            + self.relativistic_corrections.encoded_len()?
            + self.available_data().encoded_len()?
            + self.doppler_config.encoded_len()?
            + self.timestamp_noise_s.encoded_len()?
            + stochastics_vec.encoded_len()?
            + self.obstructing_body.encoded_len()?
    }

    fn encode(&self, encoder: &mut impl der::Writer) -> der::Result<()> {
        self.name.encode(encoder)?;
        self.location.encode(encoder)?;

        let msr_types_vec: Vec<MeasurementType> = self.measurement_types.iter().copied().collect();
        msr_types_vec.encode(encoder)?;

        self.light_time_correction.encode(encoder)?;
        self.relativistic_corrections.encode(encoder)?;
        self.available_data().encode(encoder)?;

        self.doppler_config.encode(encoder)?;
        self.timestamp_noise_s.encode(encoder)?;

        let stochastics_vec = self.stochastic_noises.as_ref().map(|map| {
            map.iter()
                .map(|(k, v)| MsrNoisePair {
                    msr_type: *k,
                    noise: *v,
                })
                .collect::<Vec<MsrNoisePair>>()
        });
        stochastics_vec.encode(encoder)?;

        self.obstructing_body.encode(encoder)?;

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::od::ground_station::DopplerConfig;
    use crate::od::msr::IntegrationRef;
    use crate::od::noise::{GaussMarkov, WhiteNoise};
    use anise::astro::{Location, TerrainMask};
    use anise::constants::frames::IAU_EARTH_FRAME;
    use hifitime::{TimeUnits, Unit};

    #[test]
    fn test_ground_station_asn1_roundtrip() {
        let mut measurement_types = IndexSet::new();
        measurement_types.insert(MeasurementType::Range);
        measurement_types.insert(MeasurementType::Doppler);

        let mut stochastics = IndexMap::new();
        stochastics.insert(
            MeasurementType::Range,
            StochasticNoise {
                bias: Some(GaussMarkov::new(1.days(), 5e-3).unwrap()),
                white_noise: Some(WhiteNoise::constant_white_noise(1e-3)),
            },
        );
        stochastics.insert(
            MeasurementType::Doppler,
            StochasticNoise {
                bias: Some(GaussMarkov::new(1.days(), 5e-5).unwrap()),
                white_noise: None,
            },
        );

        let mut gs = GroundStation {
            name: "Test Station".to_string(),
            location: Location {
                latitude_deg: 40.427_222,
                longitude_deg: 4.250_556,
                height_km: 0.834_939,
                frame: IAU_EARTH_FRAME.into(),
                terrain_mask: TerrainMask::from_flat_terrain(5.0),
                terrain_mask_ignored: false,
            },
            measurement_types: measurement_types.clone(),
            doppler_config: None,
            light_time_correction: true,
            timestamp_noise_s: None,
            stochastic_noises: None,
            obstructing_body: None,
            relativistic_corrections: false,
        };

        // 1. Minimal GroundStation (all optional fields None)
        let mut buf = vec![];
        gs.encode_to_vec(&mut buf).unwrap();
        let decoded = GroundStation::from_der(&buf).unwrap();
        assert_eq!(decoded, gs);

        // 2. With DopplerConfig default
        gs.doppler_config = Some(DopplerConfig::default());
        buf.clear();
        gs.encode_to_vec(&mut buf).unwrap();
        let decoded = GroundStation::from_der(&buf).unwrap();
        assert_eq!(decoded, gs);

        // 3. With DopplerConfig custom
        gs.doppler_config = Some(DopplerConfig {
            integration_time: 10 * Unit::Second,
            integration_ref: IntegrationRef::Start,
        });
        buf.clear();
        gs.encode_to_vec(&mut buf).unwrap();
        let decoded = GroundStation::from_der(&buf).unwrap();
        assert_eq!(decoded, gs);

        // 4. With timestamp_noise_s
        gs.timestamp_noise_s = Some(StochasticNoise {
            bias: Some(GaussMarkov::new(10 * Unit::Second, 1e-9).unwrap()),
            white_noise: None,
        });
        buf.clear();
        gs.encode_to_vec(&mut buf).unwrap();
        let decoded = GroundStation::from_der(&buf).unwrap();
        assert_eq!(decoded, gs);

        // 5. With stochastic_noises
        gs.stochastic_noises = Some(stochastics);
        buf.clear();
        gs.encode_to_vec(&mut buf).unwrap();
        let decoded = GroundStation::from_der(&buf).unwrap();
        assert_eq!(decoded, gs);

        // 6. With only stochastic_noises (doppler_config and timestamp_noise_s unset)
        gs.doppler_config = None;
        gs.timestamp_noise_s = None;
        buf.clear();
        gs.encode_to_vec(&mut buf).unwrap();
        let decoded = GroundStation::from_der(&buf).unwrap();
        assert_eq!(decoded, gs);

        // 7. With obstruction_body
        gs.obstructing_body = Some(IAU_EARTH_FRAME.into());
        buf.clear();
        gs.encode_to_vec(&mut buf).unwrap();
        let decoded = GroundStation::from_der(&buf).unwrap();
        assert_eq!(decoded, gs);

        // 8. With relativistic_corrections
        gs.relativistic_corrections = true;
        buf.clear();
        gs.encode_to_vec(&mut buf).unwrap();
        let decoded = GroundStation::from_der(&buf).unwrap();
        assert_eq!(decoded, gs);
    }
}
