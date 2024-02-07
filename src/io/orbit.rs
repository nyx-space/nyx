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

use either::Either;
use hifitime::Epoch;
use serde::{Deserialize, Serialize};
use serde_dhall::StaticType;

use crate::io::{epoch_from_str, epoch_to_str};
use crate::{cosmic::Frame, Orbit};

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize, StaticType)]
#[serde(transparent)]
pub struct OrbitSerde {
    #[serde(with = "either::serde_untagged")]
    inner: Either<Orbit, KeplerianOrbit>,
}

#[derive(Copy, Clone, Debug, PartialEq, Serialize, Deserialize, StaticType)]
pub struct KeplerianOrbit {
    sma_km: f64,
    ecc: f64,
    inc_deg: f64,
    raan_deg: f64,
    aop_deg: f64,
    ta_deg: f64,
    #[serde(serialize_with = "epoch_to_str", deserialize_with = "epoch_from_str")]
    epoch: Epoch,
    /// Frame contains everything we need to compute state information
    frame: Frame,
}

impl From<KeplerianOrbit> for Orbit {
    fn from(val: KeplerianOrbit) -> Self {
        Orbit::keplerian(
            val.sma_km,
            val.ecc,
            val.inc_deg,
            val.raan_deg,
            val.aop_deg,
            val.ta_deg,
            val.epoch,
            val.frame,
        )
    }
}

impl From<Orbit> for OrbitSerde {
    fn from(val: Orbit) -> Self {
        OrbitSerde {
            inner: Either::Left(val),
        }
    }
}

impl From<OrbitSerde> for Orbit {
    fn from(val: OrbitSerde) -> Orbit {
        match val.inner {
            Either::Left(orbit) => orbit,
            Either::Right(kep) => kep.into(),
        }
    }
}

#[cfg(test)]
mod ut_serde_orbit {
    use anise::constants::frames::EARTH_J2000;
    use serde_dhall;

    use super::{Epoch, Orbit, OrbitSerde};

    #[test]
    fn serde_cartesian() {
        let orbit = Orbit::keplerian(
            8159.0,
            0.001,
            38.6,
            95.1,
            82.3,
            45.6,
            Epoch::from_gregorian_utc_at_midnight(2022, 2, 29),
            EARTH_J2000,
        );
        let as_serde: OrbitSerde = orbit.into();

        let serialized = serde_dhall::serialize(&as_serde)
            .static_type_annotation()
            .to_string()
            .unwrap();
        println!("{serialized}");
    }
}
