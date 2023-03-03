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

use crate::io::{epoch_from_str, epoch_to_str, frame_from_str, frame_to_str};
use crate::{cosmic::Frame, Orbit};
use either::Either;
use hifitime::Epoch;
use serde::{Deserialize, Serialize};

#[derive(Copy, Clone, Serialize, Deserialize, Debug)]
#[serde(transparent)]
pub struct OrbitSerde {
    #[serde(with = "either::serde_untagged")]
    inner: Either<Orbit, KeplerianOrbit>,
}

#[derive(Copy, Clone, Serialize, Deserialize, Debug)]
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
    #[serde(serialize_with = "frame_to_str", deserialize_with = "frame_from_str")]
    frame: Frame,
}

impl Into<Orbit> for KeplerianOrbit {
    fn into(self) -> Orbit {
        Orbit::keplerian(
            self.sma_km,
            self.ecc,
            self.inc_deg,
            self.raan_deg,
            self.aop_deg,
            self.ta_deg,
            self.epoch,
            self.frame,
        )
    }
}

impl Into<OrbitSerde> for Orbit {
    fn into(self) -> OrbitSerde {
        OrbitSerde {
            inner: Either::Left(self),
        }
    }
}

impl Into<Orbit> for OrbitSerde {
    fn into(self) -> Orbit {
        match self.inner {
            Either::Left(orbit) => orbit,
            Either::Right(kep) => kep.into(),
        }
    }
}
