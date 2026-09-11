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

use crate::od::msr::IntegrationRef;
use der::{Decode, Encode, Reader};
use hifitime::{Duration, Unit};
use serde::{Deserialize, Serialize};

#[cfg(feature = "python")]
use pyo3::prelude::*;

#[cfg_attr(feature = "python", pyclass(from_py_object, get_all, set_all))]
#[derive(Copy, Clone, Debug, PartialEq, Eq, Serialize, Deserialize)]
pub struct DopplerConfig {
    /// Duration needed to generate a measurement (if unset, it is assumed to be instantaneous)
    pub integration_time: Duration,
    #[serde(default)]
    /// The integration reference flag as per CCSDS TDM, i.e. if the time tag of the data is
    pub integration_ref: IntegrationRef,
}

impl Default for DopplerConfig {
    fn default() -> Self {
        Self {
            integration_time: 10 * Unit::Second,
            integration_ref: IntegrationRef::Middle,
        }
    }
}

impl<'a> Decode<'a> for DopplerConfig {
    fn decode<R: Reader<'a>>(decoder: &mut R) -> der::Result<Self> {
        let integration_time_ns = decoder.decode()?;
        let integration_time = Duration::from_total_nanoseconds(integration_time_ns);
        let integration_ref = decoder.decode()?;

        Ok(DopplerConfig {
            integration_time,
            integration_ref,
        })
    }
}

impl Encode for DopplerConfig {
    fn encoded_len(&self) -> der::Result<der::Length> {
        self.integration_time.total_nanoseconds().encoded_len()?
            + self.integration_ref.encoded_len()?
    }

    fn encode(&self, encoder: &mut impl der::Writer) -> der::Result<()> {
        let integration_time_ns = self.integration_time.total_nanoseconds();
        integration_time_ns.encode(encoder)?;
        self.integration_ref.encode(encoder)?;

        Ok(())
    }
}
