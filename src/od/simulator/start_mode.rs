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

pub use crate::dynamics::{Dynamics, NyxError};
use crate::io::{epoch_from_str, epoch_to_str};
pub use crate::{cosmic::Cosm, State, TimeTagged};
use hifitime::Epoch;
use serde::Deserialize;
use serde_derive::Serialize;
use std::fmt::{Debug, Display};

/// Defines the availability methods for a tracking arc.
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
pub enum Availability {
    /// Start/stop the tracking schedule with respect to when the receiver is visible
    Visible,
    /// Start/stop the tracking schedule with respect to the provided epoch
    #[serde(serialize_with = "epoch_to_str", deserialize_with = "epoch_from_str")]
    Epoch(Epoch),
}

impl Display for Availability {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Availability::Visible => write!(f, "{self:?}"),
            Availability::Epoch(e) => write!(f, "Epoch({e})"),
        }
    }
}

impl Default for Availability {
    fn default() -> Self {
        Self::Visible
    }
}

#[test]
fn serde_startmode() {
    use core::str::FromStr;
    use serde_yaml;

    let vis: Availability = serde_yaml::from_str("!Visible").unwrap();
    assert_eq!(vis, Availability::Visible);

    let val = "!Epoch 2023-02-23T00:00:00 UTC";
    let mode: Availability = serde_yaml::from_str(val).unwrap();
    assert_eq!(
        mode,
        Availability::Epoch(Epoch::from_str("2023-02-23T00:00:00 UTC").unwrap())
    );

    assert_eq!(val, serde_yaml::to_string(&mode).unwrap().trim());
}
