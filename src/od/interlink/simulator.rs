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
use hifitime::TimeSeries;
use rand_pcg::Pcg64Mcg;

use crate::md::Trajectory;
use crate::od::prelude::TrkConfig;

use std::collections::BTreeMap;

use super::InterlinkTxSpacecraft;

/// Simulates tracking data between the transmitter spacecraft and any number of receiver spacecraft.
/// This can be used to estimate the position of the receiver spacecraft ... I think.
#[derive(Clone)]
pub struct InterlinkArcSim {
    /// Receiver spacecraft in this link
    pub rx_spacecraft: Vec<Trajectory>,
    /// Transmitter spacercaft
    pub tx_spacecraft: InterlinkTxSpacecraft,
    /// Configuration of each device
    pub configs: BTreeMap<String, TrkConfig>,
    /// Random number generator used for this tracking arc, ensures repeatability
    rng: Pcg64Mcg,
    /// Greatest common denominator time series that allows this arc to meet all of the conditions.
    time_series: TimeSeries,
}
