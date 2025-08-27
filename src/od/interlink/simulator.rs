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
use anise::almanac::Almanac;
use anise::astro::Aberration;
use anise::errors::AlmanacResult;
use anise::prelude::{Frame, Orbit};
use hifitime::{Duration, Epoch, TimeSeries, TimeUnits};
use indexmap::{IndexMap, IndexSet};
use rand_pcg::Pcg64Mcg;
use serde::{Deserialize, Serialize};
use snafu::{ensure, ResultExt};

use crate::io::ConfigRepr;
use crate::md::prelude::Traj;
use crate::md::Trajectory;
use crate::od::msr::MeasurementType;
use crate::od::noise::StochasticNoise;
use crate::od::prelude::{Measurement, NoiseNotConfiguredSnafu, ODError, TrkConfig};
use crate::od::TrackingDevice;
use crate::od::{ODAlmanacSnafu, ODTrajSnafu};
use crate::Spacecraft;
use crate::State;

use std::collections::BTreeMap;
use std::sync::Arc;

use super::InterlinkTxSpacecraft;

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
