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

use anise::constants::SPEED_OF_LIGHT_KM_S;
use anise::constants::frames::ICRF;
use anise::prelude::Almanac;
use hifitime::{Duration, Epoch, TimeUnits};
use nalgebra::Vector3;

use crate::Spacecraft;
use crate::md::prelude::Traj;
use crate::od::ground_station::GroundStation;
use crate::od::{ODAlmanacSnafu, ODError, ODTrajSnafu};
use snafu::ResultExt;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LightTimeLeg {
    pub t_start: Epoch,
    pub t_end: Epoch,
    pub r_start_icrf: Vector3<f64>,
    pub r_end_icrf: Vector3<f64>,
    pub light_time: Duration,
}

impl LightTimeLeg {
    pub fn range_km(&self) -> f64 {
        (self.r_start_icrf - self.r_end_icrf).norm()
    }
}

/// Solves the two-way tracking problem using the Picard fixed-point iteration like in JPL Moyer (2000).
///
/// This solves the problem strictly in ICRF because it is the only true inertial frame.
///
/// ```text
/// t1 (Transmit)                 t2 (Bounce)                 t3 (Receive)
/// Ground Station -------------> Spacecraft -------------> Ground Station
///               [Uplink Leg]                [Downlink Leg]
///               (Solved Second)              (Solved First)
/// ```
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct TwoWaySolution {
    pub t1_transmit: Epoch,
    pub t2_bounce: Epoch,
    pub t3_receive: Epoch,
    pub uplink: LightTimeLeg,
    pub downlink: LightTimeLeg,
}

impl TwoWaySolution {
    pub fn total_delay(&self) -> Duration {
        self.uplink.light_time + self.downlink.light_time
    }

    pub fn range_km(&self) -> f64 {
        (self.uplink.range_km() + self.downlink.range_km()) / 2.0
    }
}

/// Computes the exact two-way light-time solution using Picard fixed-point iteration.
/// IMPORTANT This does not check for elevation constraints. You must manually check that
/// at the t3_bounce epoch.
pub(crate) fn solve_two_way_picard(
    t3: Epoch,
    station: &GroundStation,
    traj: &Traj<Spacecraft>,
    almanac: &Almanac,
) -> Result<TwoWaySolution, ODError> {
    // Step 0: Anchor the reception state of the ground station in ICRF
    let gs_rx_orbit = station.to_orbit(t3, almanac).context(ODAlmanacSnafu {
        action: "building ground station orbit at t3",
    })?;
    let gs_rx_icrf = almanac
        .transform_to(gs_rx_orbit, ICRF, None)
        .context(ODAlmanacSnafu {
            action: "transforming station at t3 to ICRF",
        })?;
    let r3_icrf = gs_rx_icrf.radius_km;

    // Solve Downlink Leg (t3 -> t2)
    // Find t2 such that c * (t3 - t2) = || r_sc(t2) - r_gs(t3) ||
    let mut tau_down = Duration::ZERO;
    let mut r2_icrf = Vector3::zeros();
    let mut t2 = t3;

    // Exactly 3 iterations converge to sub-millimeter precision in ICRF
    for _ in 0..3 {
        t2 = t3 - tau_down;
        let sc_state = traj.at(t2).context(ODTrajSnafu {
            details: format!("interpolating spacecraft state at bounce epoch {t2}"),
        })?;

        // Transform spacecraft orbit to ICRF
        let sc_icrf = almanac
            .transform_to(sc_state.orbit, ICRF, None)
            .context(ODAlmanacSnafu {
                action: "transforming spacecraft at t2 to ICRF",
            })?;

        r2_icrf = sc_icrf.radius_km;
        let dist_down = (r2_icrf - r3_icrf).norm();
        tau_down = (dist_down / SPEED_OF_LIGHT_KM_S).seconds();
    }

    let downlink = LightTimeLeg {
        t_start: t2,
        t_end: t3,
        r_start_icrf: r2_icrf,
        r_end_icrf: r3_icrf,
        light_time: tau_down,
    };

    // Solve Uplink Leg (t2 -> t1)
    // Spacecraft state r_sc(t2) is now fixed.
    // Find t1 such that c * (t2 - t1) = || r_sc(t2) - r_gs(t1) ||
    let mut tau_up = tau_down; // Good initial guess
    let mut r1_icrf = Vector3::zeros();
    let mut t1 = t2 - tau_up;

    for _ in 0..3 {
        t1 = t2 - tau_up;
        let gs_tx_orbit = station.to_orbit(t1, almanac).context(ODAlmanacSnafu {
            action: "building ground station orbit at t1",
        })?;

        let gs_tx_icrf = almanac
            .transform_to(gs_tx_orbit, ICRF, None)
            .context(ODAlmanacSnafu {
                action: "transforming station at t1 to ICRF",
            })?;

        r1_icrf = gs_tx_icrf.radius_km;
        let dist_up = (r2_icrf - r1_icrf).norm();
        tau_up = (dist_up / SPEED_OF_LIGHT_KM_S).seconds();
    }

    let uplink = LightTimeLeg {
        t_start: t1,
        t_end: t2,
        r_start_icrf: r1_icrf,
        r_end_icrf: r2_icrf,
        light_time: tau_up,
    };

    // Package Two-Way Observables
    Ok(TwoWaySolution {
        t1_transmit: t1,
        t2_bounce: t2,
        t3_receive: t3,
        uplink,
        downlink,
    })
}
