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
use anise::constants::frames::{ICRF, SUN_J2000};
use anise::errors::AlmanacPhysicsSnafu;
use anise::prelude::Almanac;
use hifitime::{Duration, Epoch, TimeUnits};
use nalgebra::Vector3;

use crate::Spacecraft;
use crate::md::prelude::Traj;
use crate::od::ground_station::GroundStation;
use crate::od::{ODAlmanacSnafu, ODError, ODPlanetaryDataSnafu, ODTrajSnafu};
use snafu::ResultExt;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct LightTimeLeg {
    pub start_epoch: Epoch,
    pub end_epoch: Epoch,
    pub r_start_icrf_km: Vector3<f64>,
    pub r_end_icrf_km: Vector3<f64>,
    pub light_time: Duration,
}

impl LightTimeLeg {
    pub fn range_km(&self) -> f64 {
        (self.r_start_icrf_km - self.r_end_icrf_km).norm()
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
        0.5 * (self.t3_receive - self.t1_transmit).to_seconds() * SPEED_OF_LIGHT_KM_S
        // (self.uplink.range_km() + self.downlink.range_km()) / 2.0
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
    let sun_mu_km3_s2 = if station.relativistic_corrections == Some(true) {
        Some(
            almanac
                .frame_info(SUN_J2000)
                .context(ODPlanetaryDataSnafu {
                    action: "fetching Sun grav param for Shapiro delay",
                })?
                .mu_km3_s2()
                .context(AlmanacPhysicsSnafu {
                    action: "Sun mu not defined",
                })
                .context(ODAlmanacSnafu {
                    action: "fetching Sun grav param for Shapiro delay",
                })?,
        )
    } else {
        None
    };
    // Step 0: Anchor the reception state of the ground station in ICRF
    let gs_rx_orbit = station.to_orbit(t3, almanac).context(ODAlmanacSnafu {
        action: "building ground station orbit at t3",
    })?;
    let gs_rx_icrf = almanac
        .transform_to(gs_rx_orbit, ICRF, None)
        .context(ODAlmanacSnafu {
            action: "transforming station at t3 to ICRF",
        })?;
    let r3_icrf_km = gs_rx_icrf.radius_km;

    // Solve Downlink Leg (t3 -> t2)
    // Find t2 such that c * (t3 - t2) = || r_sc(t2) - r_gs(t3) ||
    let mut tau_down = Duration::ZERO;
    let mut r2_icrf_km = Vector3::zeros();
    let mut t2 = t3;

    // Exactly 3 iterations converge to sub-millimeter precision in ICRF
    for _ in 0..3 {
        t2 = t3 - tau_down;
        // NOTE Using with_context to lazy eval the format on the error
        let sc_state = traj.at(t2).with_context(|_| ODTrajSnafu {
            details: format!("interpolating spacecraft state at bounce epoch {t2}"),
        })?;

        // Transform spacecraft orbit to ICRF
        let sc_icrf = almanac
            .transform_to(sc_state.orbit, ICRF, None)
            .context(ODAlmanacSnafu {
                action: "transforming spacecraft at t2 to ICRF",
            })?;

        r2_icrf_km = sc_icrf.radius_km;
        let dist_down_km = (r2_icrf_km - r3_icrf_km).norm();
        // Geometric transit time + Solar Shapiro time dilation
        let tau_geometric = dist_down_km / SPEED_OF_LIGHT_KM_S;
        let tau_shapiro = if let Some(sun_mu) = sun_mu_km3_s2 {
            shapiro_delay_s(&r2_icrf_km, &r3_icrf_km, dist_down_km, sun_mu)
        } else {
            0.0
        };
        tau_down = (tau_geometric + tau_shapiro).seconds();
    }

    let downlink = LightTimeLeg {
        start_epoch: t2,
        end_epoch: t3,
        r_start_icrf_km: r2_icrf_km,
        r_end_icrf_km: r3_icrf_km,
        light_time: tau_down,
    };

    // Solve Uplink Leg (t2 -> t1)
    // Spacecraft state r_sc(t2) is now fixed.
    // Find t1 such that c * (t2 - t1) = || r_sc(t2) - r_gs(t1) ||
    let mut tau_up = tau_down; // Good initial guess
    let mut r1_icrf_km = Vector3::zeros();
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

        r1_icrf_km = gs_tx_icrf.radius_km;
        let dist_up_km = (r2_icrf_km - r1_icrf_km).norm();
        let tau_geometric = dist_up_km / SPEED_OF_LIGHT_KM_S;
        let tau_shapiro = if let Some(sun_mu) = sun_mu_km3_s2 {
            shapiro_delay_s(&r1_icrf_km, &r2_icrf_km, dist_up_km, sun_mu)
        } else {
            0.0
        };
        tau_up = (tau_geometric + tau_shapiro).seconds();
    }

    let uplink = LightTimeLeg {
        start_epoch: t1,
        end_epoch: t2,
        r_start_icrf_km: r1_icrf_km,
        r_end_icrf_km: r2_icrf_km,
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

/// Shapiro delay is computed when relativistic corrections are enabled.
/// Signals passing near a massive object take slightly longer to travel to a target and longer to return than they
/// would if the mass of the object were not present. The time delay is caused by time dilation, which increases
/// the time it takes light to travel a given distance from the perspective of an outside observer.
fn shapiro_delay_s(
    r_start_icrf_km: &Vector3<f64>,
    r_end_icrf_km: &Vector3<f64>,
    rho_km: f64,
    sun_mu_km3_s2: f64,
) -> f64 {
    let num = r_start_icrf_km.norm() + r_end_icrf_km.norm() + rho_km;
    let denom = r_start_icrf_km.norm() + r_end_icrf_km.norm() - rho_km;

    if denom <= 0.0 || num <= 0.0 {
        return 0.0;
    }

    // PPN formulation (gamma = 1): 2 * mu / c^3 * ln(...)
    let factor = (2.0 * sun_mu_km3_s2) / SPEED_OF_LIGHT_KM_S.powi(3);
    factor * (num / denom).ln()
}
