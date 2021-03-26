/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::{Orbit, StmKind};
use crate::dimensions::Matrix6;
use crate::dynamics::thrustctrl::Thruster;
use crate::utils::rss_orbit_errors;
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GuidanceMode {
    Coast,
    Thrust,
    Custom(u8),
}

/// A spacecraft state
#[derive(Clone, Copy, Debug)]
pub struct Spacecraft {
    /// Initial orbit the vehicle is in
    pub orbit: Orbit,
    /// Dry mass, i.e. mass without fuel, in kg
    pub dry_mass_kg: f64,
    /// Fuel mass (if fuel mass is negative, thrusting will fail, unless configured to break laws of physics)
    pub fuel_mass_kg: f64,
    /// in m^2
    pub srp_area_m2: f64,
    /// in m^2
    pub drag_area_m2: f64,
    /// coefficient of reflectivity, must be between 0.0 (translucent) and 2.0 (all radiation absorbed and twice the force is transmitted back).
    pub cr: f64,
    /// coefficient of drag; (spheres are between 2.0 and 2.1, use 2.2 in Earth's atmosphere).
    pub cd: f64,
    pub thruster: Option<Thruster>,
    /// Guidance mode determines whether the thruster should fire or not
    pub mode: GuidanceMode,
}

impl Spacecraft {
    /// Initialize a spacecraft state from all of its parameters
    pub fn new(
        orbit: Orbit,
        dry_mass_kg: f64,
        fuel_mass_kg: f64,
        srp_area_m2: f64,
        drag_area_m2: f64,
        cr: f64,
        cd: f64,
    ) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            fuel_mass_kg,
            srp_area_m2,
            drag_area_m2,
            cr,
            cd,
            thruster: None,
            mode: GuidanceMode::Coast,
        }
    }

    /// Initialize a spacecraft state from the SRP default 1.8 for coefficient of reflectivity (fuel mass and drag parameters nullified!)
    pub fn from_srp_defaults(orbit: Orbit, dry_mass_kg: f64, srp_area_m2: f64) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            fuel_mass_kg: 0.0,
            srp_area_m2,
            drag_area_m2: 0.0,
            cr: 1.8,
            cd: 0.0,
            thruster: None,
            mode: GuidanceMode::Coast,
        }
    }

    /// Initialize a spacecraft state from the SRP default 1.8 for coefficient of drag (fuel mass and SRP parameters nullified!)
    pub fn from_drag_defaults(orbit: Orbit, dry_mass_kg: f64, drag_area_m2: f64) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            fuel_mass_kg: 0.0,
            srp_area_m2: 0.0,
            drag_area_m2,
            cr: 0.0,
            cd: 2.2,
            thruster: None,
            mode: GuidanceMode::Coast,
        }
    }

    /// Initialize a spacecraft state from only a thruster and mass. Use this when designing control laws whilke ignoring drag and SRP.
    pub fn from_thruster(
        orbit: Orbit,
        dry_mass_kg: f64,
        fuel_mass_kg: f64,
        thruster: Thruster,
        init_mode: GuidanceMode,
    ) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            fuel_mass_kg,
            srp_area_m2: 0.0,
            drag_area_m2: 0.0,
            cr: 0.0,
            cd: 0.0,
            thruster: Some(thruster),
            mode: init_mode,
        }
    }

    /// Returns a copy of the state with a new dry mass
    pub fn with_dry_mass(self, dry_mass_kg: f64) -> Self {
        let mut me = self;
        me.dry_mass_kg = dry_mass_kg;
        me
    }

    /// Returns a copy of the state with a new fuel mass
    pub fn with_fuel_mass(self, fuel_mass_kg: f64) -> Self {
        let mut me = self;
        me.fuel_mass_kg = fuel_mass_kg;
        me
    }

    /// Returns a copy of the state with a new SRP area and CR
    pub fn with_srp(self, srp_area_m2: f64, cr: f64) -> Self {
        let mut me = self;
        me.srp_area_m2 = srp_area_m2;
        me.cr = cr;
        me
    }

    /// Returns a copy of the state with a new SRP area
    pub fn with_srp_area(self, srp_area_m2: f64) -> Self {
        let mut me = self;
        me.srp_area_m2 = srp_area_m2;
        me
    }

    /// Returns a copy of the state with a new coefficient of reflectivity
    pub fn with_cr(self, cr: f64) -> Self {
        let mut me = self;
        me.cr = cr;
        me
    }

    /// Returns a copy of the state with a new drag area and CD
    pub fn with_drag(self, drag_area_m2: f64, cd: f64) -> Self {
        let mut me = self;
        me.drag_area_m2 = drag_area_m2;
        me.cd = cd;
        me
    }

    /// Returns a copy of the state with a new SRP area
    pub fn with_drag_area(self, drag_area_m2: f64) -> Self {
        let mut me = self;
        me.drag_area_m2 = drag_area_m2;
        me
    }

    /// Returns a copy of the state with a new coefficient of drag
    pub fn with_cd(self, cd: f64) -> Self {
        let mut me = self;
        me.cd = cd;
        me
    }

    /// Returns the root sum square error between this spacecraft and the other, in kilometers for the position, kilometers per second in velocity, and kilograms in fuel
    pub fn rss(&self, other: &Self) -> (f64, f64, f64) {
        let (p, v) = rss_orbit_errors(&self.orbit, &other.orbit);
        (
            p,
            v,
            (self.fuel_mass_kg - other.fuel_mass_kg).powi(2).sqrt(),
        )
    }

    /// Sets the STM of this state of identity, which also enables computation of the STM for spacecraft navigation
    pub fn enable_stm(&mut self) {
        self.orbit.stm = Some(Matrix6::identity());
        self.orbit.stm_kind = StmKind::Step;
    }

    /// Sets the STM of this state of identity, which also enables computation of the STM for trajectory optimization
    pub fn enable_traj_stm(&mut self) {
        self.orbit.stm = Some(Matrix6::identity());
        self.orbit.stm_kind = StmKind::Traj;
    }

    /// Copies the current state but sets the STM to identity
    pub fn with_stm(self) -> Self {
        let mut me = self;
        me.enable_stm();
        me
    }

    /// Sets the STM of this state of identity
    pub fn stm_identity(&mut self) {
        self.orbit.stm = Some(Matrix6::identity());
    }

    /// Unwraps this STM, or panics if unset.
    pub fn stm(&self) -> Matrix6<f64> {
        self.orbit.stm.unwrap()
    }
}

impl PartialEq for Spacecraft {
    fn eq(&self, other: &Spacecraft) -> bool {
        let mass_tol = 1e-6; // milligram
        self.orbit == other.orbit
            && (self.dry_mass_kg - other.dry_mass_kg).abs() < mass_tol
            && (self.fuel_mass_kg - other.fuel_mass_kg).abs() < mass_tol
    }
}

impl fmt::Display for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:o}\t{} kg",
            self.orbit,
            self.dry_mass_kg + self.fuel_mass_kg
        )
    }
}

impl fmt::LowerExp for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:o}\t{:e} kg",
            self.orbit,
            self.dry_mass_kg + self.fuel_mass_kg
        )
    }
}
