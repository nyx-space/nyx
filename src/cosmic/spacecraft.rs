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

use super::{Orbit, State};
use crate::dynamics::guidance::Thruster;
use crate::errors::NyxError;
use crate::linalg::{Const, DimName, Matrix6, OMatrix, OVector};
use crate::time::Epoch;
use crate::utils::rss_orbit_errors;
use std::default::Default;
use std::fmt;
use std::ops::Add;

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
    /// Optionally stores the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)
    pub stm: Option<OMatrix<f64, Const<9>, Const<9>>>,
}

impl Default for Spacecraft {
    fn default() -> Self {
        Self {
            orbit: Orbit::zeros(),
            dry_mass_kg: 0.0,
            fuel_mass_kg: 0.0,
            srp_area_m2: 0.0,
            drag_area_m2: 0.0,
            cr: 1.8,
            cd: 2.2,
            thruster: None,
            mode: GuidanceMode::Coast,
            stm: None,
        }
    }
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
            stm: orbit
                .stm
                .map(|_| OMatrix::<f64, Const<9>, Const<9>>::identity()),
            ..Default::default()
        }
    }

    /// Initialize a spacecraft state from the SRP default 1.8 for coefficient of reflectivity (fuel mass and drag parameters nullified!)
    pub fn from_srp_defaults(orbit: Orbit, dry_mass_kg: f64, srp_area_m2: f64) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            srp_area_m2,
            stm: orbit
                .stm
                .map(|_| OMatrix::<f64, Const<9>, Const<9>>::identity()),
            ..Default::default()
        }
    }

    /// Initialize a spacecraft state from the SRP default 1.8 for coefficient of drag (fuel mass and SRP parameters nullified!)
    pub fn from_drag_defaults(orbit: Orbit, dry_mass_kg: f64, drag_area_m2: f64) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            drag_area_m2,
            stm: orbit
                .stm
                .map(|_| OMatrix::<f64, Const<9>, Const<9>>::identity()),
            ..Default::default()
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
            thruster: Some(thruster),
            mode: init_mode,
            stm: orbit
                .stm
                .map(|_| OMatrix::<f64, Const<9>, Const<9>>::identity()),
            ..Default::default()
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

    /// Returns a copy of the state with a new orbit
    pub fn with_orbit(self, orbit: Orbit) -> Self {
        let mut me = self;
        me.orbit = orbit;
        me.stm = orbit
            .stm
            .map(|_| OMatrix::<f64, Const<9>, Const<9>>::identity());
        me
    }

    /// Returns a copy of the state with the provided guidance mode
    pub fn with_guidance_mode(self, mode: GuidanceMode) -> Self {
        let mut me = self;
        me.mode = mode;
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
        self.stm = Some(OMatrix::<f64, Const<9>, Const<9>>::identity());
    }

    /// Copies the current state but sets the STM to identity
    pub fn with_stm(self) -> Self {
        let mut me = self;
        me.enable_stm();
        me
    }

    /// Returns the total mass in kilograms
    pub fn mass_kg(&self) -> f64 {
        self.dry_mass_kg + self.fuel_mass_kg
    }
}

impl PartialEq for Spacecraft {
    fn eq(&self, other: &Spacecraft) -> bool {
        let mass_tol = 1e-6; // milligram
        self.orbit == other.orbit
            && (self.dry_mass_kg - other.dry_mass_kg).abs() < mass_tol
            && (self.fuel_mass_kg - other.fuel_mass_kg).abs() < mass_tol
            && (self.cr - other.cr).abs() < std::f64::EPSILON
            && (self.cd - other.cd).abs() < std::f64::EPSILON
    }
}

impl fmt::Display for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {}kg  @  {}",
            format!("{:.*}", decimals, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*}", decimals, self.orbit),
        )
    }
}

impl fmt::LowerExp for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {}kg  @  {}",
            format!("{:.*e}", decimals, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*e}", decimals, self.orbit),
        )
    }
}

impl fmt::UpperExp for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {}kg  @  {}",
            format!("{:.*E}", decimals, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*E}", decimals, self.orbit),
        )
    }
}

impl fmt::LowerHex for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {}kg  @  {}",
            format!("{:.*}", decimals, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*x}", decimals, self.orbit),
        )
    }
}

impl fmt::UpperHex for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {}kg  @  {}",
            format!("{:.*e}", decimals, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*X}", decimals, self.orbit),
        )
    }
}

impl State for Spacecraft {
    type Size = Const<9>;
    type VecLength = Const<90>;

    fn reset_stm(&mut self) {
        self.orbit.reset_stm();
        self.stm = Some(OMatrix::<f64, Const<9>, Const<9>>::identity());
    }

    fn zeros() -> Self {
        Self::default()
    }

    /// The vector is organized as such:
    /// [X, Y, Z, Vx, Vy, Vz, Cr, Cd, Fuel mass, STM(9x9)]
    fn as_vector(&self) -> Result<OVector<f64, Const<90>>, NyxError> {
        let mut vector = OVector::<f64, Const<90>>::zeros();
        // Set the orbit state info
        for (i, val) in self.orbit.to_cartesian_vec().iter().enumerate() {
            // Place the orbit state first, then skip three (Cr, Cd, Fuel), then copy orbit STM
            vector[if i < 6 { i } else { i + 3 }] = *val;
        }
        // Set the spacecraft parameters
        vector[6] = self.cr;
        vector[7] = self.cd;
        vector[8] = self.fuel_mass_kg;
        if let Some(mut stm) = self.stm {
            // Set the 6x6 of the orbit STM first
            for i in 0..6 {
                for j in 0..6 {
                    stm[(i, j)] = self.orbit.stm().unwrap()[(i, j)];
                }
            }
            for (idx, stm_val) in stm.as_slice().iter().enumerate() {
                vector[idx + Self::Size::dim()] = *stm_val;
            }
        }
        Ok(vector)
    }

    /// Vector is expected to be organized as such:
    /// [X, Y, Z, Vx, Vy, Vz, Cr, Cd, Fuel mass, STM(9x9)]
    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Const<90>>) -> Result<(), NyxError> {
        self.set_epoch(epoch);
        let sc_state =
            OVector::<f64, Self::Size>::from_column_slice(&vector.as_slice()[..Self::Size::dim()]);
        let sc_full_stm = OMatrix::<f64, Self::Size, Self::Size>::from_column_slice(
            &vector.as_slice()[Self::Size::dim()..],
        );

        if self.stm.is_some() {
            self.stm = Some(sc_full_stm);
        }

        // Extract the orbit information
        let orbit_state = sc_state.fixed_rows::<6>(0).into_owned();
        let orbit_stm = sc_full_stm.fixed_slice::<6, 6>(0, 0).into_owned();
        // Rebuild the orbit vector
        let mut orbit_vec = OVector::<f64, Const<42>>::zeros();
        orbit_vec[0] = orbit_state[0];
        orbit_vec[1] = orbit_state[1];
        orbit_vec[2] = orbit_state[2];
        orbit_vec[3] = orbit_state[3];
        orbit_vec[4] = orbit_state[4];
        orbit_vec[5] = orbit_state[5];
        for (idx, stm_val) in orbit_stm.as_slice().iter().enumerate() {
            orbit_vec[idx + 6] = *stm_val;
        }
        // And set the orbit information
        self.orbit.set(epoch, &orbit_vec)?;
        self.cr = sc_state[6];
        self.cd = sc_state[7];
        self.fuel_mass_kg = sc_state[8];
        Ok(())
    }

    /// diag(STM) = [X,Y,Z,Vx,Vy,Vz,Cr,Cd,Fuel]
    /// WARNING: Currently the STM assumes that the fuel mass is constant at ALL TIMES!
    fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, NyxError> {
        match self.stm {
            Some(stm) => Ok(stm),
            None => Err(NyxError::StateTransitionMatrixUnset),
        }
    }

    fn epoch(&self) -> Epoch {
        self.orbit.dt
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.orbit.dt = epoch
    }

    fn add(self, other: OVector<f64, Self::Size>) -> Self {
        self + other
    }
}

impl Add<OVector<f64, Const<6>>> for Spacecraft {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(self, other: OVector<f64, Const<6>>) -> Self {
        let mut me = self;
        me.orbit.x += other[0];
        me.orbit.y += other[1];
        me.orbit.z += other[2];
        me.orbit.vx += other[3];
        me.orbit.vy += other[4];
        me.orbit.vz += other[5];

        me
    }
}

impl Add<OVector<f64, Const<9>>> for Spacecraft {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(self, other: OVector<f64, Const<9>>) -> Self {
        let mut me = self;
        me.orbit.x += other[0];
        me.orbit.y += other[1];
        me.orbit.z += other[2];
        me.orbit.vx += other[3];
        me.orbit.vy += other[4];
        me.orbit.vz += other[5];
        me.cr += other[6];
        me.cd += other[7];
        me.fuel_mass_kg += other[8];

        me
    }
}
