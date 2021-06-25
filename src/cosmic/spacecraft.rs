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

use super::{Orbit, State, StmKind, TimeTagged};
use crate::dimensions::{Const, Matrix6, OMatrix, OVector};
use crate::dynamics::thrustctrl::Thruster;
use crate::errors::NyxError;
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
    /// Cr_partials only contains the column and row data needed to expand the 6x6 of the orbital STM (2x6+1 = 13), organized by row then column data
    pub cr_partials: OVector<f64, Const<13>>,
    /// Cd_partials only contains the column and row data needed to expand the 7x7 of the orbital STM + Cr partials (2x7+1 = 15), organized by row then column data
    pub cd_partials: OVector<f64, Const<15>>,
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
            cr_partials: OVector::<f64, Const<13>>::zeros(),
            cd_partials: OVector::<f64, Const<15>>::zeros(),
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
            ..Default::default()
        }
    }

    /// Initialize a spacecraft state from the SRP default 1.8 for coefficient of reflectivity (fuel mass and drag parameters nullified!)
    pub fn from_srp_defaults(orbit: Orbit, dry_mass_kg: f64, srp_area_m2: f64) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            srp_area_m2,
            ..Default::default()
        }
    }

    /// Initialize a spacecraft state from the SRP default 1.8 for coefficient of drag (fuel mass and SRP parameters nullified!)
    pub fn from_drag_defaults(orbit: Orbit, dry_mass_kg: f64, drag_area_m2: f64) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            drag_area_m2,
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

    pub fn with_orbit(self, orbit: Orbit) -> Self {
        let mut me = self;
        me.orbit = orbit;
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
            "{}\t{} kg",
            format!("{:.*x}", decimals, self.orbit),
            format!("{:.*}", decimals, self.dry_mass_kg + self.fuel_mass_kg),
        )
    }
}

impl fmt::LowerExp for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "{}\t{} kg",
            format!("{:.*e}", decimals, self.orbit),
            format!("{:.*e}", decimals, self.dry_mass_kg + self.fuel_mass_kg),
        )
    }
}

impl fmt::UpperExp for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "{}\t{} kg",
            format!("{:.*E}", decimals, self.orbit),
            format!("{:.*E}", decimals, self.dry_mass_kg + self.fuel_mass_kg),
        )
    }
}

impl fmt::LowerHex for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "{}\t{} kg",
            format!("{:.*x}", decimals, self.orbit),
            format!("{:.*}", decimals, self.dry_mass_kg + self.fuel_mass_kg),
        )
    }
}

impl fmt::UpperHex for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let decimals = f.precision().unwrap_or(6);
        write!(
            f,
            "{}\t{} kg",
            format!("{:.*X}", decimals, self.orbit),
            format!("{:.*e}", decimals, self.dry_mass_kg + self.fuel_mass_kg),
        )
    }
}

impl TimeTagged for Spacecraft {
    fn epoch(&self) -> Epoch {
        self.orbit.dt
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.orbit.dt = epoch
    }
}

impl State for Spacecraft {
    type Size = Const<8>;
    type VecLength = Const<73>;

    fn zeros() -> Self {
        Self::default()
    }

    /// The vector is organized as such:
    /// [X, Y, Z, Vx, Vy, Vz, Cr, Cd, Orbit_STM(36), Cr_partials (13), Cd_partials (15), Fuel mass ]
    fn as_vector(&self) -> Result<OVector<f64, Const<73>>, NyxError> {
        let orb_vec: OVector<f64, Const<42>> = self.orbit.as_vector()?;
        let mut vector = OVector::<f64, Const<73>>::zeros();
        for (i, val) in orb_vec.iter().enumerate() {
            // Place the orbit state first, then skip two (Cr, Cd), then copy orbit STM
            vector[if i < 6 { i } else { i + 2 }] = *val;
        }
        vector[6] = self.cr;
        vector[7] = self.cd;
        for (i, val) in self.cr_partials.iter().enumerate() {
            vector[i + 44] = *val;
        }
        for (i, val) in self.cd_partials.iter().enumerate() {
            vector[i + 57] = *val;
        }
        vector[72] = self.fuel_mass_kg;
        Ok(vector)
    }

    /// Vector is expected to be organized as such:
    /// [X, Y, Z, Vx, Vy, Vz, Cr, Cd, Orbit_STM(36), Cr_partials (13), Cd_partials (15), Fuel mass ]
    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Const<73>>) -> Result<(), NyxError> {
        self.set_epoch(epoch);
        // Rebuild the vectors
        let mut orbit_vec = OVector::<f64, Const<42>>::zeros();
        // Grab 44 elements and ignore the Cr and Cd
        for (i, val) in vector.fixed_rows::<44>(0).iter().enumerate() {
            if i == 6 || i == 7 {
                // Skip Cr and Cd
                continue;
            }
            orbit_vec[if i < 6 { i } else { i - 2 }] = *val;
        }
        self.orbit.set(epoch, &orbit_vec)?;
        self.cr = vector[6];
        self.cd = vector[7];
        // TODO: Invert the STM here based on STM Kind of orbit!!
        self.cr_partials = vector.fixed_rows::<13>(44).into_owned();
        self.cd_partials = vector.fixed_rows::<15>(57).into_owned();
        self.fuel_mass_kg = vector[72];
        Ok(())
    }

    /// diag(STM) = [X,Y,Z,Vx,Vy,Vz,Cr,Cd,Fuel]
    /// WARNING: Currently the STM assumes that the fuel mass is constant at ALL TIMES!
    fn stm(&self) -> Result<OMatrix<f64, Const<8>, Const<8>>, NyxError> {
        match self.orbit.stm {
            Some(stm) => {
                let mut rtn = OMatrix::<f64, Const<8>, Const<8>>::zeros();
                for i in 0..6 {
                    for j in 0..6 {
                        rtn[(i, j)] = stm[(i, j)];
                    }
                }
                for (i, val) in self.cr_partials.iter().enumerate() {
                    if i <= 7 {
                        // Row data
                        rtn[(6, i)] = *val;
                    } else {
                        rtn[(6 - i - 1, 6)] = *val;
                    }
                }
                for (i, val) in self.cd_partials.iter().enumerate() {
                    if i <= 8 {
                        // Row data
                        rtn[(7, i)] = *val;
                    } else {
                        rtn[(7 - i - 1, 7)] = *val;
                    }
                }
                // Fuel partial is zero for now, hence no rtn[(8,8)] = 0.0
                Ok(rtn)
            }
            None => Err(NyxError::StateTransitionMatrixUnset),
        }
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

impl Add<OVector<f64, Const<8>>> for Spacecraft {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(self, other: OVector<f64, Const<8>>) -> Self {
        let mut me = self;
        me.orbit.x += other[0];
        me.orbit.y += other[1];
        me.orbit.z += other[2];
        me.orbit.vx += other[3];
        me.orbit.vy += other[4];
        me.orbit.vz += other[5];
        me.cr += other[6];
        me.cd += other[7];

        me
    }
}
