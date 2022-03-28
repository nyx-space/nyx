/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use nalgebra::Vector3;

use super::{Orbit, State};
use crate::dynamics::guidance::Thruster;
use crate::errors::NyxError;
use crate::linalg::{Const, DimName, Matrix6, OMatrix, OVector, Vector6};
use crate::mc::MultivariateNormal;
use crate::md::StateParameter;
use crate::time::Epoch;
use crate::utils::rss_orbit_errors;
use std::default::Default;
use std::fmt;
use std::ops::Add;

/// Defines a spacecraft extension.
/// This is useful for highly specialized guidance laws that need to store additional data in the spacecraft state.
/// Most guidance laws can be implemented directly with the `Spacecraft` structure.
pub trait SpacecraftExt: Clone + Copy + Default + fmt::Debug + Send + Sync {}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GuidanceMode {
    /// Guidance is turned off and Guidance Law may switch mode to Thrust for next call
    Coast,
    /// Guidance is turned on and Guidance Law may switch mode to Coast for next call
    Thrust,
    /// Guidance is turned off and Guidance Law may not change its mode (will need to be done externally to the guidance law).
    Inhibit,
    Custom(u8),
}

impl Default for GuidanceMode {
    fn default() -> Self {
        Self::Coast
    }
}

impl SpacecraftExt for GuidanceMode {}

/// A spacecraft state
#[derive(Clone, Copy, Debug)]
pub struct BaseSpacecraft<X: SpacecraftExt> {
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
    /// Optionally stores the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)
    pub stm: Option<OMatrix<f64, Const<9>, Const<9>>>,
    /// Any extra information or extension that is needed for specific guidance laws
    pub ext: X,
}

/// A spacecraft state
pub type Spacecraft = BaseSpacecraft<GuidanceMode>;

impl<X: SpacecraftExt> Default for BaseSpacecraft<X> {
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
            stm: None,
            ext: X::default(),
        }
    }
}

impl<X: SpacecraftExt> BaseSpacecraft<X> {
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

    /// Initialize a spacecraft state from only a thruster and mass. Use this when designing guidance laws whilke ignoring drag and SRP.
    pub fn from_thruster(
        orbit: Orbit,
        dry_mass_kg: f64,
        fuel_mass_kg: f64,
        thruster: Thruster,
        ext: X,
    ) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            fuel_mass_kg,
            thruster: Some(thruster),
            ext,
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

    pub fn with_dv(self, dv: Vector3<f64>) -> Self {
        let mut me = self;
        me.orbit.apply_dv(dv);
        me
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

    /// Allows the automatic conversion of a BaseSpacecraft into a standard spacecraft.
    /// WARNING: This will IGNORE the `extra` of BaseSpacecraft
    pub fn degrade_to_spacecraft(&self) -> Spacecraft {
        Spacecraft {
            orbit: self.orbit,
            dry_mass_kg: self.dry_mass_kg,
            fuel_mass_kg: self.fuel_mass_kg,
            srp_area_m2: self.srp_area_m2,
            drag_area_m2: self.drag_area_m2,
            cr: self.cr,
            cd: self.cd,
            thruster: self.thruster,
            stm: self.stm,
            ext: GuidanceMode::Coast,
        }
    }
}

impl<X: SpacecraftExt> PartialEq for BaseSpacecraft<X> {
    fn eq(&self, other: &Self) -> bool {
        let mass_tol = 1e-6; // milligram
        self.orbit == other.orbit
            && (self.dry_mass_kg - other.dry_mass_kg).abs() < mass_tol
            && (self.fuel_mass_kg - other.fuel_mass_kg).abs() < mass_tol
            && (self.cr - other.cr).abs() < std::f64::EPSILON
            && (self.cd - other.cd).abs() < std::f64::EPSILON
    }
}

impl<X: SpacecraftExt> fmt::Display for BaseSpacecraft<X> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*}", mass_prec, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*}", orbit_prec, self.orbit),
            self.ext,
        )
    }
}

impl<X: SpacecraftExt> fmt::LowerExp for BaseSpacecraft<X> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*e}", mass_prec, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*e}", orbit_prec, self.orbit),
            self.ext,
        )
    }
}

impl<X: SpacecraftExt> fmt::UpperExp for BaseSpacecraft<X> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*E}", mass_prec, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*E}", orbit_prec, self.orbit),
            self.ext,
        )
    }
}

impl<X: SpacecraftExt> fmt::LowerHex for BaseSpacecraft<X> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*}", mass_prec, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*x}", orbit_prec, self.orbit),
            self.ext,
        )
    }
}

impl<X: SpacecraftExt> fmt::UpperHex for BaseSpacecraft<X> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*e}", mass_prec, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*X}", orbit_prec, self.orbit),
            self.ext,
        )
    }
}

impl<X: SpacecraftExt> State for BaseSpacecraft<X> {
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

    fn value(&self, param: &StateParameter) -> Result<f64, NyxError> {
        match *param {
            StateParameter::Cd => Ok(self.cd),
            StateParameter::Cr => Ok(self.cr),
            StateParameter::FuelMass => Ok(self.fuel_mass_kg),
            StateParameter::Isp => match self.thruster {
                Some(thruster) => Ok(thruster.isp_s),
                None => Err(NyxError::NoThrusterAvail),
            },
            StateParameter::Thrust => match self.thruster {
                Some(thruster) => Ok(thruster.thrust_N),
                None => Err(NyxError::NoThrusterAvail),
            },
            _ => self.orbit.value(param),
        }
    }

    fn set_value(&mut self, param: &StateParameter, val: f64) -> Result<(), NyxError> {
        match *param {
            StateParameter::Cd => self.cd = val,
            StateParameter::Cr => self.cr = val,
            StateParameter::FuelMass => self.fuel_mass_kg = val,
            StateParameter::Isp => match self.thruster {
                Some(ref mut thruster) => thruster.isp_s = val,
                None => return Err(NyxError::NoThrusterAvail),
            },
            StateParameter::Thrust => match self.thruster {
                Some(ref mut thruster) => thruster.thrust_N = val,
                None => return Err(NyxError::NoThrusterAvail),
            },
            _ => return self.orbit.set_value(param, val),
        }
        Ok(())
    }
}

impl<X: SpacecraftExt> Add<OVector<f64, Const<6>>> for BaseSpacecraft<X> {
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

impl<X: SpacecraftExt> Add<OVector<f64, Const<9>>> for BaseSpacecraft<X> {
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

impl Spacecraft {
    /// Returns a copy of the state with the provided guidance mode
    pub fn with_guidance_mode(self, mode: GuidanceMode) -> Self {
        let mut me = self;
        me.ext = mode;
        me
    }

    pub fn mode(&self) -> GuidanceMode {
        self.ext
    }

    pub fn mut_mode(&mut self, mode: GuidanceMode) {
        self.ext = mode;
    }
}
