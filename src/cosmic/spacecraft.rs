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

use nalgebra::Vector3;
use serde::{Deserialize, Serialize};

use super::eclipse::Cosm;
use super::{Orbit, State};
use crate::dynamics::guidance::Thruster;
use crate::errors::NyxError;
use crate::io::{orbit_from_str, ConfigRepr, Configurable};
use crate::linalg::{Const, DimName, Matrix6, OMatrix, OVector};
use crate::md::StateParameter;
use crate::time::Epoch;
use crate::utils::rss_orbit_errors;

use std::default::Default;
use std::fmt;
use std::ops::Add;
use std::sync::Arc;

use crate::io::ConfigError;
#[cfg(feature = "python")]
use pyo3::prelude::*;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize, Deserialize)]
#[cfg_attr(feature = "python", pyclass)]
pub enum GuidanceMode {
    /// Guidance is turned off and Guidance Law may switch mode to Thrust for next call
    Coast,
    /// Guidance is turned on and Guidance Law may switch mode to Coast for next call
    Thrust,
    /// Guidance is turned off and Guidance Law may not change its mode (will need to be done externally to the guidance law).
    Inhibit,
}

impl Default for GuidanceMode {
    fn default() -> Self {
        Self::Coast
    }
}

impl From<f64> for GuidanceMode {
    fn from(value: f64) -> Self {
        if value >= 1.0 {
            Self::Thrust
        } else if value < 0.0 {
            Self::Inhibit
        } else {
            Self::Coast
        }
    }
}

impl From<GuidanceMode> for f64 {
    fn from(mode: GuidanceMode) -> f64 {
        match mode {
            GuidanceMode::Coast => 0.0,
            GuidanceMode::Thrust => 1.0,
            GuidanceMode::Inhibit => -1.0,
        }
    }
}

/// A spacecraft state, composed of its orbit, its dry and fuel (wet) masses (in kg), its SRP configuration, its drag configuration, its thruster configuration, and its guidance mode.
///
/// Optionally, the spacecraft state can also store the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM).
#[derive(Clone, Copy, Debug, Serialize, Deserialize)]
#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.cosmic"))]
#[cfg_attr(
    feature = "python",
    pyo3(
        text_signature = "(orbit, dry_mass_kg, fuel_mass_kg, srp_area_m2, drag_area_m2, cr, cd, thruster, mode)"
    )
)]
pub struct Spacecraft {
    /// Initial orbit the vehicle is in
    #[serde(deserialize_with = "orbit_from_str")]
    pub orbit: Orbit,
    /// Dry mass, i.e. mass without fuel, in kg
    pub dry_mass_kg: f64,
    /// Fuel mass (if fuel mass is negative, thrusting will fail, unless configured to break laws of physics)
    pub fuel_mass_kg: f64,
    /// Solar Radiation Pressure configuration for this spacecraft
    #[serde(default)]
    pub srp: SrpConfig,
    #[serde(default)]
    pub drag: DragConfig,
    pub thruster: Option<Thruster>,
    /// Any extra information or extension that is needed for specific guidance laws
    #[serde(default)]
    pub mode: GuidanceMode,
    /// Optionally stores the state transition matrix from the start of the propagation until the current time (i.e. trajectory STM, not step-size STM)
    #[serde(skip)]
    pub stm: Option<OMatrix<f64, Const<9>, Const<9>>>,
}

impl Default for Spacecraft {
    fn default() -> Self {
        Self {
            orbit: Orbit::zeros(),
            dry_mass_kg: 0.0,
            fuel_mass_kg: 0.0,
            srp: SrpConfig::default(),
            drag: DragConfig::default(),
            thruster: None,
            mode: GuidanceMode::default(),
            stm: None,
        }
    }
}

#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(text_signature = "(area_m2, cr=1.8)"))]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.cosmic"))]
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
/// The Solar Radiation Pressure configuration for a spacecraft
pub struct SrpConfig {
    /// solar radiation pressure area
    pub area_m2: f64,
    /// coefficient of reflectivity, must be between 0.0 (translucent) and 2.0 (all radiation absorbed and twice the force is transmitted back), defaults to 1.8
    pub cr: f64,
}

impl SrpConfig {
    /// Initialize the SRP from the c_r default and the provided drag area
    pub fn from_area(area_m2: f64) -> Self {
        Self {
            area_m2,
            ..Default::default()
        }
    }
}

impl Default for SrpConfig {
    fn default() -> Self {
        Self {
            area_m2: 0.0,
            cr: 1.8,
        }
    }
}

#[cfg_attr(feature = "python", pyclass)]
#[cfg_attr(feature = "python", pyo3(text_signature = "(area_m2, cd=2.2)"))]
#[cfg_attr(feature = "python", pyo3(module = "nyx_space.cosmic"))]
#[derive(Copy, Clone, Debug, Serialize, Deserialize, PartialEq)]
/// The drag configuration for a spacecraft
pub struct DragConfig {
    /// drag area
    pub area_m2: f64,
    /// coefficient of drag; (spheres are between 2.0 and 2.1, use 2.2 in Earth's atmosphere (default)).
    pub cd: f64,
}

impl DragConfig {
    /// Initialize the drag from the c_d default and the provided drag area
    pub fn from_area(area_m2: f64) -> Self {
        Self {
            area_m2,
            ..Default::default()
        }
    }
}

impl Default for DragConfig {
    fn default() -> Self {
        Self {
            area_m2: 0.0,
            cd: 2.2,
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
            srp: SrpConfig {
                area_m2: srp_area_m2,
                cr,
            },
            drag: DragConfig {
                area_m2: drag_area_m2,
                cd,
            },
            stm: orbit
                .stm
                .map(|_| OMatrix::<f64, Const<9>, Const<9>>::identity()),
            ..Default::default()
        }
    }

    /// Initialize a spacecraft state from only a thruster and mass. Use this when designing guidance laws while ignoring drag and SRP.
    pub fn from_thruster(
        orbit: Orbit,
        dry_mass_kg: f64,
        fuel_mass_kg: f64,
        thruster: Thruster,
        mode: GuidanceMode,
    ) -> Self {
        Self {
            orbit,
            dry_mass_kg,
            fuel_mass_kg,
            thruster: Some(thruster),
            mode,
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
            srp: SrpConfig::from_area(srp_area_m2),
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
            drag: DragConfig::from_area(drag_area_m2),
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
        me.srp = SrpConfig {
            area_m2: srp_area_m2,
            cr,
        };

        me
    }

    /// Returns a copy of the state with a new SRP area
    pub fn with_srp_area(self, srp_area_m2: f64) -> Self {
        let mut me = self;
        me.srp.area_m2 = srp_area_m2;
        me
    }

    /// Returns a copy of the state with a new coefficient of reflectivity
    pub fn with_cr(self, cr: f64) -> Self {
        let mut me = self;
        me.srp.cr = cr;
        me
    }

    /// Returns a copy of the state with a new drag area and CD
    pub fn with_drag(self, drag_area_m2: f64, cd: f64) -> Self {
        let mut me = self;
        me.drag = DragConfig {
            area_m2: drag_area_m2,
            cd,
        };
        me
    }

    /// Returns a copy of the state with a new SRP area
    pub fn with_drag_area(self, drag_area_m2: f64) -> Self {
        let mut me = self;
        me.drag.area_m2 = drag_area_m2;
        me
    }

    /// Returns a copy of the state with a new coefficient of drag
    pub fn with_cd(self, cd: f64) -> Self {
        let mut me = self;
        me.drag.cd = cd;
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

    /// Returns a copy of the state with the provided guidance mode
    pub fn with_guidance_mode(self, mode: GuidanceMode) -> Self {
        let mut me = self;
        me.mode = mode;
        me
    }

    pub fn mode(&self) -> GuidanceMode {
        self.mode
    }

    pub fn mut_mode(&mut self, mode: GuidanceMode) {
        self.mode = mode;
    }
}

impl PartialEq for Spacecraft {
    fn eq(&self, other: &Self) -> bool {
        let mass_tol = 1e-6; // milligram
        self.orbit == other.orbit
            && (self.dry_mass_kg - other.dry_mass_kg).abs() < mass_tol
            && (self.fuel_mass_kg - other.fuel_mass_kg).abs() < mass_tol
            && self.srp == other.srp
            && self.drag == other.drag
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::Display for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*}", mass_prec, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*}", orbit_prec, self.orbit),
            self.mode,
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::LowerExp for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*e}", mass_prec, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*e}", orbit_prec, self.orbit),
            self.mode,
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::UpperExp for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*E}", mass_prec, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*E}", orbit_prec, self.orbit),
            self.mode,
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::LowerHex for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*}", mass_prec, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*x}", orbit_prec, self.orbit),
            self.mode,
        )
    }
}

#[allow(clippy::format_in_format_args)]
impl fmt::UpperHex for Spacecraft {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let mass_prec = f.precision().unwrap_or(3);
        let orbit_prec = f.precision().unwrap_or(6);
        write!(
            f,
            "total mass = {} kg @  {}  {:?}",
            format!("{:.*e}", mass_prec, self.dry_mass_kg + self.fuel_mass_kg),
            format!("{:.*X}", orbit_prec, self.orbit),
            self.mode,
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
        vector[6] = self.srp.cr;
        vector[7] = self.drag.cd;
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
        let orbit_stm = sc_full_stm.fixed_view::<6, 6>(0, 0).into_owned();
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
        self.srp.cr = sc_state[6];
        self.drag.cd = sc_state[7];
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
        self.orbit.epoch
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.orbit.epoch = epoch
    }

    fn add(self, other: OVector<f64, Self::Size>) -> Self {
        self + other
    }

    fn value(&self, param: StateParameter) -> Result<f64, NyxError> {
        match param {
            StateParameter::Cd => Ok(self.drag.cd),
            StateParameter::Cr => Ok(self.srp.cr),
            StateParameter::DryMass => Ok(self.dry_mass_kg),
            StateParameter::FuelMass => Ok(self.fuel_mass_kg),
            StateParameter::Isp => match self.thruster {
                Some(thruster) => Ok(thruster.isp_s),
                None => Err(NyxError::NoThrusterAvail),
            },
            StateParameter::Thrust => match self.thruster {
                Some(thruster) => Ok(thruster.thrust_N),
                None => Err(NyxError::NoThrusterAvail),
            },
            StateParameter::GuidanceMode => Ok(self.mode.into()),
            _ => self.orbit.value(param),
        }
    }

    fn set_value(&mut self, param: StateParameter, val: f64) -> Result<(), NyxError> {
        match param {
            StateParameter::Cd => self.drag.cd = val,
            StateParameter::Cr => self.srp.cr = val,
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

    fn unset_stm(&mut self) {
        self.stm = None;
    }
}

impl Add<OVector<f64, Const<6>>> for Spacecraft {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(self, other: OVector<f64, Const<6>>) -> Self {
        let mut me = self;
        me.orbit.x_km += other[0];
        me.orbit.y_km += other[1];
        me.orbit.z_km += other[2];
        me.orbit.vx_km_s += other[3];
        me.orbit.vy_km_s += other[4];
        me.orbit.vz_km_s += other[5];

        me
    }
}

impl Add<OVector<f64, Const<9>>> for Spacecraft {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(self, other: OVector<f64, Const<9>>) -> Self {
        let mut me = self;
        me.orbit.x_km += other[0];
        me.orbit.y_km += other[1];
        me.orbit.z_km += other[2];
        me.orbit.vx_km_s += other[3];
        me.orbit.vy_km_s += other[4];
        me.orbit.vz_km_s += other[5];
        me.srp.cr += other[6];
        me.drag.cd += other[7];
        me.fuel_mass_kg += other[8];

        me
    }
}

impl ConfigRepr for Spacecraft {}

impl Configurable for Spacecraft {
    type IntermediateRepr = Self;

    fn from_config(cfg: Self::IntermediateRepr, _cosm: Arc<Cosm>) -> Result<Self, ConfigError>
    where
        Self: Sized,
    {
        Ok(cfg)
    }

    fn to_config(&self) -> Result<Self::IntermediateRepr, ConfigError> {
        Ok(*self)
    }
}

#[test]
fn test_serde() {
    use super::Cosm;
    use serde_yaml;
    use std::str::FromStr;

    let cosm = Cosm::de438();
    let orbit = Orbit::cartesian(
        -9042.862234,
        18536.333069,
        6999.957069,
        -3.288789,
        -2.226285,
        1.646738,
        Epoch::from_str("2018-09-15T00:15:53.098 UTC").unwrap(),
        cosm.frame("EME2000"),
    );

    let sc = Spacecraft::new(orbit, 500.0, 159.0, 2.0, 2.0, 1.8, 2.2);

    let serialized_sc = serde_yaml::to_string(&sc).unwrap();
    println!("{}", serialized_sc);

    let deser_sc: Spacecraft = serde_yaml::from_str(&serialized_sc).unwrap();

    assert_eq!(sc, deser_sc);

    // Check that we can omit the thruster info entirely.
    let s = r#"
orbit:
    x_km: -9042.862234
    y_km: 18536.333069
    z_km: 6999.957069
    vx_km_s: -3.288789
    vy_km_s: -2.226285
    vz_km_s: 1.646738
    epoch: 2018-09-15T00:15:53.098000000 UTC
    frame: Earth J2000
dry_mass_kg: 500.0
fuel_mass_kg: 159.0
srp:
    area_m2: 2.0
    cr: 1.8
drag:
    area_m2: 2.0
    cd: 2.2
    "#;

    let deser_sc: Spacecraft = serde_yaml::from_str(s).unwrap();
    assert_eq!(sc, deser_sc);

    // Check that we can specify a thruster info entirely.
    let s = r#"
orbit:
    x_km: -9042.862234
    y_km: 18536.333069
    z_km: 6999.957069
    vx_km_s: -3.288789
    vy_km_s: -2.226285
    vz_km_s: 1.646738
    epoch: 2018-09-15T00:15:53.098000000 UTC
    frame: Earth J2000
dry_mass_kg: 500.0
fuel_mass_kg: 159.0
srp:
    area_m2: 2.0
    cr: 1.8
drag:
    area_m2: 2.0
    cd: 2.2
thruster:
    thrust_N: 1e-5
    isp_s: 300.5
    "#;

    let mut sc_thruster = sc;
    sc_thruster.thruster = Some(Thruster {
        isp_s: 300.5,
        thrust_N: 1e-5,
    });
    let deser_sc: Spacecraft = serde_yaml::from_str(s).unwrap();
    assert_eq!(sc_thruster, deser_sc);

    // Tests the minimum definition which will set all of the defaults too
    let s = r#"
orbit:
    x_km: -9042.862234
    y_km: 18536.333069
    z_km: 6999.957069
    vx_km_s: -3.288789
    vy_km_s: -2.226285
    vz_km_s: 1.646738
    epoch: 2018-09-15T00:15:53.098000000 UTC
    frame: Earth J2000
dry_mass_kg: 500.0
fuel_mass_kg: 159.0
"#;

    let deser_sc: Spacecraft = serde_yaml::from_str(s).unwrap();

    let sc = Spacecraft::new(orbit, 500.0, 159.0, 0.0, 0.0, 1.8, 2.2);
    assert_eq!(sc, deser_sc);
}
