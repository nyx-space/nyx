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

use super::formatter::OutputSerde;
use super::gravity::HarmonicsMem;
use super::quantity::*;
use super::rv::Distribution;
use super::serde_derive::Deserialize;
use super::ParsingError;
use crate::cosmic::{Frame, Orbit};
use crate::md::{Event, StateParameter};
use crate::time::{Duration, Epoch};
use crate::NyxError;
use std::collections::HashMap;
use std::str::FromStr;

#[derive(Deserialize)]
pub struct StateSerde {
    pub frame: Option<String>,
    pub epoch: String,
    pub x: Option<f64>,
    pub y: Option<f64>,
    pub z: Option<f64>,
    pub vx: Option<f64>,
    pub vy: Option<f64>,
    pub vz: Option<f64>,
    pub position: Option<Vec<String>>,
    pub velocity: Option<Vec<String>>,
    pub sma: Option<f64>,
    pub ecc: Option<f64>,
    pub inc: Option<f64>,
    pub raan: Option<f64>,
    pub aop: Option<f64>,
    pub ta: Option<f64>,
    pub unit_position: Option<String>,
    pub unit_velocity: Option<String>,
}

impl StateSerde {
    pub fn as_state(&self, frame: Frame) -> Result<Orbit, ParsingError> {
        let epoch = match Epoch::from_str(&self.epoch) {
            Ok(epoch) => epoch,
            Err(e) => return Err(ParsingError::Quantity(format!("{}", e))),
        };

        // Rebuild a valid state from the three different initializations
        if self.x.is_some()
            && self.y.is_some()
            && self.z.is_some()
            && self.vx.is_some()
            && self.vy.is_some()
            && self.vz.is_some()
        {
            let pos_mul = match &self.unit_position {
                Some(unit) => match unit.to_lowercase().as_str() {
                    "km" => 1.0,
                    "m" => 1e-3,
                    "cm" => 1e-5,
                    _ => {
                        return Err(ParsingError::Distance(format!(
                            "unknown unit `{}` when converting state",
                            unit
                        )))
                    }
                },
                None => 1.0,
            };

            let vel_mul = match &self.unit_velocity {
                Some(unit) => match unit.to_lowercase().as_str() {
                    "km/s" => 1.0,
                    "m/s" => 1e-3,
                    "cm/s" => 1e-5,
                    _ => {
                        return Err(ParsingError::Velocity(format!(
                            "unknown unit `{}` when converting state",
                            unit
                        )))
                    }
                },
                None => 1.0,
            };

            Ok(Orbit::cartesian(
                self.x.unwrap() * pos_mul,
                self.y.unwrap() * pos_mul,
                self.z.unwrap() * pos_mul,
                self.vx.unwrap() * vel_mul,
                self.vy.unwrap() * vel_mul,
                self.vz.unwrap() * vel_mul,
                epoch,
                frame,
            ))
        } else if self.position.is_some() && self.velocity.is_some() {
            let position = self.position.as_ref().unwrap();
            let velocity = self.velocity.as_ref().unwrap();
            if position.len() != 3 || velocity.len() != 3 {
                return Err(ParsingError::IllDefined(
                    "Orbit ill defined: position and velocity arrays must be of size 3 exactly"
                        .to_string(),
                ));
            }
            let x = parse_quantity(&position[0])?.v();
            let y = parse_quantity(&position[1])?.v();
            let z = parse_quantity(&position[2])?.v();
            let vx = parse_quantity(&velocity[0])?.v();
            let vy = parse_quantity(&velocity[1])?.v();
            let vz = parse_quantity(&velocity[2])?.v();

            Ok(Orbit::cartesian(x, y, z, vx, vy, vz, epoch, frame))
        } else if self.sma.is_some()
            && self.ecc.is_some()
            && self.inc.is_some()
            && self.raan.is_some()
            && self.aop.is_some()
            && self.ta.is_some()
        {
            Ok(Orbit::keplerian(
                self.sma.unwrap(),
                self.ecc.unwrap(),
                self.inc.unwrap(),
                self.raan.unwrap(),
                self.aop.unwrap(),
                self.ta.unwrap(),
                epoch,
                frame,
            ))
        } else {
            Err(ParsingError::IllDefined("Orbit ill defined: specify either {position, velocity}, or {x,y,z,vx,vy,vz}, or Keplerian elements".to_string()))
        }
    }
}

/// A state difference
#[derive(Deserialize)]
pub struct DeltaStateSerde {
    pub inherit: String,
    pub x: Option<f64>,
    pub y: Option<f64>,
    pub z: Option<f64>,
    pub vx: Option<f64>,
    pub vy: Option<f64>,
    pub vz: Option<f64>,
    pub position: Option<Vec<String>>,
    pub velocity: Option<Vec<String>>,
    pub sma: Option<f64>,
    pub ecc: Option<f64>,
    pub inc: Option<f64>,
    pub raan: Option<f64>,
    pub aop: Option<f64>,
    pub ta: Option<f64>,
    pub unit_position: Option<String>,
    pub unit_velocity: Option<String>,
}

impl DeltaStateSerde {
    pub fn as_state(&self, base: Orbit) -> Result<Orbit, ParsingError> {
        let frame = base.frame;
        let epoch = base.dt;
        // Rebuild a valid state from the three different initializations
        if self.x.is_some()
            || self.y.is_some()
            || self.z.is_some()
            || self.vx.is_some()
            || self.vy.is_some()
            || self.vz.is_some()
        {
            let pos_mul = match &self.unit_position {
                Some(unit) => match unit.to_lowercase().as_str() {
                    "km" => 1.0,
                    "m" => 1e-3,
                    "cm" => 1e-5,
                    _ => panic!("unknown unit `{}`", unit),
                },
                None => 1.0,
            };

            let vel_mul = match &self.unit_velocity {
                Some(unit) => match unit.to_lowercase().as_str() {
                    "km/s" => 1.0,
                    "m/s" => 1e-3,
                    "cm/s" => 1e-5,
                    _ => panic!("unknown unit `{}`", unit),
                },
                None => 1.0,
            };

            Ok(Orbit::cartesian(
                self.x.unwrap_or(0.0) * pos_mul,
                self.y.unwrap_or(0.0) * pos_mul,
                self.z.unwrap_or(0.0) * pos_mul,
                self.vx.unwrap_or(0.0) * vel_mul,
                self.vy.unwrap_or(0.0) * vel_mul,
                self.vz.unwrap_or(0.0) * vel_mul,
                epoch,
                frame,
            ) + base)
        } else if self.position.is_some() || self.velocity.is_some() {
            let (x, y, z) = match &self.position {
                Some(position) => {
                    if position.len() != 3 {
                        return Err(ParsingError::IllDefined(
                            "Orbit ill defined: position arrays must be of size 3 exactly"
                                .to_string(),
                        ));
                    }
                    let x = parse_quantity(&position[0])?.v();
                    let y = parse_quantity(&position[1])?.v();
                    let z = parse_quantity(&position[2])?.v();
                    (x, y, z)
                }
                None => (0.0, 0.0, 0.0),
            };

            let (vx, vy, vz) = match &self.velocity {
                Some(velocity) => {
                    if velocity.len() != 3 {
                        return Err(ParsingError::IllDefined(
                            "Orbit ill defined: velocity arrays must be of size 3 exactly"
                                .to_string(),
                        ));
                    }
                    let vx = parse_quantity(&velocity[0])?.v();
                    let vy = parse_quantity(&velocity[1])?.v();
                    let vz = parse_quantity(&velocity[2])?.v();
                    (vx, vy, vz)
                }
                None => (0.0, 0.0, 0.0),
            };

            Ok(Orbit::cartesian(x, y, z, vx, vy, vz, epoch, frame) + base)
        } else if self.sma.is_some()
            || self.ecc.is_some()
            || self.inc.is_some()
            || self.raan.is_some()
            || self.aop.is_some()
            || self.ta.is_some()
        {
            let sma = match self.sma {
                Some(sma) => sma + base.sma(),
                None => base.sma(),
            };
            let ecc = match self.ecc {
                Some(ecc) => ecc + base.ecc(),
                None => base.ecc(),
            };
            let inc = match self.inc {
                Some(inc) => inc + base.inc(),
                None => base.inc(),
            };
            let raan = match self.raan {
                Some(raan) => raan + base.raan(),
                None => base.raan(),
            };
            let aop = match self.aop {
                Some(aop) => aop + base.aop(),
                None => base.aop(),
            };
            let ta = match self.ta {
                Some(ta) => ta + base.ta(),
                None => base.ta(),
            };
            Ok(Orbit::keplerian(sma, ecc, inc, raan, aop, ta, epoch, frame) + base)
        } else {
            Err(ParsingError::IllDefined("Delta state ill defined: specify either {position, velocity}, or {x,y,z,vx,vy,vz}, or Keplerian elements".to_string()))
        }
    }
}

#[derive(Deserialize)]
pub struct OrbitalDynamicsSerde {
    pub initial_state: String,
    /// If unspecified, the frame of the state is assumed
    pub integration_frame: Option<String>,
    pub point_masses: Option<Vec<String>>,
    pub accel_models: Option<Vec<String>>,
}

#[derive(Deserialize)]
pub struct SpacecraftSerde {
    pub dry_mass: f64,
    pub fuel_mass: Option<f64>,
    pub orbital_dynamics: String,
    pub force_models: Option<Vec<String>>,
}

#[derive(Deserialize)]
pub struct SolarPressureSerde {
    pub sc_area: f64,
    #[serde(default = "default_cr")]
    pub cr: f64,
    #[serde(default = "default_phi")]
    pub phi: f64,
}

fn default_cr() -> f64 {
    1.8
}
fn default_phi() -> f64 {
    1367.0
}

#[derive(Deserialize)]
pub struct Harmonics {
    pub frame: String,
    pub degree: usize,
    pub order: Option<usize>,
    pub file: String,
}

impl Harmonics {
    pub fn load(&self) -> Result<HarmonicsMem, NyxError> {
        let gunzipped = self.file.contains("gz");
        let order = self.order.unwrap_or(0);
        if self.file.contains("cof") {
            HarmonicsMem::from_cof(self.file.as_str(), self.degree, order, gunzipped)
        } else if self.file.contains("sha") {
            HarmonicsMem::from_shadr(self.file.as_str(), self.degree, order, gunzipped)
        } else if self.file.contains("EGM") {
            HarmonicsMem::from_egm(self.file.as_str(), self.degree, order, gunzipped)
        } else {
            panic!("could not guess file format from name");
        }
    }
}

#[derive(Deserialize)]
pub struct AccelModel {
    pub harmonics: HashMap<String, Harmonics>,
}

#[derive(Deserialize)]
pub struct ForceModel {
    pub srp: HashMap<String, SolarPressureSerde>,
}

#[derive(Deserialize)]
#[serde(rename_all = "camelCase")]
pub enum PropagatorKind {
    Dormand45,
    Dormand78,
    Fehlberg45,
    CashKarp45,
    Rk89,
    Rk4,
    Verner56,
}

#[derive(Deserialize)]
pub struct PropagatorSerde {
    /// Name of the string being associated with this propagator
    pub dynamics: String,
    /// Name of the stoping condition used
    pub stop_cond: String,
    /// The tolerance of the propagator (default to 1e-12)
    pub tolerance: Option<f64>,
    pub output: Option<String>,
    /// If no kind is specified, an RK89 will be used
    pub kind: Option<PropagatorKind>,
}

#[derive(Deserialize)]
pub struct OdpSerde {
    /// Name of the navigation propagator
    pub navigation_prop: String,
    /// Name of the measurement generator/reader
    pub measurements: String,
    /// Name of the initial estimate name
    pub initial_estimate: String,
    /// Measurement noise
    pub msr_noise: Vec<f64>,
    /// Optional SNC
    pub snc: Option<Vec<f64>>,
    /// Optional SNC turn-off time
    pub snc_disable: Option<String>,
    /// Optional SNC exponential decay
    pub snc_decay: Option<Vec<String>>,
    /// Set the number of measurements to switch to an EKF
    pub ekf_msr_trigger: Option<usize>,
    /// Set the acceptable time between measurements
    pub ekf_disable_time: Option<Duration>,
    /// An optional output of a NavSolution
    pub output: Option<String>,
}

#[derive(Deserialize)]
pub struct EstimateSerde {
    /// The initial state estimate, must be a reference to a state
    pub state: String,
    /// The full covariance matrix, must be 6x6
    pub covar_mat: Option<Vec<Vec<f64>>>,
    /// The diagonal of the covariance matrix (a single vector of 6 values)
    pub covar_diag: Option<Vec<f64>>,
}

#[derive(Deserialize)]
pub struct MeasurementSerde {
    /// Name of the truth propagator, if unspecified then the output must be specified
    pub propagator: Option<String>,
    /// Names of the measurement devices to use
    pub msr_device: Vec<String>,
    /// Name of the output file to store the measurements
    pub output: Option<String>,
    /// Optionally specify whether to use the file if it exists
    pub use_file_if_available: Option<bool>,
}

#[derive(Deserialize)]
pub struct StationSerde {
    pub elevation: f64,
    pub range_noise: f64,
    pub range_rate_noise: f64,
    pub latitude: Option<f64>,
    pub longitude: Option<f64>,
    pub height: Option<f64>,
    pub inherit: Option<String>,
}

#[derive(Deserialize)]
pub struct ScenarioSerde {
    pub sequence: Vec<String>,
    pub propagator: HashMap<String, PropagatorSerde>,
    pub state: HashMap<String, StateSerde>,
    pub delta_state: Option<HashMap<String, DeltaStateSerde>>,
    pub orbital_dynamics: HashMap<String, OrbitalDynamicsSerde>,
    pub spacecraft: HashMap<String, SpacecraftSerde>,
    pub accel_models: Option<HashMap<String, AccelModel>>,
    pub force_models: Option<HashMap<String, ForceModel>>,
    pub output: HashMap<String, OutputSerde>,
    pub distr: Option<HashMap<String, Distribution>>,
    pub odp: Option<HashMap<String, OdpSerde>>,
    pub estimate: Option<HashMap<String, EstimateSerde>>,
    pub measurements: Option<HashMap<String, MeasurementSerde>>,
    pub stations: Option<HashMap<String, StationSerde>>,
    pub conditions: Option<HashMap<String, ConditionSerde>>,
}

#[derive(Clone, Deserialize)]
pub struct ConditionSerde {
    /// Pattern must be `event_name = event_value`, e.g. `TA = 159`
    pub event: String,
    pub search_until: String,
    pub hits: Option<usize>,
}

impl ConditionSerde {
    pub fn to_condition(&self) -> Event {
        let rplt = self.event.replace("=", "");
        let parts: Vec<&str> = rplt.split(' ').collect();
        let parameter = StateParameter::from_str(parts[0]).unwrap();
        let value = if parts.len() == 2 {
            match parts[1].trim().parse::<f64>() {
                Ok(val) => val,
                Err(e) => {
                    warn!(
                        "Could not understand value `{}` in parameter: {}",
                        parts[1], e
                    );
                    0.0
                }
            }
        } else {
            0.0
        };

        Event::new(parameter, value)
    }
}

#[test]
fn test_md_scenario() {
    extern crate toml;

    let scen: ScenarioSerde = toml::from_str(
        r#"
        sequence = ["prop_name"]

        [state.STATE_NAME]
        x = -2436.45
        y = -2436.45
        z = 6891.037
        vx = 5.088611
        vy = -5.088611
        vz = 0.0
        frame = "EME2000"
        epoch = "MJD 51544.5 TAI" # or "2018-09-15T00:15:53.098 UTC"
        unit_position = "km"  # Default value if unspecified
        unit_velocity = "km/s"  # Default value if unspecified

        [state.another_state]
        # Mix and match units
        position = ["-3436.45 km", "-3436.45 km", "6891.037e3 m"]
        velocity = ["5.088611 km/s", "-5.088611 km/s", "0.0 cm/s"]
        frame = "EME2000"
        epoch = "MJD 51540.5 TAI"

        [state.kepler_proud]
        sma = 7000
        inc = 27.4  # degrees
        ecc = 0.001
        raan = 45
        aop = 45
        ta = 180
        frame = "EME2000"
        epoch = "MJD 51540.5 TAI"

        [delta_state.my_delta]
        inherit = "another_state"
        x = 1.0  # Will add 1 km to the x component

        [delta_state.my_delta2]
        inherit = "another_state"
        velocity = ["1 m/s", "2 m/s", "-5.0 m/s"]

        [delta_state.my_delta3]
        inherit = "another_state"
        vx = 5e-3

        [orbital_dynamics.conf_name]
        integration_frame = "EME2000"
        initial_state = "state_name"
        point_masses = ["Sun", "Earth", "Jupiter", "Luna"]
        accel_models = ["jgm3_70x70"]

        [orbital_dynamics.two_body_dyn]
        initial_state = "state_name"

        [spacecraft.sc1]
        dry_mass = 100.0
        fuel_mass = 20.0
        orbital_dynamics = "conf_name"
        force_models = ["my_frc"]

        [propagator.prop_name]
        flavor = "rk89"  # If unspecified, the default propagator is used
        dynamics = "sc1"  # Use "sc1" spacecraft dynamics
        stop_cond = "MJD 51540.5 TAI"
        output = "my_csv"
        
        [output.my_csv]
        filename = "./data/scenario-run.csv"
        headers = ["epoch:GregorianUtc", "x", "y", "z", "vx", "vy", "vz", "rmag:Luna"]

        [propagator.simple]
        dynamics = "two_body"  # Use "sc1" spacecraft dynamics
        stop_cond = "1 * day"

        [accel_models.my_models.harmonics.jgm3_70x70]
        frame = "EME2000"
        degree = 70
        order = 70
        file = "data/JGM3.cof.gz"

        [force_models.my_frc.srp.my_srp]
        sc_area = 1.0 # in meters squared
        cr = 1.5 # Defaults to 1.8

        [condition.third_apo]
        kind = "apoapse"
        search_until = "MJD 51540.5 TAI"
        hits = 3  # Stopping condition triggered on third apoapse passage

        "#,
    )
    .unwrap();

    assert_eq!(scen.state.len(), 3);
    assert_eq!(scen.delta_state.unwrap().len(), 3);
    assert_eq!(scen.orbital_dynamics.len(), 2);
    assert_eq!(scen.propagator.len(), 2);
    assert_eq!(scen.accel_models.as_ref().unwrap().len(), 1);
    assert_eq!(
        scen.accel_models.as_ref().unwrap()["my_models"]
            .harmonics
            .len(),
        1
    );
    assert_eq!(scen.sequence.len(), 1);
}

#[test]
fn test_od_scenario() {
    extern crate toml;

    let scen: ScenarioSerde = toml::from_str(
        r#"
        sequence = ["my_flt"]

        [state.state_name]
        x = -2436.45
        y = -2436.45
        z = 6891.037
        vx = 5.088611
        vy = -5.088611
        vz = 0.0
        frame = "EME2000"
        epoch = "MJD 51544.5 TAI" # or "2018-09-15T00:15:53.098 UTC"
        unit_position = "km"  # Default value if unspecified
        unit_velocity = "km/s"  # Default value if unspecified

        [orbital_dynamics.conf_name]
        integration_frame = "EME2000"
        initial_state = "state_name"
        point_masses = ["Sun", "Earth", "Jupiter", "Luna"]
        accel_models = ["my_models"]

        [spacecraft.sc1]
        dry_mass = 100.0
        fuel_mass = 20.0
        orbital_dynamics = "conf_name"
        force_models = ["my_srp"]

        [propagator.nav_prop]
        dynamics = "sc1"
        stop_cond = "3.5 days"
        output = "my_csv"
        tolerance = 1e-9

        [propagator.truth_propagator]
        dynamics = "sc1"
        stop_cond = "3.5 days"
        output = "my_csv"

        [accel_models.my_models.harmonics.jgm3_70x70]
        frame = "EME2000"
        degree = 70
        order = 70
        file = "data/JGM3.cof.gz"

        [output.my_csv]
        filename = "./data/truth.csv"
        headers = ["epoch:GregorianUtc", "x", "y", "z", "vx", "vy", "vz", "rmag:Luna"]

        [odp.my_flt]
        navigation_prop = "nav_prop"
        initial_estimate = "my_estimate"
        msr_noise = [1e-6, 1e-3]
        snc = [1e-8, 1e-8, 1e-8]
        snc_disable = "120 * sec"
        snc_decay = ["20 * min", "20 min", "15 min"]
        measurements = "msr_sim"  # Or provide a file name
        ekf_msr_trigger = 30
        ekf_disable_time = "3600 s"  # If no measurements for an hour, disable the EKF
        output = "estimate_csv"

        [output.estimate_csv]
        filename = "./data/estimates.csv"
        headers = ["epoch:GregorianUtc", "delta_x", "delta_y", "delta_z", "delta_vx", "delta_vy", "delta_vz"]  # If unset, default will be used

        [measurements.msr_sim]
        propagator = "truth_propagator"
        msr_device = ["dss13", "dss65", "dss34"]
        output = "msr_sim.csv"
        use_file_if_available = true

        [stations.dss13]
        elevation = 10.0
        latitude = 40.427_222
        longitude = 4.250_556
        height = 0.834_939
        range_noise = 0.1
        range_rate_noise = 0.01

        [stations.dss65]
        inherit = "dss65" # Name of the station, built-in
        elevation = 10.0
        range_noise = 0.1
        range_rate_noise = 0.01

        [estimate.my_estimate]
        state = "state_name"
        #covar_mat = [[1e1, 1e1, 1e1, 1e1, 1e1, 1e1],[1e1, 1e1, 1e1, 1e1, 1e1, 1e1],[1e1, 1e1, 1e1, 1e1, 1e1, 1e1],[1e1, 1e1, 1e1, 1e1, 1e1, 1e1],[1e1, 1e1, 1e1, 1e1, 1e1, 1e1],[1e1, 1e1, 1e1, 1e1, 1e1, 1e1]]
        covar_diag = [1e1, 1e1, 1e1, 1e-2, 1e-2, 1e-2]
        "#,
    )
    .unwrap();

    assert_eq!(scen.state.len(), 1);
    assert_eq!(scen.orbital_dynamics.len(), 1);
    assert_eq!(scen.propagator.len(), 2);
    assert_eq!(scen.accel_models.as_ref().unwrap().len(), 1);
    assert_eq!(
        scen.accel_models.as_ref().unwrap()["my_models"]
            .harmonics
            .len(),
        1
    );
    assert_eq!(scen.sequence.len(), 1);
    assert_eq!(scen.odp.unwrap().len(), 1);
    assert_eq!(scen.measurements.unwrap().len(), 1);
    assert_eq!(scen.stations.unwrap().len(), 2);
    assert_eq!(scen.estimate.unwrap().len(), 1);
}
