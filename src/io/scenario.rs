use super::formatter::OutputSerde;
use super::gravity::HarmonicsMem;
use super::rv::Distribution;
use super::serde_derive::Deserialize;
use crate::celestia::{Frame, State};
use crate::time::Epoch;
use std::collections::HashMap;
use std::str::FromStr;

#[derive(Deserialize)]
pub struct StateSerde {
    pub x: f64,
    pub y: f64,
    pub z: f64,
    pub vx: f64,
    pub vy: f64,
    pub vz: f64,
    pub frame: String,
    pub epoch: String,
    pub unit_position: Option<String>,
    pub unit_velocity: Option<String>,
}

impl StateSerde {
    pub fn as_state(&self, frame: Frame) -> State {
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

        State::cartesian(
            self.x * pos_mul,
            self.y * pos_mul,
            self.z * pos_mul,
            self.vx * vel_mul,
            self.vy * vel_mul,
            self.vz * vel_mul,
            Epoch::from_str(&self.epoch).unwrap(),
            frame,
        )
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
    pub fn load(&self) -> HarmonicsMem {
        let gunzipped = self.file.contains("gz");
        let order = match self.order {
            Some(order) => order,
            None => 0,
        };
        if self.file.contains("cof") {
            HarmonicsMem::from_cof(self.file.as_str(), self.degree, order, gunzipped).unwrap()
        } else if self.file.contains("sha") {
            HarmonicsMem::from_shadr(self.file.as_str(), self.degree, order, gunzipped).unwrap()
        } else if self.file.contains("EGM") {
            HarmonicsMem::from_egm(self.file.as_str(), self.degree, order, gunzipped).unwrap()
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
    pub output: Option<String>,
    /// If no kind is specified, an RK89 will be used
    pub kind: Option<PropagatorKind>,
}

#[derive(Deserialize)]
pub struct OdpSerde {
    /// The kind of filter to use
    pub kind: String,
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
    /// Set the number of measurements to switch to an EKF
    pub ekf_msr_trigger: Option<usize>,
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
    pub output: String,
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
    pub from: Option<String>,
}

#[derive(Deserialize)]
pub struct ScenarioSerde {
    pub sequence: Vec<String>,
    pub propagator: HashMap<String, PropagatorSerde>,
    pub state: HashMap<String, StateSerde>,
    pub orbital_dynamics: HashMap<String, OrbitalDynamicsSerde>,
    pub spacecraft: HashMap<String, SpacecraftSerde>,
    pub accel_models: HashMap<String, AccelModel>,
    pub output: HashMap<String, OutputSerde>,
    pub distr: Option<HashMap<String, Distribution>>,
    pub odp: Option<HashMap<String, OdpSerde>>,
    pub estimate: Option<HashMap<String, EstimateSerde>>,
    pub measurements: Option<HashMap<String, MeasurementSerde>>,
    pub stations: Option<HashMap<String, StationSerde>>,
}

#[test]
fn test_md_scenario() {
    extern crate toml;

    let scen: ScenarioSerde = toml::from_str(
        r#"
        sequence = ["prop_name"]

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

        [state.another_state]
        x = -3436.45
        y = -3436.45
        z = 6891.037
        vx = 5.088611
        vy = -5.088611
        vz = 0.0
        frame = "EME2000"
        epoch = "MJD 51540.5 TAI"

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
        "#,
    )
    .unwrap();

    assert_eq!(scen.state.len(), 2);
    assert_eq!(scen.orbital_dynamics.len(), 2);
    assert_eq!(scen.propagator.len(), 2);
    assert_eq!(scen.accel_models.len(), 1);
    assert_eq!(scen.accel_models["my_models"].harmonics.len(), 1);
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

        [propagator.nav_prop]
        dynamics = "sc1"
        stop_cond = "3.5 days"
        output = "my_csv"

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
        kind = "cfk"  # Could be ekf or srif
        navigation_prop = "nav_prop"
        initial_estimate = "my_estimate"
        msr_noise = [1e-6, 1e-3]
        snc = [1e-8, 1e-8, 1e-8]
        snc_disable = "120 * sec"
        measurements = "msr_sim"  # Or provide a file name
        ekf_msr_trigger = 30
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
        from = "dss65" # Name of the station, built-in
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
    assert_eq!(scen.accel_models.len(), 1);
    assert_eq!(scen.accel_models["my_models"].harmonics.len(), 1);
    assert_eq!(scen.sequence.len(), 1);
    assert_eq!(scen.odp.unwrap().len(), 1);
    assert_eq!(scen.measurements.unwrap().len(), 1);
    assert_eq!(scen.stations.unwrap().len(), 2);
    assert_eq!(scen.estimate.unwrap().len(), 1);
}
