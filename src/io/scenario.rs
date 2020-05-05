use super::gravity::HarmonicsMem;
use super::output::OutputSerde;
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
                "m" => 1_000.0,
                "cm" => 100_000.0,
                _ => panic!("unknown unit `{}`", unit),
            },
            None => 1.0,
        };

        let vel_mul = match &self.unit_velocity {
            Some(unit) => match unit.to_lowercase().as_str() {
                "km/s" => 1.0,
                "m/s" => 1_000.0,
                "cm/s" => 100_000.0,
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
    /// If unspecified, the STM will not be computed
    pub stm: Option<bool>,
}

impl OrbitalDynamicsSerde {
    pub fn with_stm(&self) -> bool {
        self.stm.is_some() && self.stm.unwrap()
    }
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
pub struct ScenarioSerde {
    pub sequence: Vec<String>,
    pub propagator: HashMap<String, PropagatorSerde>,
    pub state: HashMap<String, StateSerde>,
    pub orbital_dynamics: HashMap<String, OrbitalDynamicsSerde>,
    pub accel_models: HashMap<String, AccelModel>,
    pub output: HashMap<String, OutputSerde>,
    pub distr: Option<HashMap<String, Distribution>>,
}

#[test]
fn test_deser_scenario() {
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
        stm = false
        integration_frame = "EME2000"
        initial_state = "state_name"
        point_masses = ["Sun", "Earth", "Jupiter", "Luna"]
        accel_models = ["jgm3_70x70"]

        [orbital_dynamics.two_body]
        initial_state = "state_name"

        [propagator.prop_name]
        flavor = "rk89"  # If unspecified, the default propagator is used
        dynamics = "two_body"  # Use "sc1" spacecraft dynamics
        stop_cond = "MJD 51540.5 TAI"
        output = "my_output"

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
