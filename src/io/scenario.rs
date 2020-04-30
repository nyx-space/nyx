use super::rv::Distribution;
use super::serde_derive::Deserialize;
use std::collections::HashMap;

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
pub struct Scenario {
    sequence: Vec<String>,
    propagator: HashMap<String, PropagatorSerde>,
    state: HashMap<String, StateSerde>,
    orbital_dynamics: HashMap<String, OrbitalDynamicsSerde>,
    distr: Option<HashMap<String, Distribution>>,
}

impl Scenario {
    pub fn validate(&self) {
        for seq_name in &self.sequence {
            match self.propagator.get(seq_name) {
                None => panic!("sequence refers to undefined propagator `{}` ", seq_name),
                Some(prop) => match self.orbital_dynamics.get(&prop.dynamics) {
                    None => panic!(
                        "propagator `{}` refers to undefined dynamics `{}`",
                        seq_name, prop.dynamics
                    ),
                    Some(dynamics) => {
                        if self.state.get(&dynamics.initial_state).is_none() {
                            panic!(
                                "dynamics `{}` refers to unknown state `{}`",
                                prop.dynamics, dynamics.initial_state
                            );
                        }
                    }
                },
            }
        }
    }
}

#[test]
fn test_deser_scenario() {
    extern crate toml;

    let scen: Scenario = toml::from_str(
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
        stop_cond = "one_day_prop"
        output = "my_output"

        [propagator.simple]
        dynamics = "two_body"  # Use "sc1" spacecraft dynamics
        stop_cond = "one_day_prop"
        "#,
    )
    .unwrap();

    assert_eq!(scen.state.len(), 2);
    assert_eq!(scen.orbital_dynamics.len(), 2);
    assert_eq!(scen.propagator.len(), 2);
    assert_eq!(scen.sequence.len(), 1);
}
