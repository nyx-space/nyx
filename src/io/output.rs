use super::serde_derive::Deserialize;
use super::EpochFormat;
use crate::celestia::{Cosm, Frame, State};
use std::collections::HashMap;
use std::str::FromStr;

#[derive(Deserialize)]
pub struct OutputSerde {
    pub filename: String,
    /// If not specified, the standard
    pub headers: Option<Vec<String>>,
}

impl OutputSerde {
    pub fn to_state_formatter<'a>(&self, cosm: &'a Cosm) -> StateFormatter<'a> {
        StateFormatter::new(
            self.filename.clone(),
            match &self.headers {
                Some(hdr) => hdr.clone(),
                None => vec![
                    "epoch".to_string(),
                    "x".to_string(),
                    "y".to_string(),
                    "z".to_string(),
                    "vx".to_string(),
                    "vy".to_string(),
                    "vz".to_string(),
                ],
            },
            cosm,
        )
    }
}

pub struct StateFormatter<'a> {
    pub filename: String,
    pub headers: Vec<String>,
    frames: HashMap<String, Frame>,
    cosm: &'a Cosm,
}

impl<'a> StateFormatter<'a> {
    pub fn new(filename: String, headers: Vec<String>, cosm: &'a Cosm) -> Self {
        let mut frames = HashMap::new();
        for hdr in &headers {
            let splt: Vec<&str> = hdr.split(':').collect();
            if splt[0].to_lowercase() != "epoch" && splt.len() == 2 {
                // Get the frame
                match cosm.try_frame(splt[1]) {
                    Ok(frame) => frames.insert(splt[1].to_string(), frame),
                    Err(e) => panic!("unknown frame `{}` in header ({})", splt[1], e),
                };
            }
        }

        Self {
            filename,
            headers: headers.iter().map(|h| h.to_lowercase()).collect(),
            frames,
            cosm,
        }
    }

    pub fn format(&self, state: &State) -> Vec<String> {
        // Start by computing the state in all of the frames needed
        let mut mapped = HashMap::new();
        for (name, frame) in &self.frames {
            mapped.insert(name.to_lowercase(), self.cosm.frame_chg(state, *frame));
        }
        let mut formatted = Vec::new();

        for hdr in &self.headers {
            let splt: Vec<&str> = hdr.split(':').collect();
            match splt[0] {
                "epoch" => {
                    let efmt = if splt.len() == 2 {
                        EpochFormat::from_str(splt[1]).unwrap()
                    } else {
                        EpochFormat::GregorianUtc
                    };
                    formatted.push(efmt.format(state.dt));
                }
                _ => {
                    // Grab the state in the other frame if needed
                    let out_state = if splt.len() == 2 {
                        &mapped[&splt[1].to_string()]
                    } else {
                        state
                    };
                    formatted.push(match splt[0] {
                        "aol" => format!("{:.9}", out_state.aol()),
                        "aop" => format!("{:.9}", out_state.aop()),
                        "apoapsis" => format!("{:.9}", out_state.apoapsis()),
                        "ea" => format!("{:.9}", out_state.ea()),
                        "ecc" => format!("{:.9}", out_state.ecc()),
                        "energy" => format!("{:.9}", out_state.energy()),
                        "evec" => format!("{:.9}", out_state.evec()),
                        "geodetic_height" => format!("{:.9}", out_state.geodetic_height()),
                        "geodetic_latitude" => format!("{:.9}", out_state.geodetic_latitude()),
                        "geodetic_longitude" => format!("{:.9}", out_state.geodetic_longitude()),
                        "hmag" => format!("{:.9}", out_state.hmag()),
                        "hvec" => format!("{:.9}", out_state.hvec()),
                        "hx" => format!("{:.9}", out_state.hx()),
                        "hy" => format!("{:.9}", out_state.hy()),
                        "hz" => format!("{:.9}", out_state.hz()),
                        "inc" => format!("{:.9}", out_state.inc()),
                        "ma" => format!("{:.9}", out_state.ma()),
                        "periapsis" => format!("{:.9}", out_state.periapsis()),
                        "period" => format!("{:.9}", out_state.period()),
                        "r_hat" => format!("{:.9}", out_state.r_hat()),
                        "raan" => format!("{:.9}", out_state.raan()),
                        "radius" => format!("{:.9}", out_state.radius()),
                        "rmag" => format!("{:.9}", out_state.rmag()),
                        "semi_parameter" => format!("{:.9}", out_state.semi_parameter()),
                        "sma" => format!("{:.9}", out_state.sma()),
                        "ta" => format!("{:.9}", out_state.ta()),
                        "tlong" => format!("{:.9}", out_state.tlong()),
                        "v_hat" => format!("{:.9}", out_state.v_hat()),
                        "velocity" => format!("{:.9}", out_state.velocity()),
                        "vmag" => format!("{:.9}", out_state.vmag()),

                        "x" => format!("{:.9}", out_state.x),
                        "y" => format!("{:.9}", out_state.y),
                        "z" => format!("{:.9}", out_state.z),
                        "vx" => format!("{:.9}", out_state.vx),
                        "vy" => format!("{:.9}", out_state.vy),
                        "vz" => format!("{:.9}", out_state.vz),
                        _ => panic!("unknown header `{}`", splt[0]),
                    });
                }
            }
        }

        formatted
    }
}
