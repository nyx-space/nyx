use super::serde::ser::SerializeSeq;
use super::serde::{Serialize, Serializer};
use super::serde_derive::Deserialize;
use super::EpochFormat;
use crate::celestia::{Cosm, Frame, State};
use crate::dimensions::U6;
use crate::od::estimate::NavSolution;
use crate::od::EstimableState;
use std::cmp::PartialEq;
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
        match &self.headers {
            Some(hdr) => StateFormatter::from_headers(hdr.to_vec(), self.filename.clone(), cosm),
            None => StateFormatter::default(self.filename.clone(), cosm),
        }
    }
}

/// Allowed headers, with an optional frame.
/// TODO: Support units
#[allow(non_camel_case_types)]
#[derive(Clone, Debug, PartialEq)]
pub enum StateHeader {
    /// Argument of Periapse (deg)
    AoL { frame: Option<String> },
    /// Argument of Latitude (deg)
    AoP { frame: Option<String> },
    /// Radius of apoapsis (km)
    apoapsis { frame: Option<String> },
    /// Eccentric anomaly (deg)
    EA { frame: Option<String> },
    /// Eccentricity (no unit)
    ECC { frame: Option<String> },
    /// The epoch in the specified format
    Epoch(EpochFormat),
    /// Specific energy
    energy { frame: Option<String> },
    /// Eccentricity vector (no unit), as [e_x,e_y,e_z]
    evec { frame: Option<String> },
    /// Geodetic height (km)
    geodetic_height { frame: Option<String> },
    /// Geodetic latitude (deg)
    geodetic_latitude { frame: Option<String> },
    /// Geodetic longitude (deg)
    geodetic_longitude { frame: Option<String> },
    /// Orbital momentum
    hmag { frame: Option<String> },
    /// Orbital momentum vector, as [e_x,e_y,e_z]
    hvec { frame: Option<String> },
    /// X component of the orbital momentum vector
    HX { frame: Option<String> },
    /// Y component of the orbital momentum vector
    HY { frame: Option<String> },
    /// Z component of the orbital momentum vector
    HZ { frame: Option<String> },
    /// Inclination (deg)
    INC { frame: Option<String> },
    /// Mean anomaly (deg)
    MA { frame: Option<String> },
    /// Radius of periapse (km)
    periapsis { frame: Option<String> },
    /// Orbital period (s)
    period { frame: Option<String> },
    /// Right ascension of the ascending node (deg)
    RAAN { frame: Option<String> },
    /// Radius vector (km), as [r_x,r_y,r_z]
    radius { frame: Option<String> },
    /// Norm of the radius vector
    rmag { frame: Option<String> },
    /// Semi parameter (km)
    semi_parameter { frame: Option<String> },
    /// Semi major axis (km)
    SMA { frame: Option<String> },
    /// True anomaly
    TA { frame: Option<String> },
    /// True longitude
    TLong { frame: Option<String> },
    /// Velocity vector (km/s), as [v_x,v_y,v_z]
    velocity { frame: Option<String> },
    /// Norm of the velocity vector (km/s)
    vmag { frame: Option<String> },
    /// X component of the radius (km)
    X { frame: Option<String> },
    /// Y component of the radius (km)
    Y { frame: Option<String> },
    /// Z component of the radius (km)
    Z { frame: Option<String> },
    /// X component of the velocity (km/s)
    VX { frame: Option<String> },
    /// Y component of the velocity (km/s)
    VY { frame: Option<String> },
    /// Z component of the velocity (km/s)
    VZ { frame: Option<String> },
}

impl Serialize for StateHeader {
    /// NOTE: This is not part of unit testing because there is no deseralization of State (yet)
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self {
            StateHeader::AoL { .. } => serializer.serialize_str("AoL"),
            StateHeader::AoP { .. } => serializer.serialize_str("AoP"),
            StateHeader::apoapsis { .. } => serializer.serialize_str("apoapsis"),
            StateHeader::EA { .. } => serializer.serialize_str("EA"),
            StateHeader::ECC { .. } => serializer.serialize_str("ECC"),
            StateHeader::energy { .. } => serializer.serialize_str("energy"),
            StateHeader::evec { .. } => serializer.serialize_str("evec"),
            StateHeader::geodetic_height { .. } => serializer.serialize_str("geodetic_height"),
            StateHeader::geodetic_latitude { .. } => serializer.serialize_str("geodetic_latitude"),
            StateHeader::geodetic_longitude { .. } => {
                serializer.serialize_str("geodetic_longitude")
            }
            StateHeader::hmag { .. } => serializer.serialize_str("hmag"),
            StateHeader::hvec { .. } => serializer.serialize_str("hvec"),
            StateHeader::HX { .. } => serializer.serialize_str("HX"),
            StateHeader::HY { .. } => serializer.serialize_str("HY"),
            StateHeader::HZ { .. } => serializer.serialize_str("HZ"),
            StateHeader::INC { .. } => serializer.serialize_str("INC"),
            StateHeader::MA { .. } => serializer.serialize_str("MA"),
            StateHeader::periapsis { .. } => serializer.serialize_str("periapsis"),
            StateHeader::period { .. } => serializer.serialize_str("period"),
            StateHeader::RAAN { .. } => serializer.serialize_str("RAAN"),
            StateHeader::radius { .. } => serializer.serialize_str("radius"),
            StateHeader::rmag { .. } => serializer.serialize_str("rmag"),
            StateHeader::semi_parameter { .. } => serializer.serialize_str("semi_parameter"),
            StateHeader::SMA { .. } => serializer.serialize_str("SMA"),
            StateHeader::TA { .. } => serializer.serialize_str("TA"),
            StateHeader::TLong { .. } => serializer.serialize_str("TLong"),
            StateHeader::velocity { .. } => serializer.serialize_str("velocity"),
            StateHeader::vmag { .. } => serializer.serialize_str("vmag"),
            StateHeader::X { .. } => serializer.serialize_str("X"),
            StateHeader::Y { .. } => serializer.serialize_str("Y"),
            StateHeader::Z { .. } => serializer.serialize_str("Z"),
            StateHeader::VX { .. } => serializer.serialize_str("VX"),
            StateHeader::VY { .. } => serializer.serialize_str("VY"),
            StateHeader::VZ { .. } => serializer.serialize_str("VZ"),
            StateHeader::Epoch(efmt) => {
                serializer.serialize_str(format!("Epoch:{:?}", efmt).as_str())
            }
        }
    }
}

/// Allowed headers, with an optional frame.
/// TODO: Support units
#[allow(non_camel_case_types)]
#[derive(Clone, Debug, PartialEq)]
pub enum NavSolutionHeader {
    /// The epoch in the specified format
    Epoch(EpochFormat),
    /// Headers for the estimated state
    EstimatedState(Vec<StateHeader>),
    /// Headers for the nominal state
    NominalState(Vec<StateHeader>),
    /// State deviation X (km)
    Delta_x,
    /// State deviation Y (km)
    Delta_y,
    /// State deviation Z (km)
    Delta_z,
    /// State deviation VX (km/s)
    Delta_vx,
    /// State deviation VY (km/s)
    Delta_vy,
    /// State deviation VZ (km/s)
    Delta_vz,
    /// Covariance matrix [1,1]
    Cx_x,
    /// Covariance matrix [2,1]
    Cy_x,
    /// Covariance matrix [2,2]
    Cy_y,
    /// Covariance matrix [3,1]
    Cz_x,
    /// Covariance matrix [3,2]
    Cz_y,
    /// Covariance matrix [3,3]
    Cz_z,
    /// Covariance matrix [4,1]
    Cx_dot_x,
    /// Covariance matrix [4,2]
    Cx_dot_y,
    /// Covariance matrix [4,3]
    Cx_dot_z,
    /// Covariance matrix [4,4]
    Cx_dot_x_dot,
    /// Covariance matrix [5,1]
    Cy_dot_x,
    /// Covariance matrix [5,2]
    Cy_dot_y,
    /// Covariance matrix [5,3]
    Cy_dot_z,
    /// Covariance matrix [5,4]
    Cy_dot_x_dot,
    /// Covariance matrix [5,5]
    Cy_dot_y_dot,
    /// Covariance matrix [6,1]
    Cz_dot_x,
    /// Covariance matrix [6,2]
    Cz_dot_y,
    /// Covariance matrix [6,3]
    Cz_dot_z,
    /// Covariance matrix [6,4]
    Cz_dot_x_dot,
    /// Covariance matrix [6,5]
    Cz_dot_y_dot,
    /// Covariance matrix [6,6]
    Cz_dot_z_dot,
}

impl Serialize for NavSolutionHeader {
    /// NOTE: This is not part of unit testing because there is no deseralization of State (yet)
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self {
            NavSolutionHeader::Epoch(efmt) => {
                serializer.serialize_str(format!("Epoch:{:?}", efmt).as_str())
            }
            NavSolutionHeader::EstimatedState(hdr) => {
                let mut seq = serializer.serialize_seq(Some(hdr.len()))?;
                for element in hdr {
                    seq.serialize_element(element)?;
                }
                seq.end()
            }
            NavSolutionHeader::NominalState(hdr) => {
                let mut seq = serializer.serialize_seq(Some(hdr.len()))?;
                for element in hdr {
                    seq.serialize_element(element)?;
                }
                seq.end()
            }
            NavSolutionHeader::Delta_x => serializer.serialize_str("delta_x"),
            NavSolutionHeader::Delta_y => serializer.serialize_str("delta_y"),
            NavSolutionHeader::Delta_z => serializer.serialize_str("delta_z"),
            NavSolutionHeader::Delta_vx => serializer.serialize_str("delta_vx"),
            NavSolutionHeader::Delta_vy => serializer.serialize_str("delta_vy"),
            NavSolutionHeader::Delta_vz => serializer.serialize_str("delta_vz"),
            NavSolutionHeader::Cx_x => serializer.serialize_str("cx_x"),
            NavSolutionHeader::Cy_x => serializer.serialize_str("cy_x"),
            NavSolutionHeader::Cy_y => serializer.serialize_str("cy_y"),
            NavSolutionHeader::Cz_x => serializer.serialize_str("cz_x"),
            NavSolutionHeader::Cz_y => serializer.serialize_str("cz_y"),
            NavSolutionHeader::Cz_z => serializer.serialize_str("cz_z"),
            NavSolutionHeader::Cx_dot_x => serializer.serialize_str("cx_dot_x"),
            NavSolutionHeader::Cx_dot_y => serializer.serialize_str("cx_dot_y"),
            NavSolutionHeader::Cx_dot_z => serializer.serialize_str("cx_dot_z"),
            NavSolutionHeader::Cx_dot_x_dot => serializer.serialize_str("cx_dot_x_dot"),
            NavSolutionHeader::Cy_dot_x => serializer.serialize_str("cy_dot_x"),
            NavSolutionHeader::Cy_dot_y => serializer.serialize_str("cy_dot_y"),
            NavSolutionHeader::Cy_dot_z => serializer.serialize_str("cy_dot_z"),
            NavSolutionHeader::Cy_dot_x_dot => serializer.serialize_str("cy_dot_x_dot"),
            NavSolutionHeader::Cy_dot_y_dot => serializer.serialize_str("cy_dot_y_dot"),
            NavSolutionHeader::Cz_dot_x => serializer.serialize_str("cz_dot_x"),
            NavSolutionHeader::Cz_dot_y => serializer.serialize_str("cz_dot_y"),
            NavSolutionHeader::Cz_dot_z => serializer.serialize_str("cz_dot_z"),
            NavSolutionHeader::Cz_dot_x_dot => serializer.serialize_str("cz_dot_x_dot"),
            NavSolutionHeader::Cz_dot_y_dot => serializer.serialize_str("cz_dot_y_dot"),
            NavSolutionHeader::Cz_dot_z_dot => serializer.serialize_str("cz_dot_z_dot"),
        }
    }
}

/// A formatter for states
pub struct StateFormatter<'a> {
    pub filename: String,
    pub headers: Vec<StateHeader>,
    frames: HashMap<String, Frame>,
    cosm: &'a Cosm,
}

impl<'a> StateFormatter<'a> {
    pub fn from_headers(headers: Vec<String>, filename: String, cosm: &'a Cosm) -> Self {
        let mut frames = HashMap::new();
        let mut hdrs = Vec::with_capacity(20);
        // Rebuild the header tokens
        for hdr in &headers {
            let splt: Vec<&str> = hdr.split(':').collect();

            match splt[0] {
                "epoch" => {
                    hdrs.push(StateHeader::Epoch(if splt.len() == 2 {
                        EpochFormat::from_str(splt[1]).unwrap()
                    } else {
                        EpochFormat::GregorianUtc
                    }));
                }
                _ => {
                    let frame_name = if splt.len() == 2 {
                        Some(splt[1].to_owned())
                    } else {
                        None
                    };
                    hdrs.push(match splt[0].to_lowercase().as_str() {
                        "aol" => StateHeader::AoL { frame: frame_name },
                        "aop" => StateHeader::AoP { frame: frame_name },
                        "apoapsis" => StateHeader::apoapsis { frame: frame_name },
                        "ea" => StateHeader::EA { frame: frame_name },
                        "ecc" => StateHeader::ECC { frame: frame_name },
                        "energy" => StateHeader::energy { frame: frame_name },
                        "evec" => StateHeader::evec { frame: frame_name },
                        "geodetic_height" => StateHeader::geodetic_height { frame: frame_name },
                        "geodetic_latitude" => StateHeader::geodetic_latitude { frame: frame_name },
                        "geodetic_longitude" => {
                            StateHeader::geodetic_longitude { frame: frame_name }
                        }
                        "hmag" => StateHeader::hmag { frame: frame_name },
                        "hvec" => StateHeader::hvec { frame: frame_name },
                        "hx" => StateHeader::HX { frame: frame_name },
                        "hy" => StateHeader::HY { frame: frame_name },
                        "hz" => StateHeader::HZ { frame: frame_name },
                        "inc" => StateHeader::INC { frame: frame_name },
                        "ma" => StateHeader::MA { frame: frame_name },
                        "periapsis" => StateHeader::periapsis { frame: frame_name },
                        "period" => StateHeader::period { frame: frame_name },
                        "raan" => StateHeader::RAAN { frame: frame_name },
                        "radius" => StateHeader::radius { frame: frame_name },
                        "rmag" => StateHeader::rmag { frame: frame_name },
                        "semi_parameter" => StateHeader::semi_parameter { frame: frame_name },
                        "sma" => StateHeader::SMA { frame: frame_name },
                        "ta" => StateHeader::TA { frame: frame_name },
                        "tlong" => StateHeader::TLong { frame: frame_name },
                        "velocity" => StateHeader::velocity { frame: frame_name },
                        "vmag" => StateHeader::vmag { frame: frame_name },
                        "x" => StateHeader::X { frame: frame_name },
                        "y" => StateHeader::Y { frame: frame_name },
                        "z" => StateHeader::Z { frame: frame_name },
                        "vx" => StateHeader::VX { frame: frame_name },
                        "vy" => StateHeader::VY { frame: frame_name },
                        "vz" => StateHeader::VZ { frame: frame_name },
                        _ => panic!("unknown header `{}`", splt[0]),
                    });
                }
            }

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
            headers: hdrs,
            frames,
            cosm,
        }
    }

    /// Default headers are [Epoch (GregorianTai), X, Y, Z, VX, VY, VZ], where position is in km and velocity in km/s.
    pub fn default(filename: String, cosm: &'a Cosm) -> Self {
        Self {
            filename,
            headers: vec![
                StateHeader::Epoch(EpochFormat::GregorianTai),
                StateHeader::X { frame: None },
                StateHeader::Y { frame: None },
                StateHeader::Z { frame: None },
                StateHeader::VX { frame: None },
                StateHeader::VY { frame: None },
                StateHeader::VZ { frame: None },
            ],
            frames: HashMap::new(),
            cosm,
        }
    }

    pub fn fmt(&self, state: &State) -> Vec<String> {
        // Start by computing the state in all of the frames needed
        let mut mapped = HashMap::new();
        for (name, frame) in &self.frames {
            mapped.insert(name.to_lowercase(), self.cosm.frame_chg(state, *frame));
        }
        let mut formatted = Vec::new();

        for hdr in &self.headers {
            match hdr {
                StateHeader::Epoch(efmt) => formatted.push(efmt.format(state.dt)),
                StateHeader::AoL { frame }
                | StateHeader::AoP { frame }
                | StateHeader::apoapsis { frame }
                | StateHeader::EA { frame }
                | StateHeader::ECC { frame }
                | StateHeader::energy { frame }
                | StateHeader::evec { frame }
                | StateHeader::geodetic_height { frame }
                | StateHeader::geodetic_latitude { frame }
                | StateHeader::geodetic_longitude { frame }
                | StateHeader::hmag { frame }
                | StateHeader::hvec { frame }
                | StateHeader::HX { frame }
                | StateHeader::HY { frame }
                | StateHeader::HZ { frame }
                | StateHeader::INC { frame }
                | StateHeader::MA { frame }
                | StateHeader::periapsis { frame }
                | StateHeader::period { frame }
                | StateHeader::RAAN { frame }
                | StateHeader::radius { frame }
                | StateHeader::rmag { frame }
                | StateHeader::semi_parameter { frame }
                | StateHeader::SMA { frame }
                | StateHeader::TA { frame }
                | StateHeader::TLong { frame }
                | StateHeader::velocity { frame }
                | StateHeader::vmag { frame }
                | StateHeader::X { frame }
                | StateHeader::Y { frame }
                | StateHeader::Z { frame }
                | StateHeader::VX { frame }
                | StateHeader::VY { frame }
                | StateHeader::VZ { frame } => {
                    // Grab the state in the other frame if needed
                    let out_state = if frame.is_some() {
                        &mapped[&frame.as_ref().unwrap().to_lowercase()]
                    } else {
                        state
                    };

                    formatted.push(match hdr {
                        StateHeader::AoL { .. } => format!("{:.16e}", out_state.aol()),
                        StateHeader::AoP { .. } => format!("{:.16e}", out_state.aop()),
                        StateHeader::apoapsis { .. } => format!("{:.16e}", out_state.apoapsis()),
                        StateHeader::EA { .. } => format!("{:.16e}", out_state.ea()),
                        StateHeader::ECC { .. } => format!("{:.16e}", out_state.ecc()),
                        StateHeader::energy { .. } => format!("{:.16e}", out_state.energy()),
                        StateHeader::evec { .. } => format!(
                            "[{:.16e},{:.16e},{:.16e}]",
                            out_state.evec()[0],
                            out_state.evec()[1],
                            out_state.evec()[2]
                        ),
                        StateHeader::geodetic_height { .. } => {
                            format!("{:.16e}", out_state.geodetic_height())
                        }
                        StateHeader::geodetic_latitude { .. } => {
                            format!("{:.16e}", out_state.geodetic_latitude())
                        }
                        StateHeader::geodetic_longitude { .. } => {
                            format!("{:.16e}", out_state.geodetic_longitude())
                        }
                        StateHeader::hmag { .. } => format!("{:.16e}", out_state.hmag()),
                        StateHeader::hvec { .. } => format!(
                            "[{:.16e},{:.16e},{:.16e}]",
                            out_state.hvec()[0],
                            out_state.hvec()[1],
                            out_state.hvec()[2]
                        ),
                        StateHeader::HX { .. } => format!("{:.16e}", out_state.hx()),
                        StateHeader::HY { .. } => format!("{:.16e}", out_state.hy()),
                        StateHeader::HZ { .. } => format!("{:.16e}", out_state.hz()),
                        StateHeader::INC { .. } => format!("{:.16e}", out_state.inc()),
                        StateHeader::MA { .. } => format!("{:.16e}", out_state.ma()),
                        StateHeader::periapsis { .. } => format!("{:.16e}", out_state.periapsis()),
                        StateHeader::period { .. } => format!("{:.16e}", out_state.period()),
                        StateHeader::RAAN { .. } => format!("{:.16e}", out_state.raan()),
                        StateHeader::radius { .. } => format!(
                            "[{:.16e},{:.16e},{:.16e}]",
                            out_state.radius()[0],
                            out_state.radius()[1],
                            out_state.radius()[2]
                        ),
                        StateHeader::rmag { .. } => format!("{:.16e}", out_state.rmag()),
                        StateHeader::semi_parameter { .. } => {
                            format!("{:.16e}", out_state.semi_parameter())
                        }
                        StateHeader::SMA { .. } => format!("{:.16e}", out_state.sma()),
                        StateHeader::TA { .. } => format!("{:.16e}", out_state.ta()),
                        StateHeader::TLong { .. } => format!("{:.16e}", out_state.tlong()),
                        StateHeader::velocity { .. } => format!(
                            "[{:.16e},{:.16e},{:.16e}]",
                            out_state.velocity()[0],
                            out_state.velocity()[1],
                            out_state.velocity()[2]
                        ),
                        StateHeader::vmag { .. } => format!("{:.16e}", out_state.vmag()),
                        StateHeader::X { .. } => format!("{:.16e}", out_state.x),
                        StateHeader::Y { .. } => format!("{:.16e}", out_state.y),
                        StateHeader::Z { .. } => format!("{:.16e}", out_state.z),
                        StateHeader::VX { .. } => format!("{:.16e}", out_state.vx),
                        StateHeader::VY { .. } => format!("{:.16e}", out_state.vy),
                        StateHeader::VZ { .. } => format!("{:.16e}", out_state.vz),
                        _ => panic!("unsupported header `{:?}`", hdr),
                    });
                }
            };
        }

        formatted
    }
}

/// A formatter for navigation solution
pub struct NavSolutionFormatter<'a> {
    pub filename: String,
    pub headers: Vec<NavSolutionHeader>,
    pub estimated_headers: StateFormatter<'a>,
    pub nominal_headers: StateFormatter<'a>,
}

impl<'a> NavSolutionFormatter<'a> {
    pub fn from_headers(headers: Vec<String>, filename: String, cosm: &'a Cosm) -> Self {
        let mut frames = HashMap::new();
        let mut hdrs = Vec::with_capacity(40);
        let mut est_hdrs = Vec::with_capacity(20);
        let mut nom_hdrs = Vec::with_capacity(20);
        // Rebuild the header tokens
        for hdr in &headers {
            let lowered = hdr.to_lowercase();
            let splt: Vec<&str> = lowered.split(':').collect();

            match splt[0] {
                "epoch" => {
                    hdrs.push(NavSolutionHeader::Epoch(if splt.len() == 2 {
                        EpochFormat::from_str(splt[1]).unwrap()
                    } else {
                        EpochFormat::GregorianUtc
                    }));
                }
                "delta_x" => hdrs.push(NavSolutionHeader::Delta_x),
                "delta_y" => hdrs.push(NavSolutionHeader::Delta_y),
                "delta_z" => hdrs.push(NavSolutionHeader::Delta_z),
                "delta_vx" => hdrs.push(NavSolutionHeader::Delta_vx),
                "delta_vy" => hdrs.push(NavSolutionHeader::Delta_vy),
                "delta_vz" => hdrs.push(NavSolutionHeader::Delta_vz),
                "cx_x" => hdrs.push(NavSolutionHeader::Cx_x),
                "cy_x" => hdrs.push(NavSolutionHeader::Cy_x),
                "cy_y" => hdrs.push(NavSolutionHeader::Cy_y),
                "cz_x" => hdrs.push(NavSolutionHeader::Cz_x),
                "cz_y" => hdrs.push(NavSolutionHeader::Cz_y),
                "cz_z" => hdrs.push(NavSolutionHeader::Cz_z),
                "cx_dot_x" => hdrs.push(NavSolutionHeader::Cx_dot_x),
                "cx_dot_y" => hdrs.push(NavSolutionHeader::Cx_dot_y),
                "cx_dot_z" => hdrs.push(NavSolutionHeader::Cx_dot_z),
                "cx_dot_x_dot" => hdrs.push(NavSolutionHeader::Cx_dot_x_dot),
                "cy_dot_x" => hdrs.push(NavSolutionHeader::Cy_dot_x),
                "cy_dot_y" => hdrs.push(NavSolutionHeader::Cy_dot_y),
                "cy_dot_z" => hdrs.push(NavSolutionHeader::Cy_dot_z),
                "cy_dot_x_dot" => hdrs.push(NavSolutionHeader::Cy_dot_x_dot),
                "cy_dot_y_dot" => hdrs.push(NavSolutionHeader::Cy_dot_y_dot),
                "cz_dot_x" => hdrs.push(NavSolutionHeader::Cz_dot_x),
                "cz_dot_y" => hdrs.push(NavSolutionHeader::Cz_dot_y),
                "cz_dot_z" => hdrs.push(NavSolutionHeader::Cz_dot_z),
                "cz_dot_x_dot" => hdrs.push(NavSolutionHeader::Cz_dot_x_dot),
                "cz_dot_y_dot" => hdrs.push(NavSolutionHeader::Cz_dot_y_dot),
                "cz_dot_z_dot" => hdrs.push(NavSolutionHeader::Cz_dot_z_dot),
                "estimate" | "nominal" => {
                    let frame_name = if splt.len() == 3 {
                        Some(splt[2].to_owned())
                    } else {
                        None
                    };
                    let state_hdr = match splt[1] {
                        "aol" => StateHeader::AoL { frame: frame_name },
                        "aop" => StateHeader::AoP { frame: frame_name },
                        "apoapsis" => StateHeader::apoapsis { frame: frame_name },
                        "ea" => StateHeader::EA { frame: frame_name },
                        "ecc" => StateHeader::ECC { frame: frame_name },
                        "energy" => StateHeader::energy { frame: frame_name },
                        "evec" => StateHeader::evec { frame: frame_name },
                        "geodetic_height" => StateHeader::geodetic_height { frame: frame_name },
                        "geodetic_latitude" => StateHeader::geodetic_latitude { frame: frame_name },
                        "geodetic_longitude" => {
                            StateHeader::geodetic_longitude { frame: frame_name }
                        }
                        "hmag" => StateHeader::hmag { frame: frame_name },
                        "hvec" => StateHeader::hvec { frame: frame_name },
                        "hx" => StateHeader::HX { frame: frame_name },
                        "hy" => StateHeader::HY { frame: frame_name },
                        "hz" => StateHeader::HZ { frame: frame_name },
                        "inc" => StateHeader::INC { frame: frame_name },
                        "ma" => StateHeader::MA { frame: frame_name },
                        "periapsis" => StateHeader::periapsis { frame: frame_name },
                        "period" => StateHeader::period { frame: frame_name },
                        "raan" => StateHeader::RAAN { frame: frame_name },
                        "radius" => StateHeader::radius { frame: frame_name },
                        "rmag" => StateHeader::rmag { frame: frame_name },
                        "semi_parameter" => StateHeader::semi_parameter { frame: frame_name },
                        "sma" => StateHeader::SMA { frame: frame_name },
                        "ta" => StateHeader::TA { frame: frame_name },
                        "tlong" => StateHeader::TLong { frame: frame_name },
                        "velocity" => StateHeader::velocity { frame: frame_name },
                        "vmag" => StateHeader::vmag { frame: frame_name },
                        "x" => StateHeader::X { frame: frame_name },
                        "y" => StateHeader::Y { frame: frame_name },
                        "z" => StateHeader::Z { frame: frame_name },
                        "vx" => StateHeader::VX { frame: frame_name },
                        "vy" => StateHeader::VY { frame: frame_name },
                        "vz" => StateHeader::VZ { frame: frame_name },
                        _ => panic!("unknown header `{}`", splt[0]),
                    };

                    if splt[0] == "estimate" {
                        est_hdrs.push(state_hdr);
                    } else {
                        nom_hdrs.push(state_hdr);
                    }
                }
                _ => panic!("unknown header `{}`", splt[0]),
            }

            if splt[0] != "epoch" && splt.len() == 2 {
                // Get the frame
                match cosm.try_frame(splt[1]) {
                    Ok(frame) => frames.insert(splt[1].to_string(), frame),
                    Err(e) => panic!("unknown frame `{}` in header ({})", splt[1], e),
                };
            }
        }

        // Add the nominal and estimate headers (needed to add the header row)
        hdrs.push(NavSolutionHeader::EstimatedState(est_hdrs.clone()));
        hdrs.push(NavSolutionHeader::NominalState(nom_hdrs.clone()));

        Self {
            filename,
            headers: hdrs,
            nominal_headers: StateFormatter {
                filename: "file_should_not_exist".to_owned(),
                headers: nom_hdrs,
                frames: HashMap::new(),
                cosm: &cosm,
            },
            estimated_headers: StateFormatter {
                filename: "file_should_not_exist".to_owned(),
                headers: est_hdrs,
                frames: HashMap::new(),
                cosm: &cosm,
            },
        }
    }

    /// Default headers are [Epoch (GregorianTai), X, Y, Z, VX, VY, VZ], where position is in km and velocity in km/s.
    pub fn default(filename: String, cosm: &'a Cosm) -> Self {
        let est_hdrs = vec![
            StateHeader::X { frame: None },
            StateHeader::Y { frame: None },
            StateHeader::Z { frame: None },
            StateHeader::VX { frame: None },
            StateHeader::VY { frame: None },
            StateHeader::VZ { frame: None },
        ];
        Self {
            filename,
            headers: vec![
                NavSolutionHeader::Epoch(EpochFormat::GregorianTai),
                NavSolutionHeader::Delta_x,
                NavSolutionHeader::Delta_y,
                NavSolutionHeader::Delta_z,
                NavSolutionHeader::Delta_vx,
                NavSolutionHeader::Delta_vy,
                NavSolutionHeader::Delta_vz,
                NavSolutionHeader::Cx_x,
                NavSolutionHeader::Cy_y,
                NavSolutionHeader::Cz_z,
                NavSolutionHeader::Cx_dot_x_dot,
                NavSolutionHeader::Cy_dot_y_dot,
                NavSolutionHeader::Cz_dot_z_dot,
                NavSolutionHeader::EstimatedState(est_hdrs.clone()),
            ],
            nominal_headers: StateFormatter {
                filename: "file_should_not_exist".to_owned(),
                headers: Vec::new(),
                frames: HashMap::new(),
                cosm: &cosm,
            },
            estimated_headers: StateFormatter {
                filename: "file_should_not_exist".to_owned(),
                headers: est_hdrs,
                frames: HashMap::new(),
                cosm: &cosm,
            },
        }
    }

    pub fn fmt<T: EstimableState<U6>, S: NavSolution<T>>(&self, sol: &S) -> Vec<String> {
        let mut formatted = Vec::new();

        for hdr in &self.headers {
            match hdr {
                NavSolutionHeader::EstimatedState(_) => {
                    // The formatter is already initialized
                    for fmtval in self.estimated_headers.fmt(&sol.orbital_state()) {
                        formatted.push(fmtval);
                    }
                }
                NavSolutionHeader::NominalState(_) => {
                    // The formatter is already initialized
                    for fmtval in self.nominal_headers.fmt(&sol.orbital_state()) {
                        formatted.push(fmtval);
                    }
                }
                NavSolutionHeader::Epoch(efmt) => formatted.push(efmt.format(sol.epoch())),
                NavSolutionHeader::Delta_x => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[0]))
                }
                NavSolutionHeader::Delta_y => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[1]))
                }
                NavSolutionHeader::Delta_z => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[2]))
                }
                NavSolutionHeader::Delta_vx => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[3]))
                }
                NavSolutionHeader::Delta_vy => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[4]))
                }
                NavSolutionHeader::Delta_vz => {
                    formatted.push(format!("{:.16e}", sol.state_deviation()[5]))
                }
                NavSolutionHeader::Cx_x => formatted.push(format!("{:.16e}", sol.covar()[(0, 0)])),
                NavSolutionHeader::Cy_x => formatted.push(format!("{:.16e}", sol.covar()[(1, 0)])),
                NavSolutionHeader::Cy_y => formatted.push(format!("{:.16e}", sol.covar()[(1, 1)])),
                NavSolutionHeader::Cz_x => formatted.push(format!("{:.16e}", sol.covar()[(2, 0)])),
                NavSolutionHeader::Cz_y => formatted.push(format!("{:.16e}", sol.covar()[(2, 1)])),
                NavSolutionHeader::Cz_z => formatted.push(format!("{:.16e}", sol.covar()[(2, 2)])),
                NavSolutionHeader::Cx_dot_x => {
                    formatted.push(format!("{:.16e}", sol.covar()[(3, 0)]))
                }
                NavSolutionHeader::Cx_dot_y => {
                    formatted.push(format!("{:.16e}", sol.covar()[(3, 1)]))
                }
                NavSolutionHeader::Cx_dot_z => {
                    formatted.push(format!("{:.16e}", sol.covar()[(3, 2)]))
                }
                NavSolutionHeader::Cx_dot_x_dot => {
                    formatted.push(format!("{:.16e}", sol.covar()[(3, 3)]))
                }
                NavSolutionHeader::Cy_dot_x => {
                    formatted.push(format!("{:.16e}", sol.covar()[(4, 0)]))
                }
                NavSolutionHeader::Cy_dot_y => {
                    formatted.push(format!("{:.16e}", sol.covar()[(4, 1)]))
                }
                NavSolutionHeader::Cy_dot_z => {
                    formatted.push(format!("{:.16e}", sol.covar()[(4, 2)]))
                }
                NavSolutionHeader::Cy_dot_x_dot => {
                    formatted.push(format!("{:.16e}", sol.covar()[(4, 3)]))
                }
                NavSolutionHeader::Cy_dot_y_dot => {
                    formatted.push(format!("{:.16e}", sol.covar()[(4, 4)]))
                }
                NavSolutionHeader::Cz_dot_x => {
                    formatted.push(format!("{:.16e}", sol.covar()[(5, 0)]))
                }
                NavSolutionHeader::Cz_dot_y => {
                    formatted.push(format!("{:.16e}", sol.covar()[(5, 1)]))
                }
                NavSolutionHeader::Cz_dot_z => {
                    formatted.push(format!("{:.16e}", sol.covar()[(5, 2)]))
                }
                NavSolutionHeader::Cz_dot_x_dot => {
                    formatted.push(format!("{:.16e}", sol.covar()[(5, 3)]))
                }
                NavSolutionHeader::Cz_dot_y_dot => {
                    formatted.push(format!("{:.16e}", sol.covar()[(5, 4)]))
                }
                NavSolutionHeader::Cz_dot_z_dot => {
                    formatted.push(format!("{:.16e}", sol.covar()[(5, 5)]))
                }
            };
        }

        formatted
    }
}
