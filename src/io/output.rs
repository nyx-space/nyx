use super::serde::{Serialize, Serializer};
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
        match &self.headers {
            Some(hdr) => StateFormatter::from_headers(hdr.to_vec(), self.filename.clone(), cosm),
            None => StateFormatter::default(self.filename.clone(), cosm),
        }
    }
}

/// Allowed headers, with an optional frame.
/// TODO: Support units
#[allow(non_camel_case_types)]
#[derive(Debug)]
pub enum Header {
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

impl Serialize for Header {
    /// NOTE: This is not part of unit testing because there is no deseralization of State (yet)
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self {
            Header::AoL { .. } => serializer.serialize_str("AoL"),
            Header::AoP { .. } => serializer.serialize_str("AoP"),
            Header::apoapsis { .. } => serializer.serialize_str("apoapsis"),
            Header::EA { .. } => serializer.serialize_str("EA"),
            Header::ECC { .. } => serializer.serialize_str("ECC"),
            Header::energy { .. } => serializer.serialize_str("energy"),
            Header::evec { .. } => serializer.serialize_str("evec"),
            Header::geodetic_height { .. } => serializer.serialize_str("geodetic_height"),
            Header::geodetic_latitude { .. } => serializer.serialize_str("geodetic_latitude"),
            Header::geodetic_longitude { .. } => serializer.serialize_str("geodetic_longitude"),
            Header::hmag { .. } => serializer.serialize_str("hmag"),
            Header::hvec { .. } => serializer.serialize_str("hvec"),
            Header::HX { .. } => serializer.serialize_str("HX"),
            Header::HY { .. } => serializer.serialize_str("HY"),
            Header::HZ { .. } => serializer.serialize_str("HZ"),
            Header::INC { .. } => serializer.serialize_str("INC"),
            Header::MA { .. } => serializer.serialize_str("MA"),
            Header::periapsis { .. } => serializer.serialize_str("periapsis"),
            Header::period { .. } => serializer.serialize_str("period"),
            Header::RAAN { .. } => serializer.serialize_str("RAAN"),
            Header::radius { .. } => serializer.serialize_str("radius"),
            Header::rmag { .. } => serializer.serialize_str("rmag"),
            Header::semi_parameter { .. } => serializer.serialize_str("semi_parameter"),
            Header::SMA { .. } => serializer.serialize_str("SMA"),
            Header::TA { .. } => serializer.serialize_str("TA"),
            Header::TLong { .. } => serializer.serialize_str("TLong"),
            Header::velocity { .. } => serializer.serialize_str("velocity"),
            Header::vmag { .. } => serializer.serialize_str("vmag"),
            Header::X { .. } => serializer.serialize_str("X"),
            Header::Y { .. } => serializer.serialize_str("Y"),
            Header::Z { .. } => serializer.serialize_str("Z"),
            Header::VX { .. } => serializer.serialize_str("VX"),
            Header::VY { .. } => serializer.serialize_str("VY"),
            Header::VZ { .. } => serializer.serialize_str("VZ"),
            Header::Epoch(efmt) => serializer.serialize_str(format!("Epoch:{:?}", efmt).as_str()),
        }
    }
}

pub struct StateFormatter<'a> {
    pub filename: String,
    pub headers: Vec<Header>,
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
                    hdrs.push(Header::Epoch(if splt.len() == 2 {
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
                        "aol" => Header::AoL { frame: frame_name },
                        "aop" => Header::AoP { frame: frame_name },
                        "apoapsis" => Header::apoapsis { frame: frame_name },
                        "ea" => Header::EA { frame: frame_name },
                        "ecc" => Header::ECC { frame: frame_name },
                        "energy" => Header::energy { frame: frame_name },
                        "evec" => Header::evec { frame: frame_name },
                        "geodetic_height" => Header::geodetic_height { frame: frame_name },
                        "geodetic_latitude" => Header::geodetic_latitude { frame: frame_name },
                        "geodetic_longitude" => Header::geodetic_longitude { frame: frame_name },
                        "hmag" => Header::hmag { frame: frame_name },
                        "hvec" => Header::hvec { frame: frame_name },
                        "hx" => Header::HX { frame: frame_name },
                        "hy" => Header::HY { frame: frame_name },
                        "hz" => Header::HZ { frame: frame_name },
                        "inc" => Header::INC { frame: frame_name },
                        "ma" => Header::MA { frame: frame_name },
                        "periapsis" => Header::periapsis { frame: frame_name },
                        "period" => Header::period { frame: frame_name },
                        "raan" => Header::RAAN { frame: frame_name },
                        "radius" => Header::radius { frame: frame_name },
                        "rmag" => Header::rmag { frame: frame_name },
                        "semi_parameter" => Header::semi_parameter { frame: frame_name },
                        "sma" => Header::SMA { frame: frame_name },
                        "ta" => Header::TA { frame: frame_name },
                        "tlong" => Header::TLong { frame: frame_name },
                        "velocity" => Header::velocity { frame: frame_name },
                        "vmag" => Header::vmag { frame: frame_name },
                        "x" => Header::X { frame: frame_name },
                        "y" => Header::Y { frame: frame_name },
                        "z" => Header::Z { frame: frame_name },
                        "vx" => Header::VX { frame: frame_name },
                        "vy" => Header::VY { frame: frame_name },
                        "vz" => Header::VZ { frame: frame_name },
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
                Header::Epoch(EpochFormat::GregorianTai),
                Header::X { frame: None },
                Header::Y { frame: None },
                Header::Z { frame: None },
                Header::VX { frame: None },
                Header::VY { frame: None },
                Header::VZ { frame: None },
            ],
            frames: HashMap::new(),
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
            match hdr {
                Header::Epoch(efmt) => formatted.push(efmt.format(state.dt)),
                Header::AoL { frame }
                | Header::AoP { frame }
                | Header::apoapsis { frame }
                | Header::EA { frame }
                | Header::ECC { frame }
                | Header::energy { frame }
                | Header::evec { frame }
                | Header::geodetic_height { frame }
                | Header::geodetic_latitude { frame }
                | Header::geodetic_longitude { frame }
                | Header::hmag { frame }
                | Header::hvec { frame }
                | Header::HX { frame }
                | Header::HY { frame }
                | Header::HZ { frame }
                | Header::INC { frame }
                | Header::MA { frame }
                | Header::periapsis { frame }
                | Header::period { frame }
                | Header::RAAN { frame }
                | Header::radius { frame }
                | Header::rmag { frame }
                | Header::semi_parameter { frame }
                | Header::SMA { frame }
                | Header::TA { frame }
                | Header::TLong { frame }
                | Header::velocity { frame }
                | Header::vmag { frame }
                | Header::X { frame }
                | Header::Y { frame }
                | Header::Z { frame }
                | Header::VX { frame }
                | Header::VY { frame }
                | Header::VZ { frame } => {
                    // Grab the state in the other frame if needed
                    let out_state = if frame.is_some() {
                        &mapped[&frame.as_ref().unwrap().to_lowercase()]
                    } else {
                        state
                    };

                    formatted.push(match hdr {
                        Header::AoL { .. } => format!("{:.9}", out_state.aol()),
                        Header::AoP { .. } => format!("{:.9}", out_state.aop()),
                        Header::apoapsis { .. } => format!("{:.9}", out_state.apoapsis()),
                        Header::EA { .. } => format!("{:.9}", out_state.ea()),
                        Header::ECC { .. } => format!("{:.9}", out_state.ecc()),
                        Header::energy { .. } => format!("{:.9}", out_state.energy()),
                        Header::evec { .. } => format!(
                            "[{:.9},{:.9},{:.9}]",
                            out_state.evec()[0],
                            out_state.evec()[1],
                            out_state.evec()[2]
                        ),
                        Header::geodetic_height { .. } => {
                            format!("{:.9}", out_state.geodetic_height())
                        }
                        Header::geodetic_latitude { .. } => {
                            format!("{:.9}", out_state.geodetic_latitude())
                        }
                        Header::geodetic_longitude { .. } => {
                            format!("{:.9}", out_state.geodetic_longitude())
                        }
                        Header::hmag { .. } => format!("{:.9}", out_state.hmag()),
                        Header::hvec { .. } => format!(
                            "[{:.9},{:.9},{:.9}]",
                            out_state.hvec()[0],
                            out_state.hvec()[1],
                            out_state.hvec()[2]
                        ),
                        Header::HX { .. } => format!("{:.9}", out_state.hx()),
                        Header::HY { .. } => format!("{:.9}", out_state.hy()),
                        Header::HZ { .. } => format!("{:.9}", out_state.hz()),
                        Header::INC { .. } => format!("{:.9}", out_state.inc()),
                        Header::MA { .. } => format!("{:.9}", out_state.ma()),
                        Header::periapsis { .. } => format!("{:.9}", out_state.periapsis()),
                        Header::period { .. } => format!("{:.9}", out_state.period()),
                        Header::RAAN { .. } => format!("{:.9}", out_state.raan()),
                        Header::radius { .. } => format!(
                            "[{:.9},{:.9},{:.9}]",
                            out_state.radius()[0],
                            out_state.radius()[1],
                            out_state.radius()[2]
                        ),
                        Header::rmag { .. } => format!("{:.9}", out_state.rmag()),
                        Header::semi_parameter { .. } => {
                            format!("{:.9}", out_state.semi_parameter())
                        }
                        Header::SMA { .. } => format!("{:.9}", out_state.sma()),
                        Header::TA { .. } => format!("{:.9}", out_state.ta()),
                        Header::TLong { .. } => format!("{:.9}", out_state.tlong()),
                        Header::velocity { .. } => format!(
                            "[{:.9},{:.9},{:.9}]",
                            out_state.velocity()[0],
                            out_state.velocity()[1],
                            out_state.velocity()[2]
                        ),
                        Header::vmag { .. } => format!("{:.9}", out_state.vmag()),
                        Header::X { .. } => format!("{:.9}", out_state.x),
                        Header::Y { .. } => format!("{:.9}", out_state.y),
                        Header::Z { .. } => format!("{:.9}", out_state.z),
                        Header::VX { .. } => format!("{:.9}", out_state.vx),
                        Header::VY { .. } => format!("{:.9}", out_state.vy),
                        Header::VZ { .. } => format!("{:.9}", out_state.vz),
                        _ => panic!("unsupported header `{:?}`", hdr),
                    });
                }
            };
        }

        formatted
    }
}
