use crate::NyxError;
use std::convert::TryFrom;

/// Defines the default celestial bodies in the provided de438 XB.
#[derive(Copy, Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub enum Bodies {
    SSB,
    Sun,
    MercuryBarycenter,
    Mercury,
    VenusBarycenter,
    Venus,
    EarthBarycenter,
    Earth,
    Luna,
    MarsBarycenter,
    JupiterBarycenter,
    SaturnBarycenter,
    UranusBarycenter,
    NeptuneBarycenter,
    PlutoBarycenter,
}

impl Bodies {
    /// Returns the path in the standard de438 XB
    pub fn ephem_path(&self) -> &'static [usize] {
        match *self {
            Self::SSB => &[],
            Self::Sun => &[0],
            Self::MercuryBarycenter => &[1],
            Self::Mercury => &[1],
            Self::VenusBarycenter => &[2],
            Self::Venus => &[2],
            Self::EarthBarycenter => &[3],
            Self::Earth => &[3, 0],
            Self::Luna => &[3, 1],
            Self::MarsBarycenter => &[4],
            Self::JupiterBarycenter => &[5],
            Self::SaturnBarycenter => &[6],
            Self::UranusBarycenter => &[7],
            Self::NeptuneBarycenter => &[8],
            Self::PlutoBarycenter => &[9],
        }
    }

    /// Returns the human name
    pub fn name(&self) -> String {
        match *self {
            Self::SSB => "Solar System Barycenter".to_string(),
            Self::Sun => "Sun".to_string(),
            Self::MercuryBarycenter => "Mercury".to_string(),
            Self::Mercury => "Mercury".to_string(),
            Self::VenusBarycenter => "Venus".to_string(),
            Self::Venus => "Venus".to_string(),
            Self::EarthBarycenter => "Earth Moon Barycenter".to_string(),
            Self::Earth => "Earth".to_string(),
            Self::Luna => "Moon".to_string(),
            Self::MarsBarycenter => "Mars".to_string(),
            Self::JupiterBarycenter => "Jupiter Barycenter".to_string(),
            Self::SaturnBarycenter => "Saturn Barycenter".to_string(),
            Self::UranusBarycenter => "Uranus Barycenter".to_string(),
            Self::NeptuneBarycenter => "Neptune Barycenter".to_string(),
            Self::PlutoBarycenter => "Pluto Barycenter".to_string(),
        }
    }
}

impl TryFrom<String> for Bodies {
    type Error = NyxError;

    fn try_from(name: String) -> Result<Self, Self::Error> {
        match name.to_lowercase().as_str() {
            "solar system barycenter" | "ssb" => Ok(Self::SSB),
            "sun" => Ok(Self::Sun),
            "mercury" => Ok(Self::Mercury),
            "venus" => Ok(Self::Venus),
            "earth moon barycenter" => Ok(Self::EarthBarycenter),
            "earth" => Ok(Self::Earth),
            "moon" | "luna" => Ok(Self::Luna),
            "mars" | "mars barycenter" => Ok(Self::MarsBarycenter),
            "jupiter" | "jupiter barycenter" => Ok(Self::JupiterBarycenter),
            "saturn" | "saturn barycenter" => Ok(Self::SaturnBarycenter),
            "uranus" | "uranus barycenter" => Ok(Self::UranusBarycenter),
            "neptune" | "neptune barycenter" => Ok(Self::NeptuneBarycenter),
            "pluto" | "pluto barycenter" => Ok(Self::PlutoBarycenter),
            _ => Err(NyxError::ObjectNotFound(name)),
        }
    }
}

impl TryFrom<Vec<usize>> for Bodies {
    type Error = NyxError;

    fn try_from(ephem_path: Vec<usize>) -> Result<Self, Self::Error> {
        match ephem_path.len() {
            0 => Ok(Self::SSB),
            1 => match ephem_path[0] {
                0 => Ok(Self::Sun),
                1 => Ok(Self::Mercury),
                2 => Ok(Self::Venus),
                3 => Ok(Self::EarthBarycenter),
                4 => Ok(Self::MarsBarycenter),
                5 => Ok(Self::JupiterBarycenter),
                6 => Ok(Self::SaturnBarycenter),
                7 => Ok(Self::UranusBarycenter),
                8 => Ok(Self::NeptuneBarycenter),
                9 => Ok(Self::PlutoBarycenter),
                _ => Err(NyxError::ObjectNotFound(format!("{:?}", ephem_path))),
            },
            2 if ephem_path[0] == 3 => match ephem_path[1] {
                // This only support the Earth system
                0 => Ok(Self::Earth),
                1 => Ok(Self::Luna),
                _ => Err(NyxError::ObjectNotFound(format!("{:?}", ephem_path))),
            },
            _ => Err(NyxError::ObjectNotFound(format!("{:?}", ephem_path))),
        }
    }
}
