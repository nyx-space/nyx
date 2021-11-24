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

extern crate flate2;
extern crate regex;
extern crate serde;
extern crate serde_derive;

use crate::errors::NyxError;
use crate::time::Epoch;
use std::convert::From;
use std::fmt;
use std::str::FromStr;

/// Handles loading of gravity models using files of NASA PDS and GMAT COF. Several gunzipped files are provided with nyx.
pub mod gravity;

/// Handles writing to an XYZV file
pub mod cosmo;

/// Handles reading from frames defined in input files
pub mod frame_serde;

/// Handles reading random variables
pub mod rv;

pub mod scenario;

pub mod odp;

pub mod formatter;

pub mod quantity;

pub mod ccsds;

#[allow(clippy::upper_case_acronyms)]
#[derive(Debug)]
pub enum ParsingError {
    MD(String),
    OD(String),
    UseOdInstead,
    UseMdInstead,
    EpochFormat,
    CovarFormat,
    FileNotFound(String),
    FileNotUTF8(String),
    SequenceNotFound(String),
    LoadingError(String),
    PropagatorNotFound(String),
    Duration(String),
    Quantity(String),
    Distance(String),
    Velocity(String),
    IllDefined(String),
    ExecutionError(NyxError),
}

impl From<NyxError> for ParsingError {
    fn from(error: NyxError) -> Self {
        Self::ExecutionError(error)
    }
}

/// Specifies the format of the Epoch during serialization
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EpochFormat {
    /// Default is GregorianUtc, as defined in [hifitime](https://docs.rs/hifitime/).
    GregorianUtc,
    GregorianTai,
    MjdTai,
    MjdTt,
    MjdUtc,
    JdeEt,
    JdeTai,
    JdeTt,
    JdeUtc,
    /// Seconds past a provided TAI Epoch
    TaiSecs(f64),
    /// Days past a provided TAI Epoch
    TaiDays(f64),
}

impl EpochFormat {
    pub fn format(&self, dt: Epoch) -> String {
        match *self {
            EpochFormat::GregorianUtc => dt.as_gregorian_utc_str(),
            EpochFormat::GregorianTai => dt.as_gregorian_tai_str(),
            EpochFormat::MjdTai => format!("{:.9}", dt.as_mjd_tai_days()),
            EpochFormat::MjdTt => format!("{:.9}", dt.as_mjd_tt_days()),
            EpochFormat::MjdUtc => format!("{:.9}", dt.as_mjd_utc_days()),
            EpochFormat::JdeEt => format!("{:.9}", dt.as_jde_et_days()),
            EpochFormat::JdeTai => format!("{:.9}", dt.as_jde_tai_days()),
            EpochFormat::JdeTt => format!("{:.9}", dt.as_jde_tt_days()),
            EpochFormat::JdeUtc => format!("{:.9}", dt.as_jde_utc_days()),
            EpochFormat::TaiSecs(e) => format!("{:.9}", dt.as_tai_seconds() - e),
            EpochFormat::TaiDays(e) => format!("{:.9}", dt.as_tai_days() - e),
        }
    }
}

impl fmt::Display for EpochFormat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            EpochFormat::GregorianUtc => write!(f, "Gregorian UTC"),
            EpochFormat::GregorianTai => write!(f, "Gregorian TAI"),
            EpochFormat::MjdTai => write!(f, "MJD TAI"),
            EpochFormat::MjdTt => write!(f, "MJD TT"),
            EpochFormat::MjdUtc => write!(f, "MJD UTC"),
            EpochFormat::JdeEt => write!(f, "JDE ET"),
            EpochFormat::JdeTai => write!(f, "JDE TAI"),
            EpochFormat::JdeTt => write!(f, "JDE TT"),
            EpochFormat::JdeUtc => write!(f, "JDE UTC"),
            EpochFormat::TaiSecs(_) => write!(f, "TAI+ s"),
            EpochFormat::TaiDays(_) => write!(f, "TAI+ days"),
        }
    }
}

impl FromStr for EpochFormat {
    type Err = ParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().replace(" ", "").as_str() {
            "gregorianutc" => Ok(EpochFormat::GregorianUtc),
            "gregoriantai" => Ok(EpochFormat::GregorianTai),
            "mjdtai" => Ok(EpochFormat::MjdTai),
            "mjdtt" => Ok(EpochFormat::MjdTt),
            "mjdutc" => Ok(EpochFormat::MjdUtc),
            "jdeet" => Ok(EpochFormat::JdeEt),
            "jdetai" => Ok(EpochFormat::JdeTai),
            "jdett" => Ok(EpochFormat::JdeTt),
            "jdeutc" => Ok(EpochFormat::JdeUtc),
            _ => Err(ParsingError::EpochFormat),
        }
    }
}

/// Specifies the format of the covariance during serialization
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum CovarFormat {
    /// Default: allows plotting the variance of the elements instead of the covariance
    Sqrt,
    /// Keeps the covariance as computed, i.e. one sigma (~68%), causes e.g. positional elements in km^2.
    Sigma1,
    /// Three sigma covers about 99.7% of the distribution
    Sigma3,
    /// Allows specifying a custom multiplication factor of each element of the covariance.
    MulSigma(f64),
}

impl fmt::Display for CovarFormat {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            CovarFormat::Sqrt => write!(f, "exptd_val_"),
            CovarFormat::Sigma1 => write!(f, "covar_"),
            CovarFormat::Sigma3 => write!(f, "3sig_covar"),
            CovarFormat::MulSigma(x) => write!(f, "{}sig_covar", x),
        }
    }
}

impl FromStr for CovarFormat {
    type Err = ParsingError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().replace(" ", "").as_str() {
            "sqrt" => Ok(CovarFormat::Sqrt),
            "1sigma" | "sigma1" => Ok(CovarFormat::Sigma1),
            "3sigma" | "sigma3" => Ok(CovarFormat::Sigma3),
            _ => Err(ParsingError::CovarFormat),
        }
    }
}
