/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

use super::regex::Regex;
use crate::time::{SECONDS_PER_DAY, SECONDS_PER_HOUR, SECONDS_PER_MINUTE};

pub use super::ParsingError;

#[derive(Copy, Clone, Debug)]
pub enum Quantity {
    /// Stores the distance in kilometers
    Distance(f64),
    /// Stores the velocity in km/s
    Velocity(f64),
    /// Stores the duration in seconds
    Duration(f64),
}

impl Quantity {
    /// Returns the value of this quantity in kilometers for distances and kilometers per second for velocities.
    pub fn v(self) -> f64 {
        match self {
            Self::Distance(v) => v,
            Self::Velocity(v) => v,
            Self::Duration(v) => v,
        }
    }
}

/// Parse a duration
///
/// ```
/// extern crate nyx_space as nyx;
///
/// use nyx::io::quantity::parse_duration;
/// use std::f64::EPSILON;
///
/// assert!((parse_duration("1 * days").unwrap().v() - 86_400.0).abs() < EPSILON);
/// assert!((parse_duration("1 days").unwrap().v() - 86_400.0).abs() < EPSILON);
/// assert!((parse_duration("1* day").unwrap().v() - 86_400.0).abs() < EPSILON);
/// assert!((parse_duration("1 *d").unwrap().v() - 86_400.0).abs() < EPSILON);
/// assert!((parse_duration("1 * h").unwrap().v() - 3_600.0).abs() < EPSILON);
/// assert!((parse_duration("1 h").unwrap().v() - 3_600.0).abs() < EPSILON);
/// assert!((parse_duration("1.0000 * hour").unwrap().v() - 3_600.0).abs() < EPSILON);
/// assert!((parse_duration("1.0000 * hours").unwrap().v() - 3_600.0).abs() < EPSILON);
/// assert!((parse_duration("1.0 * min").unwrap().v() - 60.0).abs() < EPSILON);
/// assert!((parse_duration("1. * s").unwrap().v() - 1.0).abs() < EPSILON);
/// assert!((parse_duration("1 * s").unwrap().v() - 1.0).abs() < EPSILON);
/// assert!((parse_duration("1 s").unwrap().v() - 1.0).abs() < EPSILON);
/// ```
pub fn parse_duration(duration: &str) -> Result<Quantity, ParsingError> {
    let reg = Regex::new(r"^(\d+\.?\d*)\W*(\w+)$").unwrap();
    match reg.captures(duration) {
        Some(cap) => {
            let mut time_s = cap[1].to_owned().parse::<f64>().unwrap();
            match cap[2].to_owned().to_lowercase().as_str() {
                "days" | "day" | "d" => time_s *= SECONDS_PER_DAY,
                "hours" | "hour" | "h" => time_s *= SECONDS_PER_HOUR,
                "min" | "mins" | "minute" | "minutes" | "m" => time_s *= SECONDS_PER_MINUTE,
                "s" | "sec" | "secs" => time_s *= 1.0,
                _ => {
                    return Err(ParsingError::Duration(format!(
                        "unknown duration unit in `{}`",
                        duration
                    )))
                }
            }
            Ok(Quantity::Duration(time_s))
        }
        None => Err(ParsingError::Duration(format!(
            "Could not parse stopping condition: `{}`",
            duration
        ))),
    }
}

/// Parse a distance or velocity
///
/// ```
/// extern crate nyx_space as nyx;
/// use nyx::io::quantity::parse_quantity;
/// use std::f64::EPSILON;
///
/// assert!((parse_quantity("1.0 km").unwrap().v() - 1.0).abs() < EPSILON);
/// assert!((parse_quantity("-1.3 mm").unwrap().v() - -1.3e-6).abs() < EPSILON);
/// assert!((parse_quantity("3.4e3 m/s").unwrap().v() - 3.4).abs() < EPSILON);
/// assert!((parse_quantity("3.4e0 km/s").unwrap().v() - 3.4).abs() < EPSILON);
/// assert!((parse_quantity("3.4e-3 Mm/s").unwrap().v() - 3.4).abs() < EPSILON);
/// assert!((parse_quantity("7 m").unwrap().v() - 7e-3).abs() < EPSILON);
/// assert!((parse_quantity("-6 km/h").unwrap().v() - -6.0/3_600.0).abs() < EPSILON);
/// ```
pub fn parse_quantity(input: &str) -> Result<Quantity, ParsingError> {
    let reg =
        Regex::new(r#"(-?\d+\.?\d*(?:e-?\d+\.?\d*)?)\W*([G|M|k|m|u|n]?)m/?([h|s])?"#).unwrap();

    match reg.captures(input) {
        Some(cap) => {
            let mut value = cap[1].to_owned().parse::<f64>().unwrap();
            // The second group can be empty, in which case the input was in meters.
            match cap[2].to_owned().as_str() {
                "G" => value *= 1e6,
                "M" => value *= 1e3,
                "k" => value *= 1.0,
                "" => value *= 1e-3,
                "m" => value *= 1e-6,
                "u" => value *= 1e-9,
                "n" => value *= 1e-12,
                _ => {
                    return Err(ParsingError::Distance(format!(
                        "unknown distance multiplier unit in `{}`",
                        input
                    )))
                }
            }
            if let Some(time_div) = cap.get(3) {
                // This is a velocity
                match time_div.as_str().to_lowercase().as_str() {
                    "h" => value /= SECONDS_PER_HOUR,
                    "s" => value *= 1.0,
                    _ => {
                        return Err(ParsingError::Velocity(format!(
                            "unknown time divisor unit in `{}`",
                            input
                        )))
                    }
                }
                Ok(Quantity::Velocity(value))
            } else {
                Ok(Quantity::Distance(value))
            }
        }
        None => {
            // Try to parse as a duration
            match parse_duration(input) {
                Ok(v) => Ok(v),
                Err(_) => Err(ParsingError::Quantity(format!(
                    "Could not understand quantity: `{}`",
                    input
                ))),
            }
        }
    }
}
