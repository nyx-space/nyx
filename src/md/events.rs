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

use crate::celestia::{Cosm, Frame, Orbit};
use crate::dimensions::allocator::Allocator;
use crate::dimensions::DefaultAllocator;
use crate::errors::NyxError;
use crate::time::{Duration, TimeUnit};
use crate::utils::between_pm_x;
use crate::SpacecraftState;
use crate::State;
use std::fmt;
use std::str::FromStr;
use std::sync::Arc;

fn angled_value(cur_angle: f64, desired_angle: f64) -> f64 {
    if between_pm_x(cur_angle, desired_angle) > 0.0 {
        cur_angle - desired_angle
    } else {
        cur_angle + 2.0 * desired_angle
    }
}

/// A trait to specify how a specific event must be evaluated.
pub trait EventEvaluator<S: State>: fmt::Display + Send + Sync
where
    DefaultAllocator: Allocator<f64, S::Size>,
{
    // Evaluation of event crossing, must return whether the condition happened between between both states.
    fn eval_crossing(&self, prev_state: &S, next_state: &S) -> bool {
        self.eval(prev_state) * self.eval(next_state) < 0.0
    }

    // Evaluation of the event, must return a value corresponding to whether the state is before or after the event
    fn eval(&self, state: &S) -> f64;

    fn epoch_precision(&self) -> Duration;
    fn value_precision(&self) -> f64;
}

#[derive(Clone, Debug)]
pub struct Event {
    pub parameter: StateParameter,
    pub epoch_precision: TimeUnit,
    pub value_precision: f64,
    pub in_frame: Option<(Frame, Arc<Cosm>)>,
}

impl fmt::Display for Event {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.parameter)?;
        if let Some((frame, _)) = self.in_frame {
            write!(f, "in frame {}", frame)?;
        }
        fmt::Result::Ok(())
    }
}

impl Event {
    /// Match a specific event for the parameter to hit the specified value.
    /// By default, the time precision is 1 millisecond and the value precision is 1e-3 of whatever
    /// unit is the default for that parameter. For example, a radius event will seek the requested
    /// value at the meter level, and an angle event will seek it at the thousands of a degree.
    pub fn new(parameter: StateParameter) -> Self {
        Self {
            parameter,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 1e-3,
            in_frame: None,
        }
    }

    /// Allows only reporting the current value
    pub fn status(parameter: StateParameter) -> Self {
        Self {
            parameter,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 0.0,
            in_frame: None,
        }
    }

    /// Match a specific event in another frame, using the default epoch precision and value.
    pub fn in_frame(parameter: StateParameter, target_frame: Frame, cosm: Arc<Cosm>) -> Self {
        Self {
            parameter,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 1e-3,
            in_frame: Some((target_frame, cosm)),
        }
    }

    /// Report the status in a specific frame
    pub fn status_in_frame(
        parameter: StateParameter,
        target_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self {
            parameter,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 0.0,
            in_frame: Some((target_frame, cosm)),
        }
    }
}

/// Allowed headers, with an optional frame.
#[allow(non_camel_case_types)]
#[derive(Clone, Debug, PartialEq)]
pub enum StateParameter {
    /// Argument of Latitude (deg)
    AoL(f64),
    /// Argument of Periapse (deg)
    AoP(f64),
    /// Apoapsis, shortcut for TA == 180.0
    Apoapsis,
    /// Radius of apoapsis (km)
    ApoapsisRadius(f64),
    /// Eccentric anomaly (deg)
    EccentricAnomaly(f64),
    /// Eccentricity (no unit)
    Eccentricity(f64),
    /// Specific energy
    Energy(f64),
    /// fuel mass in kilograms
    FuelMass(f64),
    /// Geodetic height (km)
    GeodeticHeight(f64),
    /// Geodetic latitude (deg)
    GeodeticLatitude(f64),
    /// Geodetic longitude (deg)
    GeodeticLongitude(f64),
    /// Orbital momentum
    Hmag(f64),
    /// X component of the orbital momentum vector
    HX(f64),
    /// Y component of the orbital momentum vector
    HY(f64),
    /// Z component of the orbital momentum vector
    HZ(f64),
    /// Inclination (deg)
    Inclination(f64),
    /// Mean anomaly (deg)
    MeanAnomaly(f64),
    /// Periapsis, shortcut for TA == 0.0
    Periapsis,
    /// Radius of periapse (km)
    PeriapsisRadius(f64),
    /// Orbital period (s)
    Period(f64),
    /// Right ascension of the ascending node (deg)
    RAAN(f64),
    /// Norm of the radius vector
    Rmag(f64),
    /// Semi parameter (km)
    SemiParameter(f64),
    /// Semi major axis (km)
    SMA(f64),
    /// True anomaly
    TrueAnomaly(f64),
    /// True longitude
    TrueLongitude(f64),
    /// Norm of the velocity vector (km/s)
    Vmag(f64),
    /// X component of the radius (km)
    X(f64),
    /// Y component of the radius (km)
    Y(f64),
    /// Z component of the radius (km)
    Z(f64),
    /// X component of the velocity (km/s)
    VX(f64),
    /// Y component of the velocity (km/s)
    VY(f64),
    /// Z component of the velocity (km/s)
    VZ(f64),
    /// Allows creating new custom events by matching on the value of the field
    Custom(usize, f64),
}

impl FromStr for StateParameter {
    type Err = NyxError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let rplt = s.to_lowercase().replace("=", "");
        let parts: Vec<&str> = rplt.split(' ').collect();
        if parts.len() == 1 {
            match parts[0].trim() {
                "apoapsis" => Ok(Self::Apoapsis),
                "periapsis" => Ok(Self::Periapsis),
                _ => Err(NyxError::LoadingError(format!(
                    "Unknown event state parameter without value: {}",
                    s
                ))),
            }
        } else if parts.len() == 2 {
            // let value: f64 = parts[1].trim().parse().unwrap();
            match parts[1].trim().parse::<f64>() {
                Ok(value) => match parts[0].trim() {
                    "aol" => Ok(Self::AoL(value)),
                    "aop" => Ok(Self::AoP(value)),
                    "apoapsisradius" => Ok(Self::ApoapsisRadius(value)),
                    "eccentricanomaly" => Ok(Self::EccentricAnomaly(value)),
                    "eccentricity" => Ok(Self::Eccentricity(value)),
                    "energy" => Ok(Self::Energy(value)),
                    "fuelmass" => Ok(Self::FuelMass(value)),
                    "geodeticheight" => Ok(Self::GeodeticHeight(value)),
                    "geodeticlatitude" => Ok(Self::GeodeticLatitude(value)),
                    "geodeticlongitude" => Ok(Self::GeodeticLongitude(value)),
                    "hmag" => Ok(Self::Hmag(value)),
                    "hx" => Ok(Self::HX(value)),
                    "hy" => Ok(Self::HY(value)),
                    "hz" => Ok(Self::HZ(value)),
                    "inclination" => Ok(Self::Inclination(value)),
                    "meananomaly" => Ok(Self::MeanAnomaly(value)),
                    "periapsisradius" => Ok(Self::PeriapsisRadius(value)),
                    "period" => Ok(Self::Period(value)),
                    "raan" => Ok(Self::RAAN(value)),
                    "rmag" => Ok(Self::Rmag(value)),
                    "semiparameter" => Ok(Self::SemiParameter(value)),
                    "sma" => Ok(Self::SMA(value)),
                    "trueanomaly" => Ok(Self::TrueAnomaly(value)),
                    "truelongitude" => Ok(Self::TrueLongitude(value)),
                    "vmag" => Ok(Self::Vmag(value)),
                    "x" => Ok(Self::X(value)),
                    "y" => Ok(Self::Y(value)),
                    "z" => Ok(Self::Z(value)),
                    "vx" => Ok(Self::VX(value)),
                    "vy" => Ok(Self::VY(value)),
                    "vz" => Ok(Self::VZ(value)),
                    _ => Err(NyxError::LoadingError(format!(
                        "Unknown event state parameter: {}",
                        s
                    ))),
                },
                Err(_) => Err(NyxError::LoadingError(format!(
                    "Could not parse value in event: {}",
                    s
                ))),
            }
        } else {
            Err(NyxError::LoadingError(format!(
                "Could not parse event: {}",
                s
            )))
        }
    }
}

impl EventEvaluator<Orbit> for Event {
    #[allow(clippy::identity_op)]
    fn epoch_precision(&self) -> Duration {
        1 * self.epoch_precision
    }

    fn value_precision(&self) -> f64 {
        self.value_precision
    }

    fn eval(&self, state: &Orbit) -> f64 {
        // Transform the state if needed
        let state = if let Some((frame, cosm)) = &self.in_frame {
            cosm.frame_chg(state, *frame)
        } else {
            *state
        };

        // Return the parameter centered around the desired value
        match self.parameter {
            StateParameter::AoL(value) => angled_value(state.aol(), value),
            StateParameter::AoP(value) => angled_value(state.aop(), value),
            StateParameter::Apoapsis => angled_value(state.ta(), 180.0),
            StateParameter::ApoapsisRadius(value) => state.apoapsis() - value,
            StateParameter::EccentricAnomaly(value) => angled_value(state.ea(), value),
            StateParameter::Eccentricity(value) => state.ecc() - value,
            StateParameter::Energy(value) => state.energy() - value,
            StateParameter::GeodeticHeight(value) => state.geodetic_height() - value,
            StateParameter::GeodeticLatitude(value) => state.geodetic_latitude() - value,
            StateParameter::GeodeticLongitude(value) => state.geodetic_longitude() - value,
            StateParameter::Hmag(value) => state.hmag() - value,
            StateParameter::HX(value) => state.hx() - value,
            StateParameter::HY(value) => state.hy() - value,
            StateParameter::HZ(value) => state.hz() - value,
            StateParameter::Inclination(value) => angled_value(state.inc(), value),
            StateParameter::MeanAnomaly(value) => angled_value(state.ma(), value),
            StateParameter::Periapsis => angled_value(state.ta(), 0.0),
            StateParameter::PeriapsisRadius(value) => state.periapsis() - value,
            StateParameter::Period(value) => state.period().in_seconds() - value,
            StateParameter::RAAN(value) => angled_value(state.raan(), value),
            StateParameter::Rmag(value) => state.rmag() - value,
            StateParameter::SemiParameter(value) => state.semi_parameter() - value,
            StateParameter::SMA(value) => state.sma() - value,
            StateParameter::TrueAnomaly(value) => angled_value(state.ta(), value),
            StateParameter::TrueLongitude(value) => angled_value(state.tlong(), value),
            StateParameter::Vmag(value) => state.vmag() - value,
            StateParameter::X(value) => state.x - value,
            StateParameter::Y(value) => state.y - value,
            StateParameter::Z(value) => state.z - value,
            StateParameter::VX(value) => state.vx - value,
            StateParameter::VY(value) => state.vy - value,
            StateParameter::VZ(value) => state.vz - value,
            _ => unimplemented!(),
        }
    }
}

impl EventEvaluator<SpacecraftState> for Event {
    fn eval(&self, state: &SpacecraftState) -> f64 {
        match self.parameter {
            StateParameter::FuelMass(value) => state.fuel_mass_kg - value,
            _ => self.eval(&state.orbit),
        }
    }

    #[allow(clippy::identity_op)]
    fn epoch_precision(&self) -> Duration {
        1 * self.epoch_precision
    }

    fn value_precision(&self) -> f64 {
        self.value_precision
    }
}

// /// Built-in events, will likely be expanded as development continues.
// #[derive(Clone, Copy, Debug)]
// pub enum EventKind {
//     Sma(f64),
//     Ecc(f64),
//     Inc(f64),
//     Raan(f64),
//     Aop(f64),
//     TA(f64),
//     Periapse,
//     Apoapse,
//     Fuel(f64),
// }

// impl fmt::Display for EventKind {
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         write!(f, "{:?}", self)
//     }
// }

// /// An orbital event, in the same frame or in another frame.
// #[derive(Debug)]
// pub struct OrbitalEvent {
//     pub kind: EventKind,
//     pub tgt: Option<Frame>,
//     pub cosm: Option<Arc<Cosm>>,
// }

// impl OrbitalEvent {
//     pub fn new(kind: EventKind) -> Box<Self> {
//         Box::new(OrbitalEvent {
//             kind,
//             tgt: None,
//             cosm: None,
//         })
//     }
//     pub fn in_frame(kind: EventKind, tgt: Frame, cosm: Arc<Cosm>) -> Box<Self> {
//         Box::new(OrbitalEvent {
//             kind,
//             tgt: Some(tgt),
//             cosm: Some(cosm),
//         })
//     }
// }

// impl Event<Orbit> for OrbitalEvent {
//     fn eval(&self, state: &Orbit) -> f64 {
//         let state = match self.tgt {
//             Some(tgt) => self.cosm.unwrap().frame_chg(state, tgt),
//             None => *state,
//         };

//         match self.kind {
//             EventKind::Sma(sma) => state.sma() - sma,
//             EventKind::Ecc(ecc) => state.ecc() - ecc,
//             EventKind::Inc(inc) => state.inc() - inc,
//             EventKind::Raan(raan) => state.raan() - raan,
//             EventKind::Aop(aop) => state.aop() - aop,
//             EventKind::TA(angle) => {
//                 if between_pm_x(state.ta(), angle) > 0.0 {
//                     state.ta() - angle
//                 } else {
//                     state.ta() + 2.0 * angle
//                 }
//             }
//             EventKind::Periapse => between_pm_180(state.ta()),
//             EventKind::Apoapse => {
//                 if between_pm_180(state.ta()) > 0.0 {
//                     state.ta() - 180.0
//                 } else {
//                     state.ta() + 360.0
//                 }
//             }
//             _ => panic!("event {:?} not supported", self.kind),
//         }
//     }

//     fn eval_crossing(&self, prev_state: &Orbit, next_state: &Orbit) -> bool {
//         let prev_val = self.eval(prev_state);
//         let next_val = self.eval(next_state);
//         match self.kind {
//             // XXX: Should this condition be applied to all angles?
//             EventKind::Periapse => prev_val < 0.0 && next_val >= 0.0,
//             EventKind::Apoapse => prev_val > 0.0 && next_val <= 0.0,
//             _ => prev_val * next_val <= 0.0,
//         }
//     }
// }

// #[derive(Debug)]
// pub struct SCEvent {
//     pub kind: EventKind,
//     pub orbital: Option<Box<OrbitalEvent>>,
// }

// impl SCEvent {
//     pub fn fuel_mass(mass: f64) -> Box<Self> {
//         Box::new(Self {
//             kind: EventKind::Fuel(mass),
//             orbital: None,
//         })
//     }
//     pub fn orbital(event: Box<OrbitalEvent>) -> Box<Self> {
//         Box::new(Self {
//             kind: event.kind,
//             orbital: Some(event),
//         })
//     }
// }

// impl Event<SpacecraftState> for SCEvent {
//     fn eval(&self, state: &SpacecraftState) -> f64 {
//         match self.kind {
//             EventKind::Fuel(mass) => state.fuel_mass_kg - mass,
//             _ => self.orbital.as_ref().unwrap().eval(&state.orbit),
//         }
//     }

//     fn eval_crossing(&self, prev_state: &SpacecraftState, next_state: &SpacecraftState) -> bool {
//         match self.kind {
//             EventKind::Fuel(mass) => {
//                 prev_state.fuel_mass_kg <= mass && next_state.fuel_mass_kg > mass
//             }
//             _ => self
//                 .orbital
//                 .as_ref()
//                 .unwrap()
//                 .eval_crossing(&prev_state.orbit, &next_state.orbit),
//         }
//     }
// }

// #[test]
// fn demo() {
//     let event = Event::new(StateParameter::TA, value: f64);
//     traj.find_minmax(StateParameter::TA, TimeUnit::Second);
// }
