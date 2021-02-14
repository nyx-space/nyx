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
    pub value: f64,
    pub epoch_precision: TimeUnit,
    pub value_precision: f64,
    pub in_frame: Option<(Frame, Arc<Cosm>)>,
}

impl fmt::Display for Event {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self.parameter)?;
        if self.value.abs() > std::f64::EPSILON {
            write!(f, "within value of {:.6}", self.value)?;
        }
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
    pub fn new(parameter: StateParameter, value: f64) -> Self {
        Self {
            parameter,
            value,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 1e-3,
            in_frame: None,
        }
    }

    /// Allows only reporting the current value
    pub fn status(parameter: StateParameter) -> Self {
        Self {
            parameter,
            value: 0.0,
            epoch_precision: TimeUnit::Millisecond,
            value_precision: 0.0,
            in_frame: None,
        }
    }

    /// Match a specific event in another frame, using the default epoch precision and value.
    pub fn in_frame(
        parameter: StateParameter,
        value: f64,
        target_frame: Frame,
        cosm: Arc<Cosm>,
    ) -> Self {
        Self {
            parameter,
            value,
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
            value: 0.0,
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
    AoL,
    /// Argument of Periapse (deg)
    AoP,
    /// Apoapsis, shortcut for TA == 180.0
    Apoapsis,
    /// Radius of apoapsis (km)
    ApoapsisRadius,
    /// Eccentric anomaly (deg)
    EccentricAnomaly,
    /// Eccentricity (no unit)
    Eccentricity,
    /// Specific energy
    Energy,
    /// fuel mass in kilograms
    FuelMass,
    /// Geodetic height (km)
    GeodeticHeight,
    /// Geodetic latitude (deg)
    GeodeticLatitude,
    /// Geodetic longitude (deg)
    GeodeticLongitude,
    /// Orbital momentum
    Hmag,
    /// X component of the orbital momentum vector
    HX,
    /// Y component of the orbital momentum vector
    HY,
    /// Z component of the orbital momentum vector
    HZ,
    /// Inclination (deg)
    Inclination,
    /// Mean anomaly (deg)
    MeanAnomaly,
    /// Periapsis, shortcut for TA == 0.0
    Periapsis,
    /// Radius of periapse (km)
    PeriapsisRadius,
    /// Orbital period (s)
    Period,
    /// Right ascension of the ascending node (deg)
    RAAN,
    /// Norm of the radius vector
    Rmag,
    /// Semi parameter (km)
    SemiParameter,
    /// Semi major axis (km)
    SMA,
    /// True anomaly
    TrueAnomaly,
    /// True longitude
    TrueLongitude,
    /// Norm of the velocity vector (km/s)
    Vmag,
    /// X component of the radius (km)
    X,
    /// Y component of the radius (km)
    Y,
    /// Z component of the radius (km)
    Z,
    /// X component of the velocity (km/s)
    VX,
    /// Y component of the velocity (km/s)
    VY,
    /// Z component of the velocity (km/s)
    VZ,
    /// Allows creating new custom events by matching on the value of the field
    Custom(usize),
}

impl FromStr for StateParameter {
    type Err = NyxError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s.to_lowercase().as_str() {
            "aol" => Ok(Self::AoL),
            "aop" => Ok(Self::AoP),
            "apoapsis" => Ok(Self::Apoapsis),
            "apoapsisradius" => Ok(Self::ApoapsisRadius),
            "eccentricanomaly" => Ok(Self::EccentricAnomaly),
            "eccentricity" => Ok(Self::Eccentricity),
            "energy" => Ok(Self::Energy),
            "fuelmass" => Ok(Self::FuelMass),
            "geodeticheight" => Ok(Self::GeodeticHeight),
            "geodeticlatitude" => Ok(Self::GeodeticLatitude),
            "geodeticlongitude" => Ok(Self::GeodeticLongitude),
            "hmag" => Ok(Self::Hmag),
            "hx" => Ok(Self::HX),
            "hy" => Ok(Self::HY),
            "hz" => Ok(Self::HZ),
            "inclination" => Ok(Self::Inclination),
            "meananomaly" => Ok(Self::MeanAnomaly),
            "periapsis" => Ok(Self::Periapsis),
            "periapsisradius" => Ok(Self::PeriapsisRadius),
            "period" => Ok(Self::Period),
            "raan" => Ok(Self::RAAN),
            "rmag" => Ok(Self::Rmag),
            "semiparameter" => Ok(Self::SemiParameter),
            "sma" => Ok(Self::SMA),
            "trueanomaly" => Ok(Self::TrueAnomaly),
            "truelongitude" => Ok(Self::TrueLongitude),
            "vmag" => Ok(Self::Vmag),
            "x" => Ok(Self::X),
            "y" => Ok(Self::Y),
            "z" => Ok(Self::Z),
            "vx" => Ok(Self::VX),
            "vy" => Ok(Self::VY),
            "vz" => Ok(Self::VZ),
            _ => Err(NyxError::LoadingError(format!(
                "Unknown event state parameter: {}",
                s
            ))),
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
            StateParameter::AoL => angled_value(state.aol(), self.value),
            StateParameter::AoP => angled_value(state.aop(), self.value),
            StateParameter::Apoapsis => angled_value(state.ta(), 180.0),
            StateParameter::ApoapsisRadius => state.apoapsis() - self.value,
            StateParameter::EccentricAnomaly => angled_value(state.ea(), self.value),
            StateParameter::Eccentricity => state.ecc() - self.value,
            StateParameter::Energy => state.energy() - self.value,
            StateParameter::GeodeticHeight => state.geodetic_height() - self.value,
            StateParameter::GeodeticLatitude => state.geodetic_latitude() - self.value,
            StateParameter::GeodeticLongitude => state.geodetic_longitude() - self.value,
            StateParameter::Hmag => state.hmag() - self.value,
            StateParameter::HX => state.hx() - self.value,
            StateParameter::HY => state.hy() - self.value,
            StateParameter::HZ => state.hz() - self.value,
            StateParameter::Inclination => angled_value(state.inc(), self.value),
            StateParameter::MeanAnomaly => angled_value(state.ma(), self.value),
            StateParameter::Periapsis => angled_value(state.ta(), 0.0),
            StateParameter::PeriapsisRadius => state.periapsis() - self.value,
            StateParameter::Period => state.period().in_seconds() - self.value,
            StateParameter::RAAN => angled_value(state.raan(), self.value),
            StateParameter::Rmag => state.rmag() - self.value,
            StateParameter::SemiParameter => state.semi_parameter() - self.value,
            StateParameter::SMA => state.sma() - self.value,
            StateParameter::TrueAnomaly => angled_value(state.ta(), self.value),
            StateParameter::TrueLongitude => angled_value(state.tlong(), self.value),
            StateParameter::Vmag => state.vmag() - self.value,
            StateParameter::X => state.x - self.value,
            StateParameter::Y => state.y - self.value,
            StateParameter::Z => state.z - self.value,
            StateParameter::VX => state.vx - self.value,
            StateParameter::VY => state.vy - self.value,
            StateParameter::VZ => state.vz - self.value,
            _ => unimplemented!(),
        }
    }
}

impl EventEvaluator<SpacecraftState> for Event {
    fn eval(&self, state: &SpacecraftState) -> f64 {
        match self.parameter {
            StateParameter::FuelMass => state.fuel_mass_kg - self.value,
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
