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

pub use super::{Bodies, Cosm, Frame, LightTimeCalc, Orbit, Spacecraft};
use crate::md::EventEvaluator;
use crate::time::{Duration, Unit};
use std::cmp::{Eq, Ord, Ordering, PartialOrd};
use std::convert::Into;
use std::fmt;
use std::sync::Arc;

/// Stores the eclipse state
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EclipseState {
    Umbra,
    /// The f64 is between ]0; 1[ and corresponds to the percentage of penumbra: the closer to 1, the more light is seen.
    Penumbra(f64),
    Visibilis,
}

#[allow(clippy::from_over_into)]
impl Into<f64> for EclipseState {
    fn into(self) -> f64 {
        match self {
            Self::Umbra => 0.0,
            Self::Visibilis => 1.0,
            Self::Penumbra(val) => val,
        }
    }
}

impl Eq for EclipseState {}

impl Ord for EclipseState {
    /// Orders eclipse states to the greatest eclipse.
    ///
    /// *Examples*
    ///
    /// ```
    /// extern crate nyx_space as nyx;
    /// use nyx::cosmic::eclipse::EclipseState;
    /// assert!(EclipseState::Umbra == EclipseState::Umbra);
    /// assert!(EclipseState::Visibilis == EclipseState::Visibilis);
    /// assert!(EclipseState::Penumbra(0.5) == EclipseState::Penumbra(0.5));
    /// assert!(EclipseState::Umbra > EclipseState::Penumbra(0.1));
    /// assert!(EclipseState::Umbra > EclipseState::Penumbra(0.9));
    /// assert!(EclipseState::Penumbra(0.1) < EclipseState::Umbra);
    /// assert!(EclipseState::Penumbra(0.9) < EclipseState::Umbra);
    /// assert!(EclipseState::Penumbra(0.9) > EclipseState::Visibilis);
    /// assert!(EclipseState::Visibilis < EclipseState::Penumbra(0.9));
    /// ```
    fn cmp(&self, other: &Self) -> Ordering {
        match *self {
            EclipseState::Umbra => {
                if *other == EclipseState::Umbra {
                    Ordering::Equal
                } else {
                    Ordering::Greater
                }
            }
            EclipseState::Visibilis => {
                if *other == EclipseState::Visibilis {
                    Ordering::Equal
                } else {
                    Ordering::Less
                }
            }
            EclipseState::Penumbra(s) => match *other {
                EclipseState::Penumbra(o) => {
                    if s > o {
                        Ordering::Greater
                    } else {
                        Ordering::Less
                    }
                }
                EclipseState::Visibilis => Ordering::Greater,
                EclipseState::Umbra => Ordering::Less,
            },
        }
    }
}

impl PartialOrd for EclipseState {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

impl fmt::Display for EclipseState {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Self::Umbra => write!(f, "Umbra"),
            Self::Visibilis => write!(f, "Visibilis"),
            Self::Penumbra(v) => write!(f, "Penumbra {:.2}%", v * 100.0),
        }
    }
}

#[derive(Clone)]
pub struct EclipseLocator {
    pub light_source: Frame,
    pub shadow_bodies: Vec<Frame>,
    pub cosm: Arc<Cosm>,
}

impl fmt::Display for EclipseLocator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let shadow_bodies = self
            .shadow_bodies
            .iter()
            .map(|b| format!("{}, ", b))
            .collect::<String>();
        write!(
            f,
            "light-source: {}, shadows casted by: {}",
            self.light_source, shadow_bodies
        )
    }
}

impl EclipseLocator {
    /// Compute the visibility/eclipse between an observer and an observed state
    pub fn compute(&self, observer: &Orbit) -> EclipseState {
        let mut state = EclipseState::Visibilis;
        for eclipsing_body in &self.shadow_bodies {
            let this_state =
                eclipse_state(observer, self.light_source, *eclipsing_body, &self.cosm);
            if this_state > state {
                state = this_state;
            }
        }
        state
    }

    /// Creates an umbra event from this eclipse locator
    pub fn to_umbra_event(&self) -> UmbraEvent {
        UmbraEvent {
            e_loc: self.clone(),
        }
    }

    /// Creates a penumbra event from this eclipse locator
    pub fn to_penumbra_event(&self) -> PenumbraEvent {
        PenumbraEvent {
            e_loc: self.clone(),
        }
    }
}

/// An event to find the darkest eclipse state (more than 98% shadow)
pub struct UmbraEvent {
    e_loc: EclipseLocator,
}

impl fmt::Display for UmbraEvent {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "umbra event {}", self.e_loc)
    }
}

impl EventEvaluator<Orbit> for UmbraEvent {
    // Evaluation of the event, returns 0.0 for umbra, 1.0 for visibility and some value in between for penumbra
    fn eval(&self, observer: &Orbit) -> f64 {
        match self.e_loc.compute(observer) {
            EclipseState::Umbra => 0.0,
            EclipseState::Visibilis => 1.0,
            EclipseState::Penumbra(val) => val,
        }
    }

    /// Stop searching when the time has converged to less than 0.1 seconds
    fn epoch_precision(&self) -> Duration {
        0.1 * Unit::Second
    }
    /// Finds the darkest part of an eclipse within 2% of penumbra (i.e. 98% in shadow)
    fn value_precision(&self) -> f64 {
        0.02
    }
}

impl EventEvaluator<Spacecraft> for UmbraEvent {
    // Evaluation of the event, returns 0.0 for umbra, 1.0 for visibility and some value in between for penumbra
    fn eval(&self, sc: &Spacecraft) -> f64 {
        match self.e_loc.compute(&sc.orbit) {
            EclipseState::Umbra => 0.0,
            EclipseState::Visibilis => 1.0,
            EclipseState::Penumbra(val) => val,
        }
    }

    /// Stop searching when the time has converged to less than 0.1 seconds
    fn epoch_precision(&self) -> Duration {
        0.1 * Unit::Second
    }
    /// Finds the darkest part of an eclipse within 2% of penumbra (i.e. 98% in shadow)
    fn value_precision(&self) -> f64 {
        0.02
    }
}

/// An event to find the start of a penumbra
pub struct PenumbraEvent {
    e_loc: EclipseLocator,
}

impl fmt::Display for PenumbraEvent {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "penumbra event {}", self.e_loc)
    }
}

impl EventEvaluator<Orbit> for PenumbraEvent {
    fn eval(&self, observer: &Orbit) -> f64 {
        match self.e_loc.compute(observer) {
            EclipseState::Umbra => 0.0,
            EclipseState::Visibilis => 1.0,
            EclipseState::Penumbra(val) => val - 1.0,
        }
    }

    /// Stop searching when the time has converged to less than 0.1 seconds
    fn epoch_precision(&self) -> Duration {
        0.1 * Unit::Second
    }
    /// Finds the slightest penumbra within 2%(i.e. 98% in visibility)
    fn value_precision(&self) -> f64 {
        0.02
    }
}

impl EventEvaluator<Spacecraft> for PenumbraEvent {
    fn eval(&self, sc: &Spacecraft) -> f64 {
        match self.e_loc.compute(&sc.orbit) {
            EclipseState::Umbra => 0.0,
            EclipseState::Visibilis => 1.0,
            EclipseState::Penumbra(val) => val - 1.0,
        }
    }

    /// Stop searching when the time has converged to less than 0.1 seconds
    fn epoch_precision(&self) -> Duration {
        0.1 * Unit::Second
    }
    /// Finds the slightest penumbra within 2%(i.e. 98% in visibility)
    fn value_precision(&self) -> f64 {
        0.02
    }
}

/// Computes the umbra/visibilis/penumbra state between between two states accounting for eclipsing of the providing geoid.
pub fn eclipse_state(
    observer: &Orbit,
    light_source: Frame,
    eclipsing_body: Frame,
    cosm: &Cosm,
) -> EclipseState {
    // If the light source's radius is zero, just call the line of sight algorithm

    assert!(light_source.is_geoid() || light_source.is_celestial());
    assert!(eclipsing_body.is_geoid());

    if light_source.equatorial_radius() < std::f64::EPSILON {
        let observed = cosm.celestial_state(
            &light_source.ephem_path(),
            observer.dt,
            observer.frame,
            LightTimeCalc::None,
        );
        return line_of_sight(observer, &observed, eclipsing_body, cosm);
    }
    // All of the computations happen with the observer as the center.
    // `eb` stands for eclipsing body; `ls` stands for light source.
    // Get the radius vector of the spacecraft to the eclipsing body
    let r_eb = cosm.frame_chg(observer, eclipsing_body).radius();

    // Get the radius vector of the light source to the spacecraft
    let r_ls = -cosm.frame_chg(observer, light_source).radius();

    // Compute the apparent radii of the light source and eclipsing body (preventing any NaN)
    let r_ls_prime = if light_source.equatorial_radius() >= r_ls.norm() {
        light_source.equatorial_radius()
    } else {
        (light_source.equatorial_radius() / r_ls.norm()).asin()
    };
    let r_eb_prime = if eclipsing_body.equatorial_radius() >= r_eb.norm() {
        eclipsing_body.equatorial_radius()
    } else {
        (eclipsing_body.equatorial_radius() / r_eb.norm()).asin()
    };

    // Compute the apparent separation of both circles
    let d_prime = (-(r_ls.dot(&r_eb)) / (r_eb.norm() * r_ls.norm())).acos();

    if d_prime - r_ls_prime > r_eb_prime {
        // If the closest point where the apparent radius of the light source _starts_ is further
        // away than the furthest point where the eclipsing body's shadow can reach, then the light
        // source is totally visible.
        EclipseState::Visibilis
    } else if r_eb_prime > d_prime + r_ls_prime {
        // The light source is fully hidden by the eclipsing body, hence we're in total eclipse.
        EclipseState::Umbra
    } else if (r_ls_prime - r_eb_prime).abs() < d_prime && d_prime < r_ls_prime + r_eb_prime {
        // If we have reached this point, we're in penumbra.
        // Both circles, which represent the light source projected onto the plane and the eclipsing geoid,
        // now overlap creating an asymmetrial lens.
        // The following math comes from http://mathworld.wolfram.com/Circle-CircleIntersection.html
        // and https://stackoverflow.com/questions/3349125/circle-circle-intersection-points .

        // Compute the distances between the center of the eclipsing geoid and the line crossing the intersection
        // points of both circles.
        let d1 = (d_prime.powi(2) - r_ls_prime.powi(2) + r_eb_prime.powi(2)) / (2.0 * d_prime);
        let d2 = (d_prime.powi(2) + r_ls_prime.powi(2) - r_eb_prime.powi(2)) / (2.0 * d_prime);

        let shadow_area = circ_seg_area(r_eb_prime, d1) + circ_seg_area(r_ls_prime, d2);
        if shadow_area.is_nan() {
            warn!(
                "Shadow area is NaN! Please file a bug with initial states, eclipsing bodies, etc."
            );
            return EclipseState::Umbra;
        }
        // Compute the nominal area of the light source
        let nominal_area = std::f64::consts::PI * r_ls_prime.powi(2);
        // And return the percentage (between 0 and 1) of the eclipse.
        EclipseState::Penumbra(1.0 - shadow_area / nominal_area)
    } else {
        // Annular eclipse.
        // If r_eb_prime is very small, then the fraction is very small: however, we note a penumbra close to 1.0 as near full light source visibility, so let's subtract one from this.
        EclipseState::Penumbra(1.0 - r_eb_prime.powi(2) / r_ls_prime.powi(2))
    }
}

// Compute the area of the circular segment of radius r and chord length d
fn circ_seg_area(r: f64, d: f64) -> f64 {
    r.powi(2) * (d / r).acos() - d * (r.powi(2) - d.powi(2)).sqrt()
}

/// Computes the light of sight the provided time between two states accounting for eclipsing of the providing geoid.
/// This works for visibility between spacecraft and a ground station. For eclipsing and penumbras, use `eclipse_state`.
///
/// Source: Algorithm 35 of Vallado, 4th edition, page 308.
pub fn line_of_sight(
    observer: &Orbit,
    observed: &Orbit,
    eclipsing_body: Frame,
    cosm: &Cosm,
) -> EclipseState {
    if observer == observed {
        return EclipseState::Visibilis;
    }

    // Convert the states to the same frame as the eclipsing body (ensures we're in the same frame)
    let r1 = &cosm.frame_chg(observed, eclipsing_body).radius();
    let r2 = &cosm.frame_chg(observer, eclipsing_body).radius();

    let r1sq = r1.dot(r1);
    let r2sq = r2.dot(r2);
    let r1dotr2 = r1.dot(r2);

    let tau = (r1sq - r1dotr2) / (r1sq + r2sq - 2.0 * r1dotr2);
    if !(0.0..=1.0).contains(&tau)
        || (1.0 - tau) * r1sq + r1dotr2 * tau > eclipsing_body.equatorial_radius().powi(2)
    {
        EclipseState::Visibilis
    } else {
        EclipseState::Umbra
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hifitime::Epoch;

    #[test]
    fn los_edge_case() {
        let cosm = Cosm::de438();
        let eme2k = cosm.frame("EME2000");
        let luna = cosm.frame("Luna");

        let dt1 = Epoch::from_gregorian_tai(2020, 1, 1, 6, 7, 40, 0);
        let dt2 = Epoch::from_gregorian_tai(2020, 1, 1, 6, 7, 50, 0);
        let dt3 = Epoch::from_gregorian_tai(2020, 1, 1, 6, 8, 0, 0);

        let xmtr1 = Orbit::cartesian(
            397_477.494_485,
            -57_258.902_156,
            -62_857.909_437,
            0.230_482,
            2.331_362,
            0.615_501,
            dt1,
            eme2k,
        );
        let rcvr1 = Orbit::cartesian(
            338_335.467_589,
            -55_439.526_977,
            -13_327.354_273,
            0.197_141,
            0.944_261,
            0.337_407,
            dt1,
            eme2k,
        );
        let xmtr2 = Orbit::cartesian(
            397_479.756_900,
            -57_235.586_465,
            -62_851.758_851,
            0.222_000,
            2.331_768,
            0.614_614,
            dt2,
            eme2k,
        );
        let rcvr2 = Orbit::cartesian(
            338_337.438_860,
            -55_430.084_340,
            -13_323.980_229,
            0.197_113,
            0.944_266,
            0.337_402,
            dt2,
            eme2k,
        );
        let xmtr3 = Orbit::cartesian(
            397_481.934_480,
            -57_212.266_970,
            -62_845.617_185,
            0.213_516,
            2.332_122,
            0.613_717,
            dt3,
            eme2k,
        );
        let rcvr3 = Orbit::cartesian(
            338_339.409_858,
            -55_420.641_651,
            -13_320.606_228,
            0.197_086,
            0.944_272,
            0.337_398,
            dt3,
            eme2k,
        );

        assert_eq!(
            line_of_sight(&xmtr1, &rcvr1, luna, &cosm),
            EclipseState::Umbra
        );
        assert_eq!(
            line_of_sight(&xmtr2, &rcvr2, luna, &cosm),
            EclipseState::Umbra
        );
        assert_eq!(
            line_of_sight(&xmtr3, &rcvr3, luna, &cosm),
            EclipseState::Umbra
        );

        // Test converse

        assert_eq!(
            line_of_sight(&rcvr1, &xmtr1, luna, &cosm),
            EclipseState::Umbra
        );
        assert_eq!(
            line_of_sight(&rcvr2, &xmtr2, luna, &cosm),
            EclipseState::Umbra
        );
        assert_eq!(
            line_of_sight(&rcvr3, &xmtr3, luna, &cosm),
            EclipseState::Umbra
        );
    }

    #[test]
    fn los_earth_eclipse() {
        let cosm = Cosm::de438();
        let eme2k = cosm.frame("EME2000");

        let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

        let sma = eme2k.equatorial_radius() + 300.0;

        let sc1 = Orbit::keplerian(sma, 0.001, 0.1, 90.0, 75.0, 0.0, dt, eme2k);
        let sc2 = Orbit::keplerian(sma + 1.0, 0.001, 0.1, 90.0, 75.0, 0.0, dt, eme2k);
        let sc3 = Orbit::keplerian(sma, 0.001, 0.1, 90.0, 75.0, 180.0, dt, eme2k);

        // Out of phase by pi.
        assert_eq!(line_of_sight(&sc1, &sc3, eme2k, &cosm), EclipseState::Umbra);

        assert_eq!(
            line_of_sight(&sc2, &sc1, eme2k, &cosm),
            EclipseState::Visibilis
        );

        // Nearly identical orbits in the same phasing
        assert_eq!(
            line_of_sight(&sc1, &sc2, eme2k, &cosm),
            EclipseState::Visibilis
        );
    }
}
