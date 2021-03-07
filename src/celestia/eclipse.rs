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

pub use super::{Cosm, Frame, LTCorr, Orbit, Spacecraft};
use crate::md::EventEvaluator;
use crate::time::{Duration, TimeUnit};
use std::cmp::{Eq, Ord, Ordering, PartialOrd};
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

impl Eq for EclipseState {}

impl Ord for EclipseState {
    /// Orders eclipse states to the greatest eclipse.
    ///
    /// *Examples*
    ///
    /// ```
    /// extern crate nyx_space as nyx;
    /// use nyx::celestia::eclipse::EclipseState;
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
    pub correction: LTCorr,
}

impl fmt::Display for EclipseLocator {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        let shadow_bodies = self
            .shadow_bodies
            .iter()
            .map(|b| format!("{} ", b))
            .collect::<String>();
        write!(
            f,
            "light-source: {}\tshadows casted by: {} (light time corr.: {:?})",
            self.light_source, shadow_bodies, self.correction
        )
    }
}

impl EclipseLocator {
    /// Compute the visibility/eclipse between an observer and an observed state
    pub fn compute(&self, observer: &Orbit) -> EclipseState {
        let mut state = EclipseState::Visibilis;
        for eclipsing_body in &self.shadow_bodies {
            let this_state = eclipse_state(
                observer,
                self.light_source,
                *eclipsing_body,
                &self.cosm,
                self.correction,
            );
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
        0.1 * TimeUnit::Second
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
        0.1 * TimeUnit::Second
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
        0.1 * TimeUnit::Second
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
        0.1 * TimeUnit::Second
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
    correction: LTCorr,
) -> EclipseState {
    // If the light source's radius is zero, just call the line of sight algorithm

    assert!(light_source.is_geoid() || light_source.is_celestial());
    assert!(eclipsing_body.is_geoid());

    if light_source.equatorial_radius() < std::f64::EPSILON {
        let observed = cosm.celestial_state(
            &light_source.ephem_path(),
            observer.dt,
            observer.frame,
            correction,
        );
        return line_of_sight(observer, &observed, eclipsing_body, &cosm);
    }
    // All of the computations happen with the observer as the center.
    // `eb` stands for eclipsing body; `ls` stands for light source.
    // Vector from spacecraft to EB
    let r_eb = -cosm.frame_chg(observer, eclipsing_body).radius();
    let r_eb_unit = r_eb / r_eb.norm();
    // Vector from EB to LS
    let r_eb_ls = cosm
        .celestial_state(
            &light_source.ephem_path(),
            observer.dt,
            eclipsing_body,
            correction,
        )
        .radius();
    let r_eb_ls_unit = r_eb_ls / r_eb_ls.norm();
    // Compute the angle between those vectors. If that angle is less than 90 degrees, then the light source
    // is in front of the eclipsing body, so we're visible.
    let beta2 = r_eb_unit.dot(&r_eb_ls_unit).acos();
    if beta2 <= std::f64::consts::FRAC_PI_2 {
        return EclipseState::Visibilis;
    }
    // Get the state of the light source with respect to the spacecraft
    let r_ls = r_eb + r_eb_ls;
    let r_ls_unit = r_ls / r_ls.norm();
    // Compute beta3, the angle between the spacecraft and the eclipsing body, and between the spacecraft and the light source.
    // We need this to project the radius of the light source onto the plane centered at the eclipsing geoid, and normal to
    // the direction to the spacecraft.
    let cos_beta3 = r_ls_unit.dot(&r_eb_unit);
    // Now compute r_ls_p, the vector from the eclipsing body to the center of the projected light source onto the plane.
    let r_ls_p = (r_eb.norm() / cos_beta3) * r_ls_unit;
    // Now let's compute the pseudo radius of the light source near the plane.
    // This is the pseudo radius because we're building a right triangle between the intersection point of r_ls with the plane,
    // the normal from r_ls at that point until intersection with r_eb_ls.
    let beta1 = (-r_ls_unit).dot(&-r_eb_ls_unit).acos();
    // We can now compute the gamma angle
    let gamma = std::f64::consts::FRAC_PI_2 - beta2 - beta1;
    // Applying Thales theorem, we can compute the pseudo radius of the light source
    let pseudo_ls_radius = r_ls_p.norm() * light_source.equatorial_radius() / r_ls.norm();
    // And then project that onto the plane to compute the actual project radius of the light source onto the plane
    let ls_radius = pseudo_ls_radius / gamma.cos();
    // Compute the radius from the center of the eclipsing geoid to the intersection point
    let r_plane_ls = r_ls_p - r_eb;
    // Note that the eclipsing geoid's circle is centered at zero, so the norm of r_plane_ls is also the distance
    // between the center of the eclipsing body's shadow and the center of the light source's shadow.
    let d_plane_ls = r_plane_ls.norm();

    let eb_radius = eclipsing_body.equatorial_radius();

    if d_plane_ls - ls_radius > eb_radius {
        // If the closest point where the projected light source's circle _starts_ is further
        // away than the furthest point where the eclipsing body's shadow can reach, then the light
        // source is totally visible.
        return EclipseState::Visibilis;
    } else if eb_radius > d_plane_ls + ls_radius {
        // The light source is fully hidden by the eclipsing body, hence we're in total eclipse.
        // Note that because we test for the beta2 angle earlier, we know that the light source is behind the plane
        // from the vantage point of the spacecraft.
        return EclipseState::Umbra;
    }

    // If we have reached this point, we're in penumbra.
    // Both circles, which represent the light source projected onto the plane and the eclipsing geoid,
    // now overlap creating an asymmetrial lens.
    // The following math comes from http://mathworld.wolfram.com/Circle-CircleIntersection.html
    // and https://stackoverflow.com/questions/3349125/circle-circle-intersection-points .

    // Compute the distances between the center of the eclipsing geoid and the line crossing the intersection
    // points of both circles.
    let d1 = (d_plane_ls.powi(2) - ls_radius.powi(2) + eb_radius.powi(2)) / (2.0 * d_plane_ls);
    let d2 = (d_plane_ls.powi(2) + ls_radius.powi(2) - eb_radius.powi(2)) / (2.0 * d_plane_ls);

    let shadow_area = circ_seg_area(eb_radius, d1) + circ_seg_area(ls_radius, d2);
    if shadow_area.is_nan() {
        return EclipseState::Umbra;
    }
    // Compute the nominal area of the light source
    let nominal_area = std::f64::consts::PI * ls_radius.powi(2);
    // And return the percentage (between 0 and 1) of the eclipse.
    EclipseState::Penumbra((nominal_area - shadow_area) / nominal_area)
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

    let r1sq = r1.dot(&r1);
    let r2sq = r2.dot(&r2);
    let r1dotr2 = r1.dot(&r2);

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

    #[test]
    fn eclipse_sun_eclipse() {
        let cosm = Cosm::de438();
        let sun = cosm.frame("Sun J2000");
        let eme2k = cosm.frame("EME2000");

        let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

        let sma = eme2k.equatorial_radius() + 300.0;

        let sc1 = Orbit::keplerian(sma, 0.001, 0.1, 90.0, 75.0, 25.0, dt, eme2k);
        let sc2 = Orbit::keplerian(sma, 0.001, 0.1, 90.0, 75.0, 115.0, dt, eme2k);
        let sc3 = Orbit::keplerian(sma, 0.001, 0.1, 90.0, 75.0, 77.2, dt, eme2k);

        let correction = LTCorr::None;

        assert_eq!(
            eclipse_state(&sc1, sun, eme2k, &cosm, correction),
            EclipseState::Visibilis
        );
        assert_eq!(
            eclipse_state(&sc2, sun, eme2k, &cosm, correction),
            EclipseState::Umbra
        );
        match eclipse_state(&sc3, sun, eme2k, &cosm, correction) {
            EclipseState::Penumbra(val) => assert!(val > 0.9),
            _ => panic!("should be in penumbra"),
        };
    }
}
