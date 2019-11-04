use super::{Cosm, Geoid, State};
use std::cmp::{Eq, Ord, Ordering, PartialOrd};

/// Stores the eclipse state
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum EclipseState {
    Umbra,
    /// The f64 is between ]0; 1[ and corresponds to the percentage of penumbra based on the tolerance.
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

pub struct EclipseLocator<'a> {
    pub shadow_bodies: Vec<Geoid>,
    pub cosm: &'a Cosm,
}

impl<'a> EclipseLocator<'a> {
    /// Compute the visibility/eclipse between an observer and an observed state
    pub fn compute(&self, observer: &State<Geoid>, observed: &State<Geoid>) -> EclipseState {
        let mut state = EclipseState::Visibilis;
        for eclipsing_geoid in &self.shadow_bodies {
            let this_state = line_of_sight(observer, observed, *eclipsing_geoid, self.cosm);
            if this_state > state {
                state = this_state;
            }
        }
        state
    }
}

/// Computes the light of sight the provided time between two states accounting for eclipsing of the providing geoid.
/// This works for visibility between spacecraft and a ground station. For eclipsing, use `eclipse_state`.
///
/// # Algorithm
/// This function solves the [Line-Sphere intersection](https://en.wikipedia.org/wiki/Line%E2%80%93sphere_intersection) problem.
/// 1. We create a line between the observer and observed.
/// 2. We then check if the line will intersect the geoid, represented as a sphere.
/// 3. We compute the discriminant from the quadratic formula which emanates from solving the distance at which the intersection of the line and the sphere will happen.
/// 4. If the discriminant is less than 0.0, then the distance along the starting point of the line to the sphere is an complex number, i.e. no intersection.
/// 5. If the discriminant is greater than 0.0, then there are exactly two intersection points between the line and the sphere.
/// 6. If the discriminant is exactly zero, then the line "skims" the sphere.
///
/// In practice, this code uses a tolerance of 1e-12 instead of 0.0 to account for rounding issues.
pub fn line_of_sight(
    observer: &State<Geoid>,
    observed: &State<Geoid>,
    eclipsing_geoid: Geoid,
    cosm: &Cosm,
) -> EclipseState {
    if observer == observed {
        return EclipseState::Visibilis;
    }

    // Convert the observed to the same frame as the origin (fok = frame OK)
    let observed_fok = &cosm.frame_chg(observed, eclipsing_geoid);
    let observer_fok = &cosm.frame_chg(observer, eclipsing_geoid);
    // Compute the unit direction between both
    let mut l = observed_fok.radius() - observer_fok.radius();
    l /= l.norm();

    // observer minus center point of the eclipsing body is already in the eclipsing body frame, so center is zero.
    let omc = observer_fok.radius();
    let r = eclipsing_geoid.equatorial_radius;

    // l.dot(&l) should be 1.0, within rounding error.
    let discriminant_sq = (l.dot(&omc)).powi(2) - l.dot(&l) * (omc.dot(&omc) - (r.powi(2)));

    if discriminant_sq < -1e-12 {
        // No intersection between the direction of the origin and vis_chk objects and the sphere
        EclipseState::Visibilis
    } else {
        // Compute the distance from the origin to the intersection point.
        let intersect_dist = -(l.dot(&omc)) + discriminant_sq.sqrt();
        // And the intersection point itself.
        let intersect_pt = observer_fok.radius() + intersect_dist * l;
        let dist_ori_inters = observer_fok.distance_to_point(&intersect_pt);
        let dist_oth_inters = observed_fok.distance_to_point(&intersect_pt);
        let dist_ori_oth = observer_fok.distance_to(&observed_fok);
        // If the intersection point is on the line between both, then it causes an eclipse.
        if (dist_ori_inters + dist_oth_inters - dist_ori_oth).abs() < 1e-15 {
            if discriminant_sq > 1e-12 {
                EclipseState::Umbra
            } else {
                EclipseState::Penumbra(0.5)
            }
        } else {
            EclipseState::Visibilis
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hifitime::Epoch;

    #[test]
    fn trivial() {
        let cosm = Cosm::from_xb("./de438s");
        let mut earth = cosm.geoid_from_id(399);
        earth.equatorial_radius = 1.0;

        let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

        let x1 = State::<Geoid>::from_cartesian(-1.0, -1.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let x2 = State::<Geoid>::from_cartesian(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let x3 = State::<Geoid>::from_cartesian(-2.0, 1.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let x4 = State::<Geoid>::from_cartesian(-3.0, 2.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let x5 = State::<Geoid>::from_cartesian(1.0, -2.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let x6 = State::<Geoid>::from_cartesian(-1.0, 0.99999, 0.0, 0.0, 0.0, 0.0, dt, earth);

        assert_eq!(line_of_sight(&x1, &x2, earth, &cosm), EclipseState::Umbra);
        assert_eq!(
            line_of_sight(&x1, &x3, earth, &cosm),
            EclipseState::Visibilis
        );
        assert_eq!(
            line_of_sight(&x3, &x2, earth, &cosm),
            EclipseState::Penumbra(0.5)
        );
        assert_eq!(
            line_of_sight(&x4, &x3, earth, &cosm),
            EclipseState::Visibilis
        );
        assert_eq!(
            line_of_sight(&x3, &x4, earth, &cosm),
            EclipseState::Visibilis
        );
        assert_eq!(line_of_sight(&x5, &x4, earth, &cosm), EclipseState::Umbra);
        match line_of_sight(&x6, &x2, earth, &cosm) {
            EclipseState::Penumbra(perc) => assert!((0.5 - perc).abs() < 1e-2),
            _ => panic!("bug"),
        };
    }

    #[test]
    fn earth_eclipse() {
        let cosm = Cosm::from_xb("./de438s");
        let earth = cosm.geoid_from_id(399);

        let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

        let sma = earth.equatorial_radius + 300.0;

        let sc1 = State::<Geoid>::from_keplerian(sma, 0.001, 0.1, 90.0, 75.0, 0.0, dt, earth);
        let sc2 = State::<Geoid>::from_keplerian(sma + 1.0, 0.001, 0.1, 90.0, 75.0, 0.0, dt, earth);
        let sc3 = State::<Geoid>::from_keplerian(sma, 0.001, 0.1, 90.0, 75.0, 180.0, dt, earth);

        // Out of phase by pi.
        assert_eq!(line_of_sight(&sc1, &sc3, earth, &cosm), EclipseState::Umbra);

        // Nearly identical orbits in the same phasing
        assert_eq!(
            line_of_sight(&sc1, &sc2, earth, &cosm),
            EclipseState::Visibilis
        );
    }
}
