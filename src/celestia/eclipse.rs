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
    pub tolerance: f64,
    pub cosm: &'a Cosm,
}

impl<'a> EclipseLocator<'a> {
    /// Compute the visibility/eclipse between an observer and an observed state
    pub fn compute(&self, observer: &State<Geoid>, observed: &State<Geoid>) -> EclipseState {
        let mut state = EclipseState::Visibilis;
        for eclipsing_geoid in &self.shadow_bodies {
            let this_state = eclipse_state(
                observer,
                observed,
                *eclipsing_geoid,
                self.tolerance,
                self.cosm,
            );
            if this_state > state {
                state = this_state;
            }
        }
        state
    }
}

/// Computes the state of an eclipse at the provided time between two states, at a given tolerance.
/// For perfect geometric eclipsing, the tolerance should be set to 0.0. In practice, value should be less than 0.5.
///
/// **Example**: to compute whether a spacecraft can see the Moon despite the Earth maybe in the way,
/// `observer` should be the spacecraft state, `observed` should be the Moon state (at the same time),
/// and `eclipsing_geoid` should be the Earth.
/// **Assumption:** the distance and size of both the observer and observed states is such that they respectively look like single points.
pub fn eclipse_state(
    observer: &State<Geoid>,
    observed: &State<Geoid>,
    eclipsing_geoid: Geoid,
    tolerance: f64,
    cosm: &Cosm,
) -> EclipseState {
    if observer == observed {
        return EclipseState::Visibilis;
    }

    // Convert the observed to the same frame as the origin (fok = frame OK)
    let observed_fok = &cosm.frame_chg(observed, observer.frame);
    // Compute the unit direction between both
    let mut l = observed_fok.radius() - observer.radius();
    l /= l.norm();

    // Get the eclipsing body position in the same frame (allows us to set the origin to a null vector)
    let eclipsing_geoid_fok = cosm.celestial_state(
        eclipsing_geoid.id,
        observer.dt.as_jde_et_days(),
        observer.frame.id,
    );

    let omc = observer.radius() - eclipsing_geoid_fok.radius();
    let r = eclipsing_geoid.equatorial_radius;

    // l.dot(&l) should be 1.0, within rounding error.
    let discriminant_sq = (l.dot(&omc)).powi(2) - l.dot(&l) * (omc.dot(&omc) - (r.powi(2)));
    if discriminant_sq < 0.0 {
        // No intersection between the direction of the origin and vis_chk objects and the sphere
        EclipseState::Visibilis
    } else {
        // Compute the distance from the origin to the intersection point.
        let intersect_dist = -(l.dot(&omc)) + discriminant_sq.sqrt();
        // And the intersection point itself.
        let intersect_pt = observer.radius() + intersect_dist * l;
        let dist_ori_inters = observer.distance_to_point(&intersect_pt);
        let dist_oth_inters = observed_fok.distance_to_point(&intersect_pt);
        let dist_ori_oth = observer.distance_to(&observed_fok);
        // If the intersection point is on the line between both, then it causes an eclipse.
        if (dist_ori_inters + dist_oth_inters - dist_ori_oth).abs() < 1e-15 {
            if discriminant_sq > 0.0 {
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

        assert_eq!(
            eclipse_state(&x1, &x2, earth, 0.0, &cosm),
            EclipseState::Umbra
        );
        assert_eq!(
            eclipse_state(&x1, &x3, earth, 0.0, &cosm),
            EclipseState::Visibilis
        );
        assert_eq!(
            eclipse_state(&x3, &x2, earth, 0.0, &cosm),
            EclipseState::Penumbra(0.5)
        );
        assert_eq!(
            eclipse_state(&x4, &x3, earth, 0.0, &cosm),
            EclipseState::Visibilis
        );
        assert_eq!(
            eclipse_state(&x3, &x4, earth, 0.0, &cosm),
            EclipseState::Visibilis
        );
        assert_eq!(
            eclipse_state(&x5, &x4, earth, 0.0, &cosm),
            EclipseState::Umbra
        );
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
        assert_eq!(
            eclipse_state(&sc1, &sc3, earth, 0.0, &cosm),
            EclipseState::Umbra
        );

        // Nearly identical orbits in the same phasing
        assert_eq!(
            eclipse_state(&sc1, &sc2, earth, 0.0, &cosm),
            EclipseState::Visibilis
        );
    }
}
