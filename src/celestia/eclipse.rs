use super::{Cosm, Geoid, State};
use std::cmp::{Eq, Ord, Ordering, PartialOrd};

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

pub struct EclipseLocator<'a> {
    pub light_source: Geoid,
    pub shadow_bodies: Vec<Geoid>,
    pub cosm: &'a Cosm,
}

impl<'a> EclipseLocator<'a> {
    /// Compute the visibility/eclipse between an observer and an observed state
    pub fn compute(&self, observer: &State<Geoid>) -> EclipseState {
        let mut state = EclipseState::Visibilis;
        for eclipsing_geoid in &self.shadow_bodies {
            let this_state =
                eclipse_state(observer, self.light_source, *eclipsing_geoid, self.cosm);
            if this_state > state {
                state = this_state;
            }
        }
        state
    }
}

/// Computes the umbra/visibilis/penumbra state between between two states accounting for eclipsing of the providing geoid.
///
/// # Algorithm
/// TODO: add derivation
pub fn eclipse_state(
    observer: &State<Geoid>,
    light_source: Geoid,
    eclipsing_geoid: Geoid,
    cosm: &Cosm,
) -> EclipseState {
    // If the light source's radius is zero, just call the line of sight algorithm
    if light_source.equatorial_radius < std::f64::EPSILON {
        let observed = cosm.celestial_state(
            light_source.id,
            observer.dt.as_jde_et_days(),
            observer.frame.id,
        );
        return line_of_sight(observer, &observed, eclipsing_geoid, &cosm);
    }
    // All of the computations happen with the observer as the center.
    // `eb` stands for eclipsing body; `ls` stands for light source.
    // Vector from spacecraft to EB
    let r_eb = -cosm.frame_chg(observer, eclipsing_geoid).radius();
    let r_eb_unit = r_eb / r_eb.norm();
    // Vector from EB to LS
    let r_eb_ls = cosm
        .celestial_state(
            light_source.id,
            observer.dt.as_jde_et_days(),
            eclipsing_geoid.id,
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
    let pseudo_ls_radius = r_ls_p.norm() * light_source.equatorial_radius / r_ls.norm();
    // And then project that onto the plane to compute the actual project radius of the light source onto the plane
    let ls_radius = pseudo_ls_radius / gamma.cos();
    // Compute the radius from the center of the eclipsing geoid to the intersection point
    let r_plane_ls = r_ls_p - r_eb;
    // Note that the eclipsing geoid's circle is centered at zero, so the norm of r_plane_ls is also the distance
    // between the center of the eclipsing body's shadow and the center of the light source's shadow.
    let d_plane_ls = r_plane_ls.norm();

    let eb_radius = eclipsing_geoid.equatorial_radius;

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
    fn los_trivial() {
        let cosm = Cosm::from_xb("./de438s");
        let mut earth = cosm.geoid_from_id(399);
        earth.equatorial_radius = 1.0;

        let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

        let x1 = State::<Geoid>::from_cartesian(-1.0, -1.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let x2 = State::<Geoid>::from_cartesian(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let x3 = State::<Geoid>::from_cartesian(-2.0, 1.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let x4 = State::<Geoid>::from_cartesian(-3.0, 2.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let x5 = State::<Geoid>::from_cartesian(1.0, -2.0, 0.0, 0.0, 0.0, 0.0, dt, earth);

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
    }

    #[test]
    fn los_earth_eclipse() {
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

    #[test]
    fn eclipse_sun_eclipse() {
        let cosm = Cosm::from_xb("./de438s");
        let earth = cosm.geoid_from_id(399);
        let sun = cosm.geoid_from_id(10);

        let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

        let sma = earth.equatorial_radius + 300.0;

        let sc1 = State::<Geoid>::from_keplerian(sma, 0.001, 0.1, 90.0, 75.0, 25.0, dt, earth);
        let sc2 = State::<Geoid>::from_keplerian(sma, 0.001, 0.1, 90.0, 75.0, 115.0, dt, earth);
        let sc3 = State::<Geoid>::from_keplerian(sma, 0.001, 0.1, 90.0, 75.0, 77.2, dt, earth);

        assert_eq!(
            eclipse_state(&sc1, sun, earth, &cosm),
            EclipseState::Visibilis
        );
        assert_eq!(eclipse_state(&sc2, sun, earth, &cosm), EclipseState::Umbra);
        match eclipse_state(&sc3, sun, earth, &cosm) {
            EclipseState::Penumbra(val) => assert!(val > 0.9),
            _ => panic!("should be in penumbra"),
        };
    }
}
