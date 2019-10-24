use super::{Cosm, Geoid, State};

/// Stores the eclipse state
#[derive(Clone, Copy, Debug)]
pub enum EclipseState {
    Umbra,
    Penumbra,
    Visibilis,
}

/// Computes the state of an eclipse at the provided time between two states.
///
/// Example: to compute whether a spacecraft can see the Moon despite the Earth maybe in the way,
/// origin_obj should be the spacecraft state, vis_chk_obj should be the Moon state (at the same time),
/// and eclipsing_geoid should be the Earth.
pub fn eclipse_state(origin_obj: State<Geoid>, vis_chk_obj: State<Geoid>, eclipsing_geoid: Geoid, cosm: &Cosm) -> EclipseState {
    // Convert the vis_chk_obj to the same frame as the origin.
    let vis_chk_obj_fok = &cosm.frame_chg(vis_chk_obj, origin_obj.frame);
    // Compute the unit direction between both
    let mut l = dbg!(origin_obj.radius() - vis_chk_obj_fok.radius());
    l /= l.norm();
    // Get the eclipsing body position in the same frame (allows us to set the origin to a null vector)
    let eclipsing_geoid_fok = cosm.celestial_state(eclipsing_geoid.id, origin_obj.dt.as_jde_et_days(), origin_obj.frame.id);

    let c = eclipsing_geoid_fok.radius(); // Note that we are using MINUS radius here (simplifies code)
    let r = dbg!(eclipsing_geoid.equatorial_radius);

    let discriminant_sqrt = (l.dot(&c)).powi(2) - c.dot(&c) - r.powi(2);
    if dbg!(discriminant_sqrt).is_sign_negative() {
        EclipseState::Visibilis
    } else if discriminant_sqrt.is_sign_positive() {
        EclipseState::Umbra
    } else {
        EclipseState::Penumbra
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use hifitime::Epoch;

    #[test]
    fn earth_eclipse() {
        let cosm = Cosm::from_xb("./de438s");
        let earth = dbg!(cosm.geoid_from_id(399));

        let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

        let sc1 = State::<Geoid>::from_keplerian(8000.0, 0.001, 0.1, 90.0, 75.0, 0.0, dt, earth);
        let sc2 = State::<Geoid>::from_keplerian(8001.0, 0.001, 0.1, 90.0, 75.0, 0.0, dt, earth);

        dbg!(eclipse_state(sc1, sc2, earth, &cosm));
    }
}
