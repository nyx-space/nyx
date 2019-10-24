use super::{Cosm, Geoid, State};

/// Stores the eclipse state
#[derive(Clone, Copy, Debug, PartialEq)]
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
pub fn eclipse_state(
    origin_obj: &State<Geoid>,
    vis_chk_obj: &State<Geoid>,
    eclipsing_geoid: Geoid,
    cosm: &Cosm,
) -> EclipseState {
    if origin_obj == vis_chk_obj {
        return EclipseState::Visibilis;
    }
    // Convert the vis_chk_obj to the same frame as the origin.
    let vis_chk_obj_fok = &cosm.frame_chg(vis_chk_obj, origin_obj.frame);
    // Compute the unit direction between both
    let mut l = vis_chk_obj_fok.radius() - origin_obj.radius();
    dbg!(l);
    l /= l.norm();
    // Get the eclipsing body position in the same frame (allows us to set the origin to a null vector)
    let eclipsing_geoid_fok = cosm.celestial_state(
        eclipsing_geoid.id,
        origin_obj.dt.as_jde_et_days(),
        origin_obj.frame.id,
    );

    let omc = origin_obj.radius() - eclipsing_geoid_fok.radius();
    let r = eclipsing_geoid.equatorial_radius;

    let discriminant_sqrt = dbg!(l.dot(&omc)).powi(2) - dbg!(omc.dot(&omc) - r.powi(2));
    if dbg!(discriminant_sqrt) < 0.0 {
        EclipseState::Visibilis
    } else if discriminant_sqrt > 0.0 {
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
    fn trivial() {
        let cosm = Cosm::from_xb("./de438s");
        let mut earth = cosm.geoid_from_id(399);
        earth.equatorial_radius = 1.0;

        let dt = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);

        let sc1 = State::<Geoid>::from_cartesian(-1.0, -1.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let sc2 = State::<Geoid>::from_cartesian(1.0, 1.0, 0.0, 0.0, 0.0, 0.0, dt, earth);
        let sc3 = State::<Geoid>::from_cartesian(-2.0, 1.0, 0.0, 0.0, 0.0, 0.0, dt, earth);

        assert_eq!(eclipse_state(&sc1, &sc2, earth, &cosm), EclipseState::Umbra);
        assert_eq!(
            eclipse_state(&sc1, &sc3, earth, &cosm),
            EclipseState::Visibilis
        );
        assert_eq!(
            eclipse_state(&sc3, &sc2, earth, &cosm),
            EclipseState::Penumbra
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
        assert_eq!(eclipse_state(&sc1, &sc3, earth, &cosm), EclipseState::Umbra);

        // Nearly identical orbits in the same phasing
        assert_eq!(
            eclipse_state(&sc1, &sc2, earth, &cosm),
            EclipseState::Visibilis
        );
    }
}
