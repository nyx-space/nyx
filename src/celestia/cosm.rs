extern crate bytes;
extern crate petgraph;
extern crate prost;
use self::bytes::IntoBuf;
use self::petgraph::algo::astar;
use self::petgraph::prelude::*;
use crate::hifitime::{Epoch, SECONDS_PER_DAY};
use celestia::cosm::prost::Message;
use celestia::frames::Frame;
use celestia::frames::*;
use celestia::state::State;
use celestia::xb::ephem_interp::StateData::{EqualStates, VarwindowStates};
use celestia::xb::{Ephemeris, EphemerisContainer, Identifier};
use celestia::xb::{Frame as FXBFrame, FrameContainer};
use celestia::SPEED_OF_LIGHT_KMS;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::fs::File;
pub use std::io::Error as IoError;
use std::io::Read;
use std::time::Instant;
use utils::rotv;

/// Enable or not light time correction for the computation of the celestial states
#[derive(Copy, Clone, Debug, PartialEq)]
pub enum LTCorr {
    /// No correction, i.e. assumes instantaneous propagation of photons
    None,
    /// Accounts for light-time correction. This is corresponds to CN in SPICE.
    LightTime,
    /// Accounts for light-time and stellar abberation where the solar system barycenter is the inertial frame. Corresponds to CN+S in SPICE.
    Abberation,
}

// Defines Cosm, from the Greek word for "world" or "universe".
#[derive(Clone)]
pub struct Cosm {
    ephemerides: HashMap<(i32, String), Ephemeris>,
    frames: HashMap<(i32, String), FXBFrame>,
    geoids: HashMap<(i32, String), Geoid>,
    exb_map: Graph<i32, u8, Undirected>,
}

impl fmt::Debug for Cosm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "Cosm with {} ephemerides", self.ephemerides.keys().len())
    }
}

#[derive(Debug)]
pub enum CosmError {
    ObjectIDNotFound(i32),
    ObjectNameNotFound(String),
    NoInterpolationData(i32),
    NoStateData(i32),
    /// No path was found to convert from the first center to the second
    DisjointFrameCenters(i32, i32),
    /// No path was found to convert from the first orientation to the second
    DisjointFrameOrientations(i32, i32),
}

impl fmt::Display for CosmError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "CosmError: {:?}", self)
    }
}

impl Error for CosmError {
    fn source(&self) -> Option<&(dyn Error + 'static)> {
        None
    }
}

impl Cosm {
    /// Builds a Cosm from the *XB files. Path should _not_ contain file extension. Panics if the files could not be loaded.
    pub fn from_xb(filename: &str) -> Cosm {
        Self::try_from_xb(filename)
            .unwrap_or_else(|_| panic!("could not open EXB file {}", filename))
    }

    /// Attempts to build a Cosm from the *XB files.
    pub fn try_from_xb(filename: &str) -> Result<Cosm, IoError> {
        let mut cosm = Cosm {
            ephemerides: HashMap::new(),
            frames: HashMap::new(),
            geoids: HashMap::new(),
            exb_map: Graph::new_undirected(),
        };

        // BUG: This adds the Sun, but the Sun should be added by the loop below (the de file DOES contain the Sun).
        // Further, when adding the SSB, I should probably _not_ add it to the Geoid, since it isn't one. That said, this will break things.
        // Hence, I should probably add as a "virtual geoid", and computing the approximate GM from some estimated mass of the solar system (which includes asteroids, etc.).
        // Solar System Barycenter
        let ss_mass = 1.0014; // https://en.wikipedia.org/w/index.php?title=Special:CiteThisPage&page=Solar_System&id=905437334
        let sun_gm = 132_712_440_041.939_38;
        cosm.geoids.insert(
            (0, "Solar System Barycenter".to_string()),
            Geoid::perfect_sphere(0, 0, 1, ss_mass * sun_gm),
        );

        cosm.exb_map.add_node(0); // Add the SSB

        let ephemerides = load_ephemeris(&(filename.to_string() + ".exb"))?;
        for ephem in &ephemerides {
            let id = ephem.id.as_ref().unwrap();
            let exb_tpl = (id.number, id.name.clone());

            cosm.ephemerides.insert(exb_tpl.clone(), ephem.clone());

            // Compute the exb_id and axb_id from the ref frame.
            let ref_frame_id = ephem.ref_frame.clone().unwrap().number;
            let exb_id = ref_frame_id % 100_000;
            let axb_id = ref_frame_id / 100_000;

            // Add this EXB to the map, and link it to its parent
            let this_node = cosm.exb_map.add_node(id.number);
            let parent_node = match cosm.exbid_to_map_idx(exb_id) {
                Ok(p) => p,
                Err(e) => panic!(e),
            };
            // All edges are of value 1
            cosm.exb_map.add_edge(parent_node, this_node, 1);

            // Build the Geoid frames -- assume all frames are geoids if they have a GM parameter

            // Ephemeris exists
            match ephem.parameters.get("GM") {
                Some(gm) => {
                    // It's a geoid, and we assume everything else is there
                    let flattening = match ephem.parameters.get("Flattening") {
                        Some(param) => param.value,
                        None => {
                            if id.name != "Moon" {
                                0.0
                            } else {
                                0.0012
                            }
                        }
                    };
                    let equatorial_radius = match ephem.parameters.get("Equatorial radius") {
                        Some(param) => param.value,
                        None => {
                            if id.name != "Moon" {
                                0.0
                            } else {
                                1738.1
                            }
                        }
                    };
                    let semi_major_radius = match ephem.parameters.get("Equatorial radius") {
                        Some(param) => {
                            if id.name != "Earth Barycenter" {
                                param.value
                            } else {
                                6378.1370
                            }
                        }
                        None => equatorial_radius, // assume spherical if unspecified
                    };
                    let geoid = Geoid {
                        id: id.number,
                        center_id: exb_id,
                        orientation_id: axb_id,
                        gm: gm.value,
                        flattening,
                        equatorial_radius,
                        semi_major_radius,
                    };
                    cosm.geoids.insert(exb_tpl, geoid);
                }
                None => {
                    match id.number {
                        10 => {
                            // Sun
                            let mut sun = Geoid::perfect_sphere(id.number, exb_id, axb_id, sun_gm);
                            // From https://iopscience.iop.org/article/10.1088/0004-637X/750/2/135
                            sun.equatorial_radius = 696_342.0;
                            cosm.geoids.insert(exb_tpl, sun);
                        }
                        _ => {
                            info!("no GM value for EXB ID {} (exb ID: {})", id.number, exb_id);
                        }
                    }
                }
            }
        }

        for frame in load_frames(&(filename.to_string() + ".fxb")) {
            let id = frame.id.clone().unwrap();
            cosm.frames.insert((id.number, id.name), frame.clone());
        }

        Ok(cosm)
    }

    fn exbid_to_map_idx(&self, id: i32) -> Result<NodeIndex, CosmError> {
        for (idx, node) in self.exb_map.raw_nodes().iter().enumerate() {
            if node.weight == id {
                return Ok(NodeIndex::new(idx));
            }
        }
        Err(CosmError::ObjectIDNotFound(id))
    }

    /// Returns the geoid from the loaded XB, if it is in there, else an error
    pub fn try_geoid_from_id(&self, id: i32) -> Result<Geoid, CosmError> {
        for ((geoid_id, _), geoid) in &self.geoids {
            if *geoid_id == id {
                return Ok(*geoid);
            }
        }
        Err(CosmError::ObjectIDNotFound(id))
    }

    /// Returns the geoid from the loaded XB, if it is in there, else panics!
    pub fn geoid_from_id(&self, id: i32) -> Geoid {
        self.try_geoid_from_id(id).unwrap()
    }

    /// Returns the list of loaded geoids
    pub fn geoids(&self) -> Vec<Geoid> {
        self.geoids.iter().map(|(_, g)| *g).collect()
    }

    /// Mutates the GM value for the provided geoid id. Panics if ID not found.
    pub fn mut_gm_for_geoid_id(&mut self, id: i32, gm: f64) {
        let mut key = None;
        for geoid_key in self.geoids.keys() {
            if geoid_key.0 == id {
                key = Some(geoid_key.clone());
            }
        }
        let mut geoid = self.geoids.get_mut(&key.unwrap()).unwrap();
        geoid.gm = gm;
    }

    pub fn geoid_from_name(&self, name: String) -> Result<Geoid, CosmError> {
        for ((_, geoid_name), geoid) in &self.geoids {
            if *geoid_name == name {
                return Ok(*geoid);
            }
        }
        Err(CosmError::ObjectNameNotFound(name))
    }

    pub fn exbid_from_id(&self, id: i32) -> Result<Identifier, CosmError> {
        for ((exb_id, _), ephem) in &self.ephemerides {
            if *exb_id == id {
                return Ok(ephem.id.clone().unwrap());
            }
        }
        Err(CosmError::ObjectIDNotFound(id))
    }

    /// Returns the celestial state as computed from a de4xx.{FXB,EXB} file in the original frame
    pub fn raw_celestial_state(&self, exb_id: i32, jde: f64) -> Result<State<Geoid>, CosmError> {
        let exb = self.exbid_from_id(exb_id)?;

        let ephem = self
            .ephemerides
            .get(&(exb.number, exb.name))
            .ok_or(CosmError::ObjectIDNotFound(exb.number))?;

        // Compute the position as per the algorithm from jplephem
        let interp = ephem
            .interpolator
            .as_ref()
            .ok_or(CosmError::NoInterpolationData(exb.number))?;

        // the DE file epochs are all in ET mod julian
        let start_mod_julian = ephem.start_epoch.as_ref().unwrap().value;
        let coefficient_count: usize = interp.position_degree as usize;

        let exb_states = match interp
            .state_data
            .as_ref()
            .ok_or(CosmError::NoStateData(exb.number))?
        {
            EqualStates(states) => states,
            VarwindowStates(_) => panic!("var window not yet supported by Cosm"),
        };

        let interval_length: f64 = exb_states.window_duration;

        let delta_jde = jde - start_mod_julian;
        let index_f = (delta_jde / interval_length).floor();
        let mut offset = delta_jde - index_f * interval_length;
        let mut index = index_f as usize;
        if index == exb_states.position.len() {
            index -= 1;
            offset = interval_length;
        }
        let pos_coeffs = &exb_states.position[index];

        let mut interp_t = vec![0.0; coefficient_count];
        let t1 = 2.0 * offset / interval_length - 1.0;
        interp_t[0] = 1.0;
        interp_t[1] = t1;
        for i in 2..coefficient_count {
            interp_t[i] = (2.0 * t1) * interp_t[i - 1] - interp_t[i - 2];
        }

        let mut interp_dt = vec![0.0; coefficient_count];
        interp_dt[0] = 0.0;
        interp_dt[1] = 1.0;
        interp_dt[2] = 2.0 * (2.0 * t1);
        for i in 3..coefficient_count {
            interp_dt[i] = (2.0 * t1) * interp_dt[i - 1] - interp_dt[i - 2]
                + interp_t[i - 1]
                + interp_t[i - 1];
        }
        for interp_i in &mut interp_dt {
            *interp_i *= 2.0 / interval_length;
        }

        let mut x = 0.0;
        let mut y = 0.0;
        let mut z = 0.0;
        let mut vx = 0.0;
        let mut vy = 0.0;
        let mut vz = 0.0;

        for (idx, pos_factor) in interp_t.iter().enumerate() {
            let vel_factor = interp_dt[idx];
            x += pos_factor * pos_coeffs.x[idx];
            y += pos_factor * pos_coeffs.y[idx];
            z += pos_factor * pos_coeffs.z[idx];
            vx += vel_factor * pos_coeffs.x[idx];
            vy += vel_factor * pos_coeffs.y[idx];
            vz += vel_factor * pos_coeffs.z[idx];
        }

        // Get the Geoid associated with the ephemeris frame
        let ref_frame_id = ephem.ref_frame.as_ref().unwrap().number;
        let ref_frame_exb_id = ref_frame_id % 100_000;
        let storage_geoid = self.geoid_from_id(ref_frame_exb_id);
        let dt = Epoch::from_jde_tai(jde);
        Ok(State::<Geoid>::from_cartesian(
            x,
            y,
            z,
            vx / SECONDS_PER_DAY,
            vy / SECONDS_PER_DAY,
            vz / SECONDS_PER_DAY,
            dt,
            storage_geoid,
        ))
    }

    /// Attempts to return the state of the celestial object of EXB ID `exb_id` (the target) at time `jde` `as_seen_from`
    ///
    /// The light time correction is based on SPICE's implementation: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html .
    /// Aberration computation is a conversion of the stelab function in SPICE, available here
    /// https://github.com/ChristopherRabotin/cspice/blob/26c72936fb7ff6f366803a1419b7cc3c61e0b6e5/src/cspice/stelab.c#L255
    pub fn try_celestial_state(
        &self,
        target_exb_id: i32,
        datetime: Epoch,
        as_seen_from_exb_id: i32,
        correction: LTCorr,
    ) -> Result<State<Geoid>, CosmError> {
        let as_seen_from = self.try_geoid_from_id(as_seen_from_exb_id)?;
        match correction {
            LTCorr::None => {
                let target_geoid = self.try_geoid_from_id(target_exb_id)?;
                let state = State::<Geoid>::from_cartesian(
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    0.0,
                    datetime,
                    target_geoid,
                );
                Ok(-self.try_frame_chg(&state, as_seen_from)?)
            }
            LTCorr::LightTime | LTCorr::Abberation => {
                // Get the geometric states as seen from SSB
                let obs =
                    self.try_celestial_state(as_seen_from_exb_id, datetime, 0, LTCorr::None)?;
                let mut tgt = self.try_celestial_state(target_exb_id, datetime, 0, LTCorr::None)?;
                // It will take less than three iterations to converge
                for _ in 0..3 {
                    // Compute the light time
                    let lt = (tgt - obs).rmag() / SPEED_OF_LIGHT_KMS;
                    // Compute the new target state
                    let mut lt_dt = datetime;
                    lt_dt.mut_sub_secs(lt);
                    tgt = self.celestial_state(target_exb_id, lt_dt, 0, LTCorr::None);
                }
                // Compute the correct state
                let mut state = State::<Geoid>::from_cartesian(
                    (tgt - obs).x,
                    (tgt - obs).y,
                    (tgt - obs).z,
                    (tgt - obs).vx,
                    (tgt - obs).vy,
                    (tgt - obs).vz,
                    datetime,
                    as_seen_from,
                );

                // Incluee the range-rate term in the velocity computation as explained in
                // https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/abcorr.html#Reception%20case
                let state_acc = state.velocity() / state.rmag();
                let dltdt = state.radius().dot(&state_acc) / SPEED_OF_LIGHT_KMS;

                state.vx = tgt.vx * (1.0 - dltdt) - obs.vx;
                state.vy = tgt.vy * (1.0 - dltdt) - obs.vy;
                state.vz = tgt.vz * (1.0 - dltdt) - obs.vz;

                if correction == LTCorr::Abberation {
                    // Get a unit vector that points in the direction of the object
                    let r_hat = state.r_hat();
                    // Get the velocity vector (of the observer) scaled with respect to the speed of light
                    let vbyc = obs.velocity() / SPEED_OF_LIGHT_KMS;
                    /* If the square of the length of the velocity vector is greater than or equal
                    to one, the speed of the observer is greater than or equal to the speed of light.
                    The observer speed is definitely out of range. */
                    if vbyc.dot(&vbyc) >= 1.0 {
                        warn!("observer is traveling faster than the speed of light");
                    } else {
                        let h_hat = r_hat.cross(&vbyc);
                        /* If the magnitude of the vector H is zero, the observer is moving along the line
                        of sight to the object, and no correction is required. Otherwise, rotate the
                        position of the object by phi radians about H to obtain the apparent position. */
                        if h_hat.norm() > std::f64::EPSILON {
                            let phi = h_hat.norm().asin();
                            let ab_pos = rotv(&state.radius(), &h_hat, phi);
                            state.x = ab_pos[0];
                            state.y = ab_pos[1];
                            state.z = ab_pos[2];
                        }
                    }
                }
                Ok(state)
            }
        }
    }

    /// Returns the state of the celestial object of EXB ID `exb_id` (the target) at time `jde` `as_seen_from`, or panics
    pub fn celestial_state(
        &self,
        target_exb_id: i32,
        datetime: Epoch,
        as_seen_from_exb_id: i32,
        correction: LTCorr,
    ) -> State<Geoid> {
        self.try_celestial_state(target_exb_id, datetime, as_seen_from_exb_id, correction)
            .unwrap()
    }

    /// Attempts to return the provided state in the provided frame.
    pub fn try_frame_chg(
        &self,
        state: &State<Geoid>,
        new_geoid: Geoid,
    ) -> Result<State<Geoid>, CosmError> {
        if state.frame.id() == new_geoid.id {
            return Ok(*state);
        }
        // Let's get the path between both both states.
        let path = if state.rmag() > 0.0 {
            // This is a state transformation, not a celestial position, so let's iterate backward
            // Not entirely sure why, but this works.
            self.intermediate_geoid(&new_geoid, &state.frame)?
        } else {
            self.intermediate_geoid(&state.frame, &new_geoid)?
        };

        let mut new_state = if state.frame.id() == 0 {
            // SSB, let's invert this
            -*state
        } else {
            *state
        };
        new_state.frame = new_geoid;
        let mut prev_frame_id = new_state.frame.id();
        let mut neg_needed = false;
        for body in path {
            if body.id() == 0 {
                neg_needed = !neg_needed;
                continue;
            }
            // This means the target or the origin is exactly this path.
            let mut next_state = self.raw_celestial_state(body.id(), state.dt.as_jde_et_days())?;
            if prev_frame_id == next_state.frame.id() || neg_needed {
                // Let's negate the next state prior to adding it.
                next_state = -next_state;
                neg_needed = true;
            }
            new_state = new_state + next_state;
            prev_frame_id = next_state.frame.id();
        }
        // If we started by opposing the state, let's do it again
        if state.frame.id() == 0 {
            new_state = -new_state;
        }
        Ok(new_state)
    }

    /// Return the provided state in the provided frame, or panics
    pub fn frame_chg(&self, state: &State<Geoid>, new_geoid: Geoid) -> State<Geoid> {
        self.try_frame_chg(state, new_geoid).unwrap()
    }

    /// Returns the conversion path from the target `from` as seen from `to`.
    pub fn intermediate_geoid(&self, from: &Geoid, to: &Geoid) -> Result<Vec<Geoid>, CosmError> {
        if from.orientation_id() != to.orientation_id() {
            unimplemented!("orientation changes are not yet implemented");
        }
        if from.id() == to.id() {
            // Same geoids, nothing to do
            return Ok(Vec::new());
        }

        // Check if the center of the target frame is the destination frame, or vice versa, so the path is simply one frame.
        if from.center_id() == to.id() {
            return Ok(vec![*from]);
        } else if to.center_id() == from.id() {
            return Ok(vec![*to]);
        }

        let start_idx = self.exbid_to_map_idx(to.id()).unwrap();
        let end_idx = self.exbid_to_map_idx(from.center_id()).unwrap();
        match astar(
            &self.exb_map,
            start_idx,
            |finish| finish == end_idx,
            |e| *e.weight(),
            |_| 0,
        ) {
            Some((weight, path)) => {
                // Build the path with the frames
                let mut f_path = Vec::new();
                let shared_centers = from.center_id == to.center_id;
                for idx in path {
                    let exb_id = self.exb_map[idx];
                    if !(shared_centers && to.center_id == exb_id || exb_id == 0) {
                        // Ignore going through SSB since it isn't a geoid
                        f_path.push(self.geoid_from_id(exb_id));
                    }
                }

                // Let's add the final geoid.
                f_path.push(*from);
                debug!(
                    "path from {:?} to {:?} is {:?} with cost {}",
                    from.id(),
                    to.id(),
                    f_path,
                    weight
                );
                Ok(f_path)
            }
            None => Err(CosmError::DisjointFrameCenters(
                from.center_id(),
                to.center_id(),
            )),
        }
    }
}

/// Loads the provided input_filename as an EXB
///
/// This function may panic!
pub fn load_ephemeris(input_filename: &str) -> Result<Vec<Ephemeris>, IoError> {
    let mut input_exb_buf = Vec::new();

    let mut f = File::open(input_filename)?;
    f.read_to_end(&mut input_exb_buf)
        .expect("something went wrong reading the file");

    if input_exb_buf.is_empty() {
        panic!("EXB file {} is empty (zero bytes read)", input_filename);
    }

    let decode_start = Instant::now();

    let ephcnt =
        EphemerisContainer::decode(input_exb_buf.into_buf()).expect("could not decode EXB");

    let ephemerides = ephcnt.ephemerides;
    let num_eph = ephemerides.len();
    if num_eph == 0 {
        panic!("no ephemerides found in EXB");
    }
    info!(
        "Loaded {} ephemerides in {} seconds.",
        num_eph,
        decode_start.elapsed().as_secs()
    );
    Ok(ephemerides)
}

/// Loads the provided input_filename as an FXB
///
/// This function may panic!
pub fn load_frames(input_filename: &str) -> Vec<FXBFrame> {
    let mut input_fxb_buf = Vec::new();

    File::open(input_filename)
        .unwrap_or_else(|_| panic!("could not open FXB file {}", input_filename))
        .read_to_end(&mut input_fxb_buf)
        .expect("something went wrong reading the file");

    if input_fxb_buf.is_empty() {
        panic!("FXB file {} is empty (zero bytes read)", input_filename);
    }

    let decode_start = Instant::now();

    let cnt = FrameContainer::decode(input_fxb_buf.into_buf()).expect("could not decode FXB");

    let frames = cnt.frames;
    let num_eph = frames.len();
    if num_eph == 0 {
        panic!("no frames found in FXB");
    }
    info!(
        "Loaded {} frames in {} seconds.",
        num_eph,
        decode_start.elapsed().as_secs()
    );
    frames
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestia::bodies;

    /// Tests direct transformations. Test cases generated via jplephem, hence the EPSILON precision.
    /// Note however that there is a difference between jplephem and spiceypy, cf.
    /// https://github.com/brandon-rhodes/python-jplephem/issues/33
    #[test]
    fn test_cosm_direct() {
        use std::f64::EPSILON;
        let cosm = Cosm::from_xb("./de438s");

        assert_eq!(
            cosm.intermediate_geoid(
                &cosm.geoid_from_id(bodies::EARTH_BARYCENTER),
                &cosm.geoid_from_id(bodies::EARTH_BARYCENTER),
            )
            .unwrap()
            .len(),
            0,
            "Conversions within Earth does not require any transformation"
        );

        let jde = Epoch::from_jde_et(2_452_312.5);
        let c = LTCorr::None;

        assert!(
            cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, bodies::EARTH_BARYCENTER, c)
                .rmag()
                < EPSILON
        );

        let out_state = cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, bodies::SSB, c);
        assert_eq!(out_state.frame.id(), bodies::SSB);
        assert!((out_state.x - -109_837_695.021_661_42).abs() < 1e-12);
        assert!((out_state.y - 89_798_622.194_651_56).abs() < 1e-12);
        assert!((out_state.z - 38_943_878.275_922_61).abs() < 1e-12);
        assert!((out_state.vx - -20.400_327_981_451_596).abs() < 1e-12);
        assert!((out_state.vy - -20.413_134_121_084_312).abs() < 1e-12);
        assert!((out_state.vz - -8.850_448_420_104_028).abs() < 1e-12);

        let out_state = cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, bodies::EARTH_MOON, c);
        assert_eq!(out_state.frame.id(), bodies::EARTH_MOON);
        assert!((out_state.x - 81_638.253_069_843_03).abs() < 1e-9);
        assert!((out_state.y - 345_462.617_249_631_9).abs() < 1e-9);
        assert!((out_state.z - 144_380.059_413_586_45).abs() < 1e-9);
        assert!((out_state.vx - -0.960_674_300_894_127_2).abs() < 1e-12);
        assert!((out_state.vy - 0.203_736_475_764_411_6).abs() < 1e-12);
        assert!((out_state.vz - 0.183_869_552_742_917_6).abs() < 1e-12);
        // Add the reverse test too
        let out_state = cosm.celestial_state(bodies::EARTH_MOON, jde, bodies::EARTH_BARYCENTER, c);
        assert_eq!(out_state.frame.id(), bodies::EARTH_BARYCENTER);
        assert!((out_state.x - -81_638.253_069_843_03).abs() < 1e-10);
        assert!((out_state.y - -345_462.617_249_631_9).abs() < 1e-10);
        assert!((out_state.z - -144_380.059_413_586_45).abs() < 1e-10);
        assert!((out_state.vx - 0.960_674_300_894_127_2).abs() < EPSILON);
        assert!((out_state.vy - -0.203_736_475_764_411_6).abs() < EPSILON);
        assert!((out_state.vz - -0.183_869_552_742_917_6).abs() < EPSILON);

        // The following test case comes from jplephem loaded with de438s.bsp
        let out_state = cosm.celestial_state(bodies::SUN, jde, bodies::SSB, c);
        assert_eq!(out_state.frame.id(), bodies::SSB);
        assert!((out_state.x - -182_936.040_274_732_14).abs() < EPSILON);
        assert!((out_state.y - -769_329.776_328_230_7).abs() < EPSILON);
        assert!((out_state.z - -321_490.795_782_183_1).abs() < EPSILON);
        assert!((out_state.vx - 0.014_716_178_620_115_785).abs() < EPSILON);
        assert!((out_state.vy - 0.001_242_263_392_603_425).abs() < EPSILON);
        assert!((out_state.vz - 0.000_134_043_776_253_089_48).abs() < EPSILON);

        let out_state = cosm.celestial_state(bodies::EARTH, jde, bodies::EARTH_BARYCENTER, c);
        assert_eq!(out_state.frame.id(), bodies::EARTH_BARYCENTER);
        assert!((out_state.x - 1_004.153_534_699_454_6).abs() < EPSILON);
        assert!((out_state.y - 4_249.202_979_894_305).abs() < EPSILON);
        assert!((out_state.z - 1_775.880_075_192_657_8).abs() < EPSILON);
        assert!((out_state.vx - -0.011_816_329_461_539_0).abs() < EPSILON);
        assert!((out_state.vy - 0.002_505_966_193_458_6).abs() < EPSILON);
        assert!((out_state.vz - 0.002_261_602_304_895_6).abs() < EPSILON);
    }

    #[test]
    fn test_cosm_indirect() {
        use std::f64::EPSILON;

        let jde = Epoch::from_gregorian_utc_at_midnight(2002, 2, 7);

        let cosm = Cosm::from_xb("./de438s");

        let ven2ear = cosm
            .intermediate_geoid(
                &cosm.geoid_from_id(bodies::VENUS_BARYCENTER),
                &cosm.geoid_from_id(bodies::EARTH_MOON),
            )
            .unwrap();
        assert_eq!(
            ven2ear.len(),
            3,
            "Venus -> (SSB) -> Earth Barycenter -> Earth Moon"
        );

        let c = LTCorr::None;
        // Error is sometimes up to 0.6 meters!
        // I think this is related to https://github.com/brandon-rhodes/python-jplephem/issues/33
        let tol_pos = 7e-4; // km
        let tol_vel = 5e-7; // km/s

        /*
        # Preceed all of the following python examples with
        >>> import spiceypy as sp
        >>> sp.furnsh('bsp/de438s.bsp')
        >>> et = 66312064.18493939
        */

        let ven2ear_state =
            cosm.celestial_state(bodies::VENUS_BARYCENTER, jde, bodies::EARTH_MOON, c);
        assert_eq!(ven2ear_state.frame.id(), bodies::EARTH_MOON);
        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(1, et, "J2000", "NONE", 301)[0]]
        ['2.0512621957200775e+08', '-1.3561254792308527e+08', '-6.5578399676151529e+07', '3.6051374278177832e+01', '4.8889024622170766e+01', '2.0702933800843084e+01']
        */
        assert!((ven2ear_state.x - 2.051_262_195_720_077_5e8).abs() < tol_pos);
        assert!((ven2ear_state.y - -1.356_125_479_230_852_7e8).abs() < tol_pos);
        assert!((ven2ear_state.z - -6.557_839_967_615_153e7).abs() < tol_pos);
        assert!((ven2ear_state.vx - 3.605_137_427_817_783e1).abs() < tol_vel);
        assert!((ven2ear_state.vy - 4.888_902_462_217_076_6e1).abs() < tol_vel);
        assert!((ven2ear_state.vz - 2.070_293_380_084_308_4e1).abs() < tol_vel);

        // Check that conversion via a center frame works
        let moon_from_emb =
            cosm.celestial_state(bodies::EARTH_MOON, jde, bodies::EARTH_BARYCENTER, c);
        // Check this state again, as in the direct test
        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(301, et, "J2000", "NONE", 3)[0]]
        ['-8.1576591043050896e+04', '-3.4547568914480874e+05', '-1.4439185901465410e+05', '9.6071184439702662e-01', '-2.0358322542180365e-01', '-1.8380551745739407e-01']
        */
        assert_eq!(moon_from_emb.frame.id(), bodies::EARTH_BARYCENTER);
        assert!((moon_from_emb.x - -8.157_659_104_305_09e4).abs() < tol_pos);
        assert!((moon_from_emb.y - -3.454_756_891_448_087_4e5).abs() < tol_pos);
        assert!((moon_from_emb.z - -1.443_918_590_146_541e5).abs() < tol_pos);
        assert!((moon_from_emb.vx - 9.607_118_443_970_266e-1).abs() < tol_vel);
        assert!((moon_from_emb.vy - -2.035_832_254_218_036_5e-1).abs() < tol_vel);
        assert!((moon_from_emb.vz - -1.838_055_174_573_940_7e-1).abs() < tol_vel);

        let earth_from_emb = cosm.celestial_state(bodies::EARTH, jde, bodies::EARTH_BARYCENTER, c);
        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(399, et, "J2000", "NONE", 3)[0]]
        ['1.0033950894874154e+03', '4.2493637646888546e+03', '1.7760252107225667e+03', '-1.1816791248014408e-02', '2.5040812085717632e-03', '2.2608146685133296e-03']
        */
        assert!((earth_from_emb.x - 1.003_395_089_487_415_4e3).abs() < tol_pos);
        assert!((earth_from_emb.y - 4.249_363_764_688_855e3).abs() < tol_pos);
        assert!((earth_from_emb.z - 1.776_025_210_722_566_7e3).abs() < tol_pos);
        assert!((earth_from_emb.vx - -1.181_679_124_801_440_8e-2).abs() < tol_vel);
        assert!((earth_from_emb.vy - 2.504_081_208_571_763e-3).abs() < tol_vel);
        assert!((earth_from_emb.vz - 2.260_814_668_513_329_6e-3).abs() < tol_vel);

        let moon_from_earth = cosm.celestial_state(bodies::EARTH_MOON, jde, bodies::EARTH, c);
        let earth_from_moon = cosm.celestial_state(bodies::EARTH, jde, bodies::EARTH_MOON, c);

        // assert_eq!(moon_from_emb - earth_from_emb, moon_from_earth);
        assert_eq!(earth_from_moon, -moon_from_earth);

        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(301, et, "J2000", "NONE", 399)[0]]
        ['-8.2579986132538310e+04', '-3.4972505290949758e+05', '-1.4616788422537665e+05', '9.7252863564504100e-01', '-2.0608730663037542e-01', '-1.8606633212590740e-01']
        */
        assert!((moon_from_earth.x - -8.257_998_613_253_831e4).abs() < tol_pos);
        assert!((moon_from_earth.y - -3.497_250_529_094_976e5).abs() < tol_pos);
        assert!((moon_from_earth.z - -1.461_678_842_253_766_5e5).abs() < tol_pos);
        assert!((moon_from_earth.vx - 9.725_286_356_450_41e-1).abs() < tol_vel);
        assert!((moon_from_earth.vy - -2.060_873_066_303_754_2e-1).abs() < tol_vel);
        assert!((moon_from_earth.vz - -1.860_663_321_259_074e-1).abs() < tol_vel);

        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(10, et, "J2000", "NONE", 399)[0]]
        ['1.0965506591533598e+08', '-9.0570891031525031e+07', '-3.9266577019474506e+07', '2.0426570124555724e+01', '2.0412112498804031e+01', '8.8484257849460111e+00']
        */
        let sun2ear_state = cosm.celestial_state(bodies::SUN, jde, bodies::EARTH, c);
        let emb_from_ssb = cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, bodies::SSB, c);
        let sun_from_ssb = cosm.celestial_state(bodies::SUN, jde, bodies::SSB, c);
        let delta_state = sun2ear_state + (-sun_from_ssb + emb_from_ssb + earth_from_emb);

        assert!(delta_state.radius().norm() < EPSILON);
        assert!(delta_state.velocity().norm() < EPSILON);

        assert!((sun2ear_state.x - 1.096_550_659_153_359_8e8).abs() < tol_pos);
        assert!((sun2ear_state.y - -9.057_089_103_152_503e7).abs() < tol_pos);
        assert!((sun2ear_state.z - -3.926_657_701_947_451e7).abs() < tol_pos);
        assert!((sun2ear_state.vx - 2.042_657_012_455_572_4e1).abs() < tol_vel);
        assert!((sun2ear_state.vy - 2.041_211_249_880_403e1).abs() < tol_vel);
        assert!((sun2ear_state.vz - 8.848_425_784_946_011).abs() < tol_vel);
        // And check the converse
        let ear2sun_state = cosm.celestial_state(bodies::EARTH, jde, bodies::SUN, c);
        let state_sum = ear2sun_state + sun2ear_state;
        assert!(state_sum.rmag() < EPSILON);
        assert!(state_sum.vmag() < EPSILON);
    }

    #[test]
    fn test_frame_change_earth2luna() {
        let cosm = Cosm::from_xb("./de438s");
        let earth = cosm.geoid_from_id(399);
        let luna = cosm.geoid_from_id(301);

        let jde = Epoch::from_jde_et(2_458_823.5);
        // From JPL HORIZONS
        let lro = State::<Geoid>::from_cartesian(
            4.017_685_334_718_784E5,
            2.642_441_356_763_487E4,
            -3.024_209_691_251_325E4,
            -6.168_920_999_978_097E-1,
            -6.678_258_076_726_339E-1,
            4.208_264_479_358_517E-1,
            jde,
            earth,
        );

        let lro_jpl = State::<Geoid>::from_cartesian(
            -3.692_315_939_257_387E2,
            8.329_785_181_291_3E1,
            -1.764_329_108_632_533E3,
            -5.729_048_963_901_611E-1,
            -1.558_441_873_361_044,
            4.456_498_438_933_088E-2,
            jde,
            luna,
        );

        let lro_wrt_moon = cosm.frame_chg(&lro, luna);
        println!("{}", lro_jpl);
        println!("{}", lro_wrt_moon);
        let lro_moon_earth_delta = lro_jpl - lro_wrt_moon;
        // Note that the passing conditions are large. JPL uses de431MX, but nyx uses de438s.
        assert!(lro_moon_earth_delta.rmag() < 1e-2);
        assert!(lro_moon_earth_delta.vmag() < 1e-5);
        // And the converse
        let lro_wrt_earth = cosm.frame_chg(&lro_wrt_moon, earth);
        assert!((lro_wrt_earth - lro).rmag() < std::f64::EPSILON);
        assert!((lro_wrt_earth - lro).vmag() < std::f64::EPSILON);
    }

    #[test]
    fn test_frame_change_ven2luna() {
        let cosm = Cosm::from_xb("./de438s");
        let luna = cosm.geoid_from_id(301);
        let venus = cosm.geoid_from_id(2);

        let jde = Epoch::from_jde_et(2_458_823.5);
        // From JPL HORIZONS
        let lro = State::<Geoid>::from_cartesian(
            -4.393_308_217_174_602E7,
            1.874_075_194_166_327E8,
            8.763_986_396_329_135E7,
            -5.054_051_490_556_286E1,
            -1.874_720_232_671_061E1,
            -6.518_342_268_306_54,
            jde,
            venus,
        );

        let lro_jpl = State::<Geoid>::from_cartesian(
            -3.692_315_939_257_387E2,
            8.329_785_181_291_3E1,
            -1.764_329_108_632_533E3,
            -5.729_048_963_901_611E-1,
            -1.558_441_873_361_044,
            4.456_498_438_933_088E-2,
            jde,
            luna,
        );

        let lro_wrt_moon = cosm.frame_chg(&lro, luna);
        println!("{}", lro_jpl);
        println!("{}", lro_wrt_moon);
        let lro_moon_earth_delta = lro_jpl - lro_wrt_moon;
        // Note that the passing conditions are very large. JPL uses de431MX, but nyx uses de438s.
        assert!(dbg!(lro_moon_earth_delta.rmag()) < 0.3);
        assert!(dbg!(lro_moon_earth_delta.vmag()) < 1e-5);
        // And the converse
        let lro_wrt_venus = cosm.frame_chg(&lro_wrt_moon, venus);
        assert!((lro_wrt_venus - lro).rmag() < std::f64::EPSILON);
        assert!((lro_wrt_venus - lro).vmag() < std::f64::EPSILON);
    }

    #[test]
    fn test_frame_change_ssb2luna() {
        let cosm = Cosm::from_xb("./de438s");
        let luna = cosm.geoid_from_id(301);
        let ssb = cosm.geoid_from_id(0);

        let jde = Epoch::from_jde_et(2_458_823.5);
        // From JPL HORIZONS
        let lro = State::<Geoid>::from_cartesian(
            4.227_396_973_787_854E7,
            1.305_852_533_250_192E8,
            5.657_002_470_685_254E7,
            -2.964_638_617_895_494E1,
            7.078_704_012_700_072,
            3.779_568_779_111_446,
            jde,
            ssb,
        );

        let lro_jpl = State::<Geoid>::from_cartesian(
            -3.692_315_939_257_387E2,
            8.329_785_181_291_3E1,
            -1.764_329_108_632_533E3,
            -5.729_048_963_901_611E-1,
            -1.558_441_873_361_044,
            4.456_498_438_933_088E-2,
            jde,
            luna,
        );

        let lro_wrt_moon = cosm.frame_chg(&lro, luna);
        println!("{}", lro_jpl);
        println!("{}", lro_wrt_moon);
        let lro_moon_earth_delta = lro_jpl - lro_wrt_moon;
        // Note that the passing conditions are very large. JPL uses de431MX, but nyx uses de438s.
        assert!(dbg!(lro_moon_earth_delta.rmag()) < 0.3);
        assert!(dbg!(lro_moon_earth_delta.vmag()) < 1e-5);
        // And the converse
        let lro_wrt_ssb = cosm.frame_chg(&lro_wrt_moon, ssb);
        assert!((lro_wrt_ssb - lro).rmag() < std::f64::EPSILON);
        assert!((lro_wrt_ssb - lro).vmag() < std::f64::EPSILON);
    }

    #[test]
    fn test_cosm_lt_corr() {
        let cosm = Cosm::from_xb("./de438s");

        let jde = Epoch::from_jde_et(2_452_312.500_742_881);

        let out_state = cosm.celestial_state(
            bodies::EARTH_BARYCENTER,
            jde,
            bodies::MARS_BARYCENTER,
            LTCorr::LightTime,
        );

        // Note that the following data comes from SPICE (via spiceypy).
        // There is currently a difference in computation for de438s: https://github.com/brandon-rhodes/python-jplephem/issues/33 .
        // However, in writing this test, I also checked the computed light time, which matches SPICE to 2.999058779096231e-10 seconds.
        assert!(dbg!(out_state.x - -2.577_185_470_734_315_8e8).abs() < 1e-3);
        assert!(dbg!(out_state.y - -5.814_057_247_686_307e7).abs() < 1e-3);
        assert!(dbg!(out_state.z - -2.493_960_187_215_911_6e7).abs() < 1e-3);
        assert!(dbg!(out_state.vx - -3.460_563_654_257_750_7).abs() < 1e-7);
        assert!(dbg!(out_state.vy - -3.698_207_386_702_523_5e1).abs() < 1e-7);
        assert!(dbg!(out_state.vz - -1.690_807_917_994_789_7e1).abs() < 1e-7);
    }

    #[test]
    fn test_cosm_aberration_corr() {
        let cosm = Cosm::from_xb("./de438s");

        let jde = Epoch::from_jde_et(2_452_312.500_742_881);

        let out_state = cosm.celestial_state(
            bodies::EARTH_BARYCENTER,
            jde,
            bodies::MARS_BARYCENTER,
            LTCorr::Abberation,
        );

        assert!((out_state.x - -2.577_231_712_700_484_4e8).abs() < 1e-3);
        assert!((out_state.y - -5.812_356_237_533_56e7).abs() < 1e-3);
        assert!((out_state.z - -2.493_146_410_521_204_8e7).abs() < 1e-3);
        // Reenable this test after #96 is implemented.
        dbg!(out_state.vx - -3.463_585_965_206_417);
        dbg!(out_state.vy - -3.698_169_177_803_263e1);
        dbg!(out_state.vz - -1.690_783_648_756_073e1);
    }
}
