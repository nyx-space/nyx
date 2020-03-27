extern crate bytes;
extern crate petgraph;
extern crate prost;
use self::bytes::IntoBuf;
use self::petgraph::algo::astar;
use self::petgraph::prelude::*;
use super::cosm::prost::Message;
use super::frames::*;
use super::rotations::*;
use super::state::State;
use super::xb::ephem_interp::StateData::{EqualStates, VarwindowStates};
use super::xb::{Ephemeris, EphemerisContainer, Identifier as XbId};
use super::SPEED_OF_LIGHT_KMS;
use crate::hifitime::{Epoch, SECONDS_PER_DAY};
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::fs::File;
pub use std::io::Error as IoError;
use std::io::Read;
use std::time::Instant;
use utils::rotv;

/// Mass of the solar system from https://en.wikipedia.org/w/index.php?title=Special:CiteThisPage&page=Solar_System&id=905437334
pub const SS_MASS: f64 = 1.0014;
/// Mass of the Sun
pub const SUN_GM: f64 = 132_712_440_041.939_38;

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
pub struct Cosm {
    ephemerides: HashMap<i32, Ephemeris>,
    ephemerides_names: HashMap<String, i32>,
    frames: HashMap<String, FrameInfo>,
    // frames_names: HashMap<String, i32>,
    axb_names: HashMap<String, i32>,
    exb_map: Graph<i32, u8, Undirected>,
    axb_map: Graph<i32, Box<dyn ParentRotation>, Directed>,
    j2k_nidx: NodeIndex<u32>,
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
        let mut axb_map: Graph<i32, Box<dyn ParentRotation>, Directed> = Graph::new();
        let j2k_nidx = axb_map.add_node(0);

        let mut cosm = Cosm {
            ephemerides: HashMap::new(),
            ephemerides_names: HashMap::new(),
            frames: HashMap::new(),
            axb_names: HashMap::new(),
            exb_map: Graph::new_undirected(),
            axb_map,
            j2k_nidx,
        };

        // Add J2000 as the main AXB reference frame
        let ssb2k = FrameInfo::Celestial {
            axb_id: 0,
            exb_id: 0,
            gm: SS_MASS * SUN_GM,
        };

        cosm.axb_names.insert("J2000".to_owned(), 0);
        cosm.frames
            .insert("Solar System Barycenter J2k".to_owned(), ssb2k);

        match cosm.append_xb(filename) {
            None => Ok(cosm),
            Some(err) => Err(err),
        }

        // Now, let's append the IAU frames as defined in
        // Celest Mech Dyn Astr (2018) 130:22
        // https://doi.org/10.1007/s10569-017-9805-5
    }

    pub fn append_xb(&mut self, filename: &str) -> Option<IoError> {
        match load_ephemeris(&(filename.to_string() + ".exb")) {
            Err(e) => Some(e),
            Ok(ephemerides) => {
                for ephem in &ephemerides {
                    let id = ephem.id.as_ref().unwrap();
                    let exb_tpl = (id.number, id.name.clone());

                    self.ephemerides.insert(id.number, ephem.clone());
                    self.ephemerides_names.insert(id.name.clone(), id.number);

                    // TODO: Clone all of the frames and add their IAU fixed definitions

                    // Compute the exb_id.
                    let exb_id = ephem.ref_frame.clone().unwrap().number;

                    // Add this EXB to the map, and link it to its parent
                    let this_node = self.exb_map.add_node(id.number);
                    let parent_node = match self.exbid_to_map_idx(exb_id) {
                        Ok(p) => p,
                        Err(e) => panic!(e),
                    };
                    // All edges are of value 1
                    self.exb_map.add_edge(parent_node, this_node, 1);

                    // Build the Geoid frames -- assume all frames are geoids if they have a GM parameter

                    // Ephemeris exists
                    match ephem.parameters.get("GM") {
                        Some(gm) => {
                            // It's a geoid, and we assume everything else is there
                            let flattening = match ephem.parameters.get("Flattening") {
                                Some(param) => param.value,
                                None => {
                                    if id.name == "Moon" {
                                        0.0012
                                    } else {
                                        0.0
                                    }
                                }
                            };
                            let equatorial_radius = match ephem.parameters.get("Equatorial radius")
                            {
                                Some(param) => param.value,
                                None => {
                                    if id.name == "Moon" {
                                        1738.1
                                    } else {
                                        0.0
                                    }
                                }
                            };
                            let semi_major_radius = match ephem.parameters.get("Equatorial radius")
                            {
                                Some(param) => {
                                    if id.name == "Earth Barycenter" {
                                        6378.1370
                                    } else {
                                        param.value
                                    }
                                }
                                None => equatorial_radius, // assume spherical if unspecified
                            };

                            // Let's now build the J2000 version of this body
                            let obj = FrameInfo::Geoid {
                                axb_id: 0, // TODO: Get this from the EXB
                                exb_id,
                                gm: gm.value,
                                flattening,
                                equatorial_radius,
                                semi_major_radius,
                            };

                            self.frames.insert(id.name.clone(), obj);
                        }
                        None => {
                            match id.number {
                                10 => {
                                    // Build the Sun frame in J2000
                                    let sun2k = FrameInfo::Geoid {
                                        axb_id: 0,
                                        exb_id,
                                        gm: SUN_GM,
                                        flattening: 0.0,
                                        // From https://iopscience.iop.org/article/10.1088/0004-637X/750/2/135
                                        equatorial_radius: 696_342.0,
                                        semi_major_radius: 696_342.0,
                                    };

                                    self.frames.insert(id.name.clone(), sun2k);
                                    // And now in IAU body fixed

                                    // Define the IAU_SUN frame
                                    let sun2ssb_ctx = HashMap::new();
                                    let right_asc: meval::Expr = "289.13".parse().unwrap();
                                    let declin: meval::Expr = "63.87".parse().unwrap();
                                    let w_expr: meval::Expr =
                                        "84.176 + 14.18440000*d".parse().unwrap();
                                    let sun2ssb_rot = Euler3AxisDt::from_ra_dec_w(
                                        right_asc,
                                        declin,
                                        w_expr,
                                        &sun2ssb_ctx,
                                        AngleUnit::Degrees,
                                    );
                                    // Executively state that the IAU_SUN frame has the ID 10
                                    let sun_iau = FrameInfo::Geoid {
                                        axb_id: 10,
                                        exb_id,
                                        gm: SUN_GM,
                                        flattening: 0.0,
                                        // From https://iopscience.iop.org/article/10.1088/0004-637X/750/2/135
                                        equatorial_radius: 696_342.0,
                                        semi_major_radius: 696_342.0,
                                    };

                                    let sun_iau_node = self.axb_map.add_node(10);
                                    self.axb_names.insert("IAU SUN".to_owned(), 0);
                                    // And create the edge between the IAU SUN and J2k
                                    self.axb_map.add_edge(
                                        sun_iau_node,
                                        self.j2k_nidx,
                                        Box::new(sun2ssb_rot),
                                    );
                                    self.frames.insert("IAU SUN".to_owned(), sun_iau);
                                }
                                _ => {
                                    info!(
                                        "no GM value for EXB ID {} (exb ID: {})",
                                        id.number, exb_id
                                    );
                                }
                            }
                        }
                    }
                }

                None
            }
        }
    }

    fn exbid_to_map_idx(&self, id: i32) -> Result<NodeIndex, CosmError> {
        for (idx, node) in self.exb_map.raw_nodes().iter().enumerate() {
            if node.weight == id {
                return Ok(NodeIndex::new(idx));
            }
        }
        Err(CosmError::ObjectIDNotFound(id))
    }

    fn axbid_to_map_idx(&self, id: i32) -> Result<NodeIndex, CosmError> {
        for (idx, node) in self.axb_map.raw_nodes().iter().enumerate() {
            if node.weight == id {
                return Ok(NodeIndex::new(idx));
            }
        }
        Err(CosmError::ObjectIDNotFound(id))
    }

    pub fn try_frame(&self, name: String) -> Result<FrameInfo, CosmError> {
        match self.frames.get(&name) {
            Some(f) => Ok(*f),
            None => Err(CosmError::ObjectNameNotFound(name)),
        }
    }

    pub fn try_frame_by_id(&self, id: i32) -> Result<FrameInfo, CosmError> {
        match self.ephemerides.get(&id) {
            Some(e) => {
                // Now that we have the ephem for this ID, let's get the original frame
                Ok(self.frames[&e.id.as_ref().unwrap().name])
            }
            None => Err(CosmError::ObjectIDNotFound(id)),
        }
    }

    /// Returns the geoid from the loaded XB, if it is in there, else panics!
    pub fn frame(&self, name: String) -> FrameInfo {
        self.try_frame(name).unwrap()
    }

    /// Returns the geoid from the loaded XB, if it is in there, else panics!
    pub fn frame_by_id(&self, id: i32) -> FrameInfo {
        self.try_frame_by_id(id).unwrap()
    }

    /// Returns the list of loaded geoids
    pub fn frames(&self) -> Vec<FrameInfo> {
        self.frames.iter().map(|(_, g)| *g).collect()
    }

    /// Mutates the GM value for the frame used in the ephemeris only (not for every frame at this center)
    pub fn mut_gm_for_frame_id(&mut self, exb_id: i32, new_gm: f64) {
        match self.ephemerides.get(&exb_id) {
            Some(e) => {
                // Now that we have the ephem for this ID, let's get the original frame
                let name = e.id.as_ref().unwrap().name.clone(); //.to_string();
                self.mut_gm_for_frame(name, new_gm);
            }
            None => panic!("no frame ID {}", exb_id),
        }
    }

    /// Mutates the GM value for the provided geoid id. Panics if ID not found.
    pub fn mut_gm_for_frame(&mut self, name: String, new_gm: f64) {
        match self.frames.get_mut(&name) {
            Some(ref mut frame) => match frame {
                FrameInfo::Celestial { gm: mut gm, .. } => {
                    gm = new_gm;
                }
                FrameInfo::Geoid { gm: mut gm, .. } => {
                    gm = new_gm;
                }
                _ => panic!("frame ID {} does not have a GM", name),
            },
            None => panic!("no frame ID {}", name),
        }
    }

    /// Returns the celestial state as computed from a de4xx.{FXB,EXB} file in the original frame
    pub fn raw_celestial_state(&self, exb_id: i32, jde: f64) -> Result<State, CosmError> {
        let ephem = self
            .ephemerides
            .get(&exb_id)
            .ok_or(CosmError::ObjectIDNotFound(exb_id))?;

        // Compute the position as per the algorithm from jplephem
        let interp = ephem
            .interpolator
            .as_ref()
            .ok_or(CosmError::NoInterpolationData(exb_id))?;

        // the DE file epochs are all in ET mod julian
        let start_mod_julian = ephem.start_epoch.as_ref().unwrap().value;
        let coefficient_count: usize = interp.position_degree as usize;

        let exb_states = match interp
            .state_data
            .as_ref()
            .ok_or(CosmError::NoStateData(exb_id))?
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
        let storage_geoid = self.frame_by_id(ref_frame_exb_id);
        let dt = Epoch::from_jde_tai(jde);
        Ok(State::cartesian(
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
        frame: FrameInfo,
        correction: LTCorr,
    ) -> Result<State, CosmError> {
        match correction {
            LTCorr::None => {
                let target_frame = self.try_frame_by_id(target_exb_id)?;
                let state = State::cartesian(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, datetime, target_frame);
                Ok(-self.try_frame_chg(&state, frame)?)
            }
            LTCorr::LightTime | LTCorr::Abberation => {
                // Get the geometric states as seen from SSB
                let ssb2k = self.frame_by_id(0);
                let obs =
                    self.try_celestial_state(frame.exb_id(), datetime, ssb2k, LTCorr::None)?;
                let mut tgt =
                    self.try_celestial_state(target_exb_id, datetime, ssb2k, LTCorr::None)?;
                // It will take less than three iterations to converge
                for _ in 0..3 {
                    // Compute the light time
                    let lt = (tgt - obs).rmag() / SPEED_OF_LIGHT_KMS;
                    // Compute the new target state
                    let mut lt_dt = datetime;
                    lt_dt.mut_sub_secs(lt);
                    tgt = self.celestial_state(target_exb_id, lt_dt, ssb2k, LTCorr::None);
                }
                // Compute the correct state
                let mut state = State::cartesian(
                    (tgt - obs).x,
                    (tgt - obs).y,
                    (tgt - obs).z,
                    (tgt - obs).vx,
                    (tgt - obs).vy,
                    (tgt - obs).vz,
                    datetime,
                    frame,
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
        frame: FrameInfo,
        correction: LTCorr,
    ) -> State {
        self.try_celestial_state(target_exb_id, datetime, frame, correction)
            .unwrap()
    }

    /// Attempts to return the provided state in the provided frame.
    pub fn try_frame_chg(&self, state: &State, new_frame: FrameInfo) -> Result<State, CosmError> {
        if state.frame == new_frame {
            return Ok(*state);
        }
        // Let's get the path between both both states.
        let path = if state.rmag() > 0.0 {
            // This is a state transformation, not a celestial position, so let's iterate backward
            // Not entirely sure why, but this works.
            self.intermediate_geoid(&new_frame, &state.frame)?
        } else {
            self.intermediate_geoid(&state.frame, &new_frame)?
        };

        let mut new_state = if state.frame.exb_id() == 0 {
            // SSB, let's invert this
            -*state
        } else {
            *state
        };
        new_state.frame = new_frame;
        let mut prev_frame_exb_id = new_state.frame.exb_id();
        let mut neg_needed = false;
        for body in path {
            if body.exb_id() == 0 {
                neg_needed = !neg_needed;
                continue;
            }
            // This means the target or the origin is exactly this path.
            let mut next_state =
                self.raw_celestial_state(body.exb_id(), state.dt.as_jde_et_days())?;
            if prev_frame_exb_id == next_state.frame.exb_id() || neg_needed {
                // Let's negate the next state prior to adding it.
                next_state = -next_state;
                neg_needed = true;
            }
            new_state = new_state + next_state;
            prev_frame_exb_id = next_state.frame.exb_id();
        }
        // If we started by opposing the state, let's do it again
        if state.frame.exb_id() == 0 {
            new_state = -new_state;
        }
        Ok(new_state)
    }

    /// Return the provided state in the provided frame, or panics
    pub fn frame_chg(&self, state: &State, new_frame: FrameInfo) -> State {
        self.try_frame_chg(state, new_frame).unwrap()
    }

    /// Return the provided state in the provided frame, or panics
    pub fn frame_chg_by_id(&self, state: &State, new_frame: i32) -> State {
        let frame = self.frame_by_id(new_frame);
        self.try_frame_chg(state, frame).unwrap()
    }

    /// Returns the conversion path from the target `from` as seen from `to`.
    pub fn intermediate_geoid(
        &self,
        from: &FrameInfo,
        to: &FrameInfo,
    ) -> Result<Vec<FrameInfo>, CosmError> {
        // TODO: Add orientation computation

        if from == to {
            // Same geoids, nothing to do
            return Ok(Vec::new());
        }

        /*
                        // Get the frames themselves
                        let from_frame = self
                            .frames
                            .get(&from.exb_id())
                            .ok_or_else(|| CosmError::ObjectIDNotFound(from.exb_id()))?;

                        let to_frame = self
                            .frames
                            .get(&to.exb_id())
                            .ok_or_else(|| CosmError::ObjectIDNotFound(to.exb_id()))?;

                        // Check both frames have parents, if they don't and they aren't the same frame,
                        // then we can't do anything.

                        //self.frame_map[to_frame.id()];
                        //self.frame_map.neighbors_directed(a: NodeIndex<Ix>, dir: Direction)

                        if to_frame.parent.is_none() || from_frame.parent.is_none() {
                            return Err(CosmError::DisjointFrameCenters(from.exb_id(), to.exb_id()));
                        }

                        let (to_parent_frame, _to_parent_rot) = to_frame.parent.as_ref().unwrap();
                        let (from_parent_frame, _to_parent_rot) = from_frame.parent.as_ref().unwrap();

                // Check if the center of the target frame is the destination frame, or vice versa, so the path is simply one frame.
                if &from_parent_frame.info == to {
                    return Ok(vec![*from]);
                } else if &to_parent_frame.info == from {
                    return Ok(vec![*to]);
                }
        */
        let start_exb_idx = self.exbid_to_map_idx(to.exb_id()).unwrap();
        let end_exb_idx = self.exbid_to_map_idx(from.exb_id()).unwrap();

        let shared_centers = from.exb_id() == to.exb_id();
        match astar(
            &self.exb_map,
            start_exb_idx,
            |finish| finish == end_exb_idx,
            |e| *e.weight(),
            |_| 0,
        ) {
            Some((_, path)) => {
                // Build the path with the frames
                let mut f_path = Vec::new();
                for idx in path {
                    let exb_id = self.exb_map[idx];
                    if !(shared_centers && to.exb_id() == exb_id || exb_id == 0) {
                        // Ignore going through SSB since it isn't a geoid
                        f_path.push(self.frame_by_id(exb_id));
                    }
                }

                // Let's add the final geoid.
                f_path.push(*from);
                Ok(f_path)
            }
            None => Err(CosmError::DisjointFrameCenters(from.exb_id(), to.exb_id())),
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
                &cosm.frame_by_id(bodies::EARTH_BARYCENTER),
                &cosm.frame_by_id(bodies::EARTH_BARYCENTER),
            )
            .unwrap()
            .len(),
            0,
            "Conversions within Earth does not require any transformation"
        );

        let jde = Epoch::from_jde_et(2_452_312.5);
        let c = LTCorr::None;

        let earth_bary2k = cosm.frame_by_id(bodies::EARTH_BARYCENTER);
        let ssb2k = cosm.frame_by_id(bodies::SSB);
        let earth_moon2k = cosm.frame_by_id(bodies::EARTH_MOON);

        assert!(
            cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, earth_bary2k, c)
                .rmag()
                < EPSILON
        );

        let out_state = cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, ssb2k, c);
        assert_eq!(out_state.frame.exb_id(), bodies::SSB);
        assert!((out_state.x - -109_837_695.021_661_42).abs() < 1e-12);
        assert!((out_state.y - 89_798_622.194_651_56).abs() < 1e-12);
        assert!((out_state.z - 38_943_878.275_922_61).abs() < 1e-12);
        assert!((out_state.vx - -20.400_327_981_451_596).abs() < 1e-12);
        assert!((out_state.vy - -20.413_134_121_084_312).abs() < 1e-12);
        assert!((out_state.vz - -8.850_448_420_104_028).abs() < 1e-12);

        let out_state = cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, earth_moon2k, c);
        assert_eq!(out_state.frame.exb_id(), bodies::EARTH_MOON);
        assert!((out_state.x - 81_638.253_069_843_03).abs() < 1e-9);
        assert!((out_state.y - 345_462.617_249_631_9).abs() < 1e-9);
        assert!((out_state.z - 144_380.059_413_586_45).abs() < 1e-9);
        assert!((out_state.vx - -0.960_674_300_894_127_2).abs() < 1e-12);
        assert!((out_state.vy - 0.203_736_475_764_411_6).abs() < 1e-12);
        assert!((out_state.vz - 0.183_869_552_742_917_6).abs() < 1e-12);
        // Add the reverse test too
        let out_state = cosm.celestial_state(bodies::EARTH_MOON, jde, earth_bary2k, c);
        assert_eq!(out_state.frame.exb_id(), bodies::EARTH_BARYCENTER);
        assert!((out_state.x - -81_638.253_069_843_03).abs() < 1e-10);
        assert!((out_state.y - -345_462.617_249_631_9).abs() < 1e-10);
        assert!((out_state.z - -144_380.059_413_586_45).abs() < 1e-10);
        assert!((out_state.vx - 0.960_674_300_894_127_2).abs() < EPSILON);
        assert!((out_state.vy - -0.203_736_475_764_411_6).abs() < EPSILON);
        assert!((out_state.vz - -0.183_869_552_742_917_6).abs() < EPSILON);

        // The following test case comes from jplephem loaded with de438s.bsp
        let out_state = cosm.celestial_state(bodies::SUN, jde, ssb2k, c);
        assert_eq!(out_state.frame.exb_id(), bodies::SSB);
        assert!((out_state.x - -182_936.040_274_732_14).abs() < EPSILON);
        assert!((out_state.y - -769_329.776_328_230_7).abs() < EPSILON);
        assert!((out_state.z - -321_490.795_782_183_1).abs() < EPSILON);
        assert!((out_state.vx - 0.014_716_178_620_115_785).abs() < EPSILON);
        assert!((out_state.vy - 0.001_242_263_392_603_425).abs() < EPSILON);
        assert!((out_state.vz - 0.000_134_043_776_253_089_48).abs() < EPSILON);

        let out_state = cosm.celestial_state(bodies::EARTH, jde, earth_bary2k, c);
        assert_eq!(out_state.frame.exb_id(), bodies::EARTH_BARYCENTER);
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
                &cosm.frame_by_id(bodies::VENUS_BARYCENTER),
                &cosm.frame_by_id(bodies::EARTH_MOON),
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

        let earth_moon = cosm.frame_by_id(bodies::EARTH_MOON);
        let ven2ear_state = cosm.celestial_state(bodies::VENUS_BARYCENTER, jde, earth_moon, c);
        assert_eq!(ven2ear_state.frame.exb_id(), bodies::EARTH_MOON);
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
        let earth_bary = cosm.frame_by_id(bodies::EARTH_BARYCENTER);
        let moon_from_emb = cosm.celestial_state(bodies::EARTH_MOON, jde, earth_bary, c);
        // Check this state again, as in the direct test
        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(301, et, "J2000", "NONE", 3)[0]]
        ['-8.1576591043050896e+04', '-3.4547568914480874e+05', '-1.4439185901465410e+05', '9.6071184439702662e-01', '-2.0358322542180365e-01', '-1.8380551745739407e-01']
        */
        assert_eq!(moon_from_emb.frame, earth_bary);
        assert!((moon_from_emb.x - -8.157_659_104_305_09e4).abs() < tol_pos);
        assert!((moon_from_emb.y - -3.454_756_891_448_087_4e5).abs() < tol_pos);
        assert!((moon_from_emb.z - -1.443_918_590_146_541e5).abs() < tol_pos);
        assert!((moon_from_emb.vx - 9.607_118_443_970_266e-1).abs() < tol_vel);
        assert!((moon_from_emb.vy - -2.035_832_254_218_036_5e-1).abs() < tol_vel);
        assert!((moon_from_emb.vz - -1.838_055_174_573_940_7e-1).abs() < tol_vel);

        let earth_from_emb = cosm.celestial_state(bodies::EARTH, jde, earth_bary, c);
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

        let eme2k = cosm.frame_by_id(bodies::EARTH);
        let moon_from_earth = cosm.celestial_state(bodies::EARTH_MOON, jde, eme2k, c);
        let earth_from_moon = cosm.celestial_state(bodies::EARTH, jde, earth_moon, c);

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
        let sun2ear_state = cosm.celestial_state(bodies::SUN, jde, eme2k, c);
        let ssb_frame = cosm.frame_by_id(bodies::SSB);
        let emb_from_ssb = cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, ssb_frame, c);
        let sun_from_ssb = cosm.celestial_state(bodies::SUN, jde, ssb_frame, c);
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
        let sun2k = cosm.frame_by_id(bodies::SUN);
        let ear2sun_state = cosm.celestial_state(bodies::EARTH, jde, sun2k, c);
        let state_sum = ear2sun_state + sun2ear_state;
        assert!(state_sum.rmag() < EPSILON);
        assert!(state_sum.vmag() < EPSILON);
    }

    #[test]
    fn test_frame_change_earth2luna() {
        let cosm = Cosm::from_xb("./de438s");
        let earth = cosm.frame_by_id(399);
        let luna = cosm.frame_by_id(301);

        let jde = Epoch::from_jde_et(2_458_823.5);
        // From JPL HORIZONS
        let lro = State::cartesian(
            4.017_685_334_718_784E5,
            2.642_441_356_763_487E4,
            -3.024_209_691_251_325E4,
            -6.168_920_999_978_097E-1,
            -6.678_258_076_726_339E-1,
            4.208_264_479_358_517E-1,
            jde,
            earth,
        );

        let lro_jpl = State::cartesian(
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
        let luna = cosm.frame_by_id(301);
        let venus = cosm.frame_by_id(2);

        let jde = Epoch::from_jde_et(2_458_823.5);
        // From JPL HORIZONS
        let lro = State::cartesian(
            -4.393_308_217_174_602E7,
            1.874_075_194_166_327E8,
            8.763_986_396_329_135E7,
            -5.054_051_490_556_286E1,
            -1.874_720_232_671_061E1,
            -6.518_342_268_306_54,
            jde,
            venus,
        );

        let lro_jpl = State::cartesian(
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
        let luna = cosm.frame_by_id(301);
        let ssb = cosm.frame_by_id(0);

        let jde = Epoch::from_jde_et(2_458_823.5);
        // From JPL HORIZONS
        let lro = State::cartesian(
            4.227_396_973_787_854E7,
            1.305_852_533_250_192E8,
            5.657_002_470_685_254E7,
            -2.964_638_617_895_494E1,
            7.078_704_012_700_072,
            3.779_568_779_111_446,
            jde,
            ssb,
        );

        let lro_jpl = State::cartesian(
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

        let mars2k = cosm.frame_by_id(bodies::MARS_BARYCENTER);

        let out_state =
            cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, mars2k, LTCorr::LightTime);

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

        let mars2k = cosm.frame_by_id(bodies::MARS_BARYCENTER);

        let out_state =
            cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, mars2k, LTCorr::Abberation);

        assert!((out_state.x - -2.577_231_712_700_484_4e8).abs() < 1e-3);
        assert!((out_state.y - -5.812_356_237_533_56e7).abs() < 1e-3);
        assert!((out_state.z - -2.493_146_410_521_204_8e7).abs() < 1e-3);
        // Reenable this test after #96 is implemented.
        dbg!(out_state.vx - -3.463_585_965_206_417);
        dbg!(out_state.vy - -3.698_169_177_803_263e1);
        dbg!(out_state.vz - -1.690_783_648_756_073e1);
    }

    #[test]
    fn test_cosm_append() {
        let mut cosm = Cosm::from_xb("./de438s");
        // Let's load the same XB again to demonstrate idempotency.
        cosm.append_xb("./de438s");

        let jde = Epoch::from_jde_et(2_452_312.500_742_881);

        let mars2k = cosm.frame_by_id(bodies::MARS_BARYCENTER);

        let out_state =
            cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, mars2k, LTCorr::Abberation);

        assert!((out_state.x - -2.577_231_712_700_484_4e8).abs() < 1e-3);
        assert!((out_state.y - -5.812_356_237_533_56e7).abs() < 1e-3);
        assert!((out_state.z - -2.493_146_410_521_204_8e7).abs() < 1e-3);
        // Reenable this test after #96 is implemented.
        dbg!(out_state.vx - -3.463_585_965_206_417);
        dbg!(out_state.vy - -3.698_169_177_803_263e1);
        dbg!(out_state.vz - -1.690_783_648_756_073e1);
    }

    #[test]
    fn test_rotations() {
        let cosm = Cosm::from_xb("./de438s");
        let earth = cosm.frame_by_id(301);
        dbg!(earth);
    }

    fn test_new_frames_usage() {
        /*
        // Later there will be a single XB which also contains the frames
        let cosm = Cosm::from_xb("./de438s.exb");
        // Get the frame we want
        let eme2k = cosm.frame_by_name("EME2000");
        // Define the state
        let state = State::keplerian(sma, ecc, inc, raan, aop, ta, dt, eme2k);
        */
    }
}
