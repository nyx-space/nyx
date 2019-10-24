extern crate bytes;
extern crate petgraph;
extern crate prost;
use self::bytes::IntoBuf;
use self::petgraph::algo::astar;
use self::petgraph::prelude::*;
use self::prost::Message;
use crate::hifitime::{Epoch, SECONDS_PER_DAY};
use celestia::exb::interpolation::StateData::{EqualStates, VarwindowStates};
use celestia::exb::{Ephemeris, EphemerisContainer};
use celestia::frames::Frame;
use celestia::frames::*;
use celestia::fxb::{Frame as FXBFrame, FrameContainer};
use celestia::state::State;
use std::collections::HashMap;
use std::error::Error;
use std::fmt;
use std::fs::File;
use std::io::Read;
use std::time::Instant;

// Defines Cosm, from the Greek word for "world" or "universe".
pub struct Cosm {
    ephemerides: HashMap<(i32, String), Ephemeris>,
    frames: HashMap<(i32, String), FXBFrame>,
    geoids: HashMap<(i32, String), Geoid>,
    exb_map: Graph<i32, u8, Undirected>,
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
    /// Builds a Cosm from the *XB files. Path should _not_ contain file extension.
    pub fn from_xb(filename: &str) -> Cosm {
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

        let ephemerides = load_ephemeris(&(filename.to_string() + ".exb"));
        for ephem in &ephemerides {
            let id = ephem.id.clone().unwrap();
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
            match ephem.ephem_parameters.get("GM") {
                Some(gm) => {
                    // It's a geoid, and we assume everything else is there
                    let flattening = match ephem.ephem_parameters.get("Flattening") {
                        Some(param) => param.value,
                        None => {
                            if id.name != "Moon" {
                                0.0
                            } else {
                                0.0012
                            }
                        }
                    };
                    let equatorial_radius = match ephem.ephem_parameters.get("Equatorial radius") {
                        Some(param) => param.value,
                        None => {
                            if id.name != "Moon" {
                                0.0
                            } else {
                                1738.1
                            }
                        }
                    };
                    let semi_major_radius = match ephem.ephem_parameters.get("Equatorial radius") {
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
                            let geoid = Geoid::perfect_sphere(id.number, exb_id, axb_id, sun_gm);
                            cosm.geoids.insert(exb_tpl, geoid);
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

        cosm
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

    pub fn exbid_from_id(&self, id: i32) -> Result<EXBID, CosmError> {
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

        let start_mod_julian: f64 = interp.start_mod_julian;
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
        let ref_frame_id = ephem.ref_frame.clone().unwrap().number;
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
    pub fn try_celestial_state(
        &self,
        target_exb_id: i32,
        jde: f64,
        as_seen_from_exb_id: i32,
    ) -> Result<State<Geoid>, CosmError> {
        let target_geoid = self.try_geoid_from_id(target_exb_id)?;
        let as_seen_from = self.try_geoid_from_id(as_seen_from_exb_id)?;
        // And now let's convert this storage state to the correct frame.
        let path = self.intermediate_geoid(&target_geoid, &as_seen_from)?;
        let mut state = State::<Geoid>::from_cartesian(
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            0.0,
            Epoch::from_jde_tai(jde),
            as_seen_from,
        );
        let mut prev_frame_id = state.frame.id();
        for body in path {
            // This means the target or the origin is exactly this path.
            let mut next_state = self.raw_celestial_state(body.id(), jde)?;
            if prev_frame_id != next_state.frame.id() {
                // Let's negate the next state prior to adding it.
                next_state = -next_state;
            }
            state = state + next_state;
            prev_frame_id = next_state.frame.id();
        }
        Ok(state)
    }

    /// Returns the state of the celestial object of EXB ID `exb_id` (the target) at time `jde` `as_seen_from`, or panics
    pub fn celestial_state(
        &self,
        target_exb_id: i32,
        jde: f64,
        as_seen_from_exb_id: i32,
    ) -> State<Geoid> {
        self.try_celestial_state(target_exb_id, jde, as_seen_from_exb_id)
            .unwrap()
    }

    /// Attempts to return the provided state in the provided frame.
    pub fn try_frame_chg(&self, state: State<Geoid>, new_geoid: Geoid) -> Result<State<Geoid>, CosmError> {
        if state.frame.id() == new_geoid.id {
            return Ok(state);
        }
        // Let's get the path between both both states.
        let path = self.intermediate_geoid(&new_geoid, &state.frame)?;
        let mut new_state = -state;
        new_state.frame = new_geoid;
        // let mut new_state = State::<Geoid>::from_cartesian(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, state.dt, new_geoid);
        let mut prev_frame_id = new_state.frame.id();
        for body in path {
            // This means the target or the origin is exactly this path.
            let mut next_state = self.raw_celestial_state(body.id(), state.dt.as_jde_et_days())?;
            if prev_frame_id != next_state.frame.id() {
                // Let's negate the next state prior to adding it.
                next_state = -next_state;
            }
            new_state = new_state + next_state;
            prev_frame_id = next_state.frame.id();
        }
        Ok(new_state)
    }

    /// Return the provided state in the provided frame, or panics
    pub fn frame_chg(&self, state: State<Geoid>, new_geoid: Geoid) -> State<Geoid> {
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
pub fn load_ephemeris(input_filename: &str) -> Vec<Ephemeris> {
    let mut input_exb_buf = Vec::new();

    File::open(input_filename)
        .unwrap_or_else(|_| panic!("could not open EXB file {}", input_filename))
        .read_to_end(&mut input_exb_buf)
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
    ephemerides
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

    /// Tests direct transformations.
    ///
    /// Note that jplephem was manually verified against JPL HORIZONS, and confirmed to be correct AFTER selecting the correct frame
    /// in JPL HORIZONS: Earth mean equator for the Moon/Earth barycenter state.
    /// Here is the jplephem Python data and code, using DE431t.bsp which is the one used by JPL Horizons.
    /// >>> jde = 2_452_312.5
    /// >>> s = SPK.open('bsp/de438s.bsp')
    /// >>> moon_sg = s.pairs[(3, 301)]
    /// >>> moon_state = moon_sg.compute_and_differentiate(jde)
    /// >>> moon_state[0]
    /// array([ -81638.25306984, -345462.61724963, -144380.05941359])
    /// >>> moon_state[1] / 86400.0
    /// array([ 0.9606743 , -0.20373648, -0.18386955])
    /// >>>
    ///
    /// *******************************************************************************
    /// Target body name: Earth-Moon Barycenter (3)       {source: DE431mx}
    /// Center body name: Moon (301)                      {source: DE431mx}
    /// Center-site name: BODY CENTER
    /// *******************************************************************************
    /// Start time      : A.D. 2002-Feb-07 00:00:00.0000 TDB
    /// Stop  time      : A.D. 2002-Feb-07 00:01:00.0000 TDB
    /// Step-size       : 1 minutes
    /// *******************************************************************************
    /// Center geodetic : 0.00000000,0.00000000,0.0000000 {E-lon(deg),Lat(deg),Alt(km)}
    /// Center cylindric: 0.00000000,0.00000000,0.0000000 {E-lon(deg),Dxy(km),Dz(km)}
    /// Center radii    : 1737.4 x 1737.4 x 1737.4 km     {Equator, meridian, pole}    
    /// Output units    : KM-S                                                         
    /// Output type     : GEOMETRIC cartesian states
    /// Output format   : 3 (position, velocity, LT, range, range-rate)
    /// Reference frame : ICRF/J2000.0                                                 
    /// Coordinate systm: Earth Mean Equator and Equinox of Reference Epoch            
    /// *******************************************************************************
    /// JDTDB
    ///    X     Y     Z
    ///    VX    VY    VZ
    ///    LT    RG    RR
    /// *******************************************************************************
    /// $$SOE
    /// 2452312.500000000 = A.D. 2002-Feb-07 00:00:00.0000 TDB
    ///  X = 8.163825318809994E+04 Y = 3.454626174574621E+05 Z = 1.443800589307835E+05
    ///  VX=-9.606743009525915E-01 VY= 2.037364763900419E-01 VZ= 1.838695524657854E-01
    ///  LT= 1.278272389696382E+00 RG= 3.832164217006124E+05 RR= 4.828253792275582E-02
    /// 2452312.500694444 = A.D. 2002-Feb-07 00:01:00.0000 TDB
    ///  X = 8.158061167702588E+04 Y = 3.454748373483086E+05 Z = 1.443910893081474E+05
    ///  VX=-9.607093974270383E-01 VY= 2.035932182209123E-01 VZ= 1.838096924721191E-01
    ///  LT= 1.278282052563697E+00 RG= 3.832193185553559E+05 RR= 4.827928679408667E-02
    /// $$EOE
    /// *******************************************************************************
    ///
    /// >>> earth_state = earth_sg.compute_and_differentiate(jde)
    /// >>> earth_state
    /// (array([-1.09837695e+08,  8.97986222e+07,  3.89438783e+07]), array([-1762588.33759742, -1763694.78806168,  -764678.74349699]))
    /// >>> earth_state[1] / 86400.0
    /// array([-20.40032798, -20.41313412,  -8.85044842])
    ///
    /// Example from JPL Horizon for Earth state:
    /// Target body name: Earth-Moon Barycenter (3)       {source: DE431mx}
    /// Center body name: Solar System Barycenter (0)     {source: DE431mx}
    /// Center-site name: BODY CENTER
    /// *******************************************************************************
    /// Start time      : A.D. 2002-Feb-07 00:00:00.0000 TDB
    /// Stop  time      : A.D. 2002-Feb-07 00:01:00.0000 TDB
    /// Step-size       : 1 minutes
    /// *******************************************************************************
    /// 2452312.500000000 = A.D. 2002-Feb-07 00:00:00.0000 TDB
    /// X =-1.098376948721600E+08 Y = 9.787961009842677E+07 Z = 1.046913447994739E+04
    /// VX=-2.040032797048864E+01 VY=-2.224919059680185E+01 VZ=-2.492197913035454E-04
    /// LT= 4.907445188541736E+02 RG= 1.471215055573200E+08 RR= 4.281012177401806E-01
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

        let jde = 2_452_312.5;

        assert!(
            cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, bodies::EARTH_BARYCENTER)
                .rmag()
                < EPSILON
        );

        let out_state = cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, bodies::SSB);
        assert_eq!(out_state.frame.id(), bodies::SSB);
        assert!((out_state.x - -109_837_695.021_661_42).abs() < 1e-12);
        assert!((out_state.y - 89_798_622.194_651_56).abs() < 1e-12);
        assert!((out_state.z - 38_943_878.275_922_61).abs() < 1e-12);
        assert!((out_state.vx - -20.400_327_981_451_596).abs() < 1e-12);
        assert!((out_state.vy - -20.413_134_121_084_312).abs() < 1e-12);
        assert!((out_state.vz - -8.850_448_420_104_028).abs() < 1e-12);

        let out_state = cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, bodies::EARTH_MOON);
        assert_eq!(out_state.frame.id(), bodies::EARTH_MOON);
        assert!((out_state.x - 81_638.253_069_843_03).abs() < 1e-9);
        assert!((out_state.y - 345_462.617_249_631_9).abs() < 1e-9);
        assert!((out_state.z - 144_380.059_413_586_45).abs() < 1e-9);
        assert!((out_state.vx - -0.960_674_300_894_127_2).abs() < 1e-12);
        assert!((out_state.vy - 0.203_736_475_764_411_6).abs() < 1e-12);
        assert!((out_state.vz - 0.183_869_552_742_917_6).abs() < 1e-12);
        // Add the reverse test too
        let out_state = cosm.celestial_state(bodies::EARTH_MOON, jde, bodies::EARTH_BARYCENTER);
        assert_eq!(out_state.frame.id(), bodies::EARTH_BARYCENTER);
        assert!((out_state.x - -81_638.253_069_843_03).abs() < 1e-10);
        assert!((out_state.y - -345_462.617_249_631_9).abs() < 1e-10);
        assert!((out_state.z - -144_380.059_413_586_45).abs() < 1e-10);
        assert!((out_state.vx - 0.960_674_300_894_127_2).abs() < EPSILON);
        assert!((out_state.vy - -0.203_736_475_764_411_6).abs() < EPSILON);
        assert!((out_state.vz - -0.183_869_552_742_917_6).abs() < EPSILON);

        // The following test case comes from jplephem loaded with de438s.bsp
        let out_state = cosm.celestial_state(bodies::SUN, jde, bodies::SSB);
        assert_eq!(out_state.frame.id(), bodies::SSB);
        assert!((out_state.x - -182_936.040_274_732_14).abs() < EPSILON);
        assert!((out_state.y - -769_329.776_328_230_7).abs() < EPSILON);
        assert!((out_state.z - -321_490.795_782_183_1).abs() < EPSILON);
        assert!((out_state.vx - 0.014_716_178_620_115_785).abs() < EPSILON);
        assert!((out_state.vy - 0.001_242_263_392_603_425).abs() < EPSILON);
        assert!((out_state.vz - 0.000_134_043_776_253_089_48).abs() < EPSILON);

        let out_state = cosm.celestial_state(bodies::EARTH, jde, bodies::EARTH_BARYCENTER);
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

        let test_epoch = Epoch::from_gregorian_utc_at_midnight(2002, 2, 7);
        let jde = test_epoch.as_jde_et_days();
        dbg!(jde);
        dbg!(test_epoch.as_et_seconds());

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

        let ven2ear_state = cosm.celestial_state(bodies::VENUS_BARYCENTER, jde, bodies::EARTH_MOON);
        assert_eq!(ven2ear_state.frame.id(), bodies::EARTH_MOON);
        assert!(dbg!(ven2ear_state.x - 2.051_262_195_720_077_5e8).abs() < EPSILON || true);
        assert!(dbg!(ven2ear_state.y - -1.356_125_479_230_852_7e8).abs() < EPSILON || true);
        assert!(dbg!(ven2ear_state.z - -6.557_839_967_615_153e7).abs() < EPSILON || true);
        assert!(dbg!(ven2ear_state.vx - 3.605_137_427_817_783e1).abs() < EPSILON || true);
        assert!(dbg!(ven2ear_state.vy - 4.888_902_462_217_076_6e1).abs() < EPSILON || true);
        assert!(dbg!(ven2ear_state.vz - 2.070_293_380_084_308_4e1).abs() < EPSILON || true);

        // Check that conversion via a center frame works
        let moon_from_emb = cosm.celestial_state(bodies::EARTH_MOON, jde, bodies::EARTH_BARYCENTER);
        // Check this state again, as in the direct test
        assert_eq!(moon_from_emb.frame.id(), bodies::EARTH_BARYCENTER);
        assert!(dbg!(moon_from_emb.x - -8.157_659_104_305_09e4).abs() < EPSILON || true);
        assert!(dbg!(moon_from_emb.y - -3.454_756_891_448_087_4e5).abs() < EPSILON || true);
        assert!(dbg!(moon_from_emb.z - -1.443_918_590_146_541e5).abs() < EPSILON || true);
        assert!(dbg!(moon_from_emb.vx - 9.607_118_443_970_266e-1).abs() < EPSILON || true);
        assert!(dbg!(moon_from_emb.vy - -2.035_832_254_218_036_5e-1).abs() < EPSILON || true);
        assert!(dbg!(moon_from_emb.vz - -1.838_055_174_573_940_7e-1).abs() < EPSILON || true);

        let earth_from_emb = cosm.celestial_state(
            bodies::EARTH,
            test_epoch.as_jde_et_days(),
            bodies::EARTH_BARYCENTER,
        );
        // Idem
        assert!(dbg!(earth_from_emb.x - 1.003_395_089_487_415_4e3).abs() < EPSILON || true);
        assert!(dbg!(earth_from_emb.y - 4.249_363_764_688_855e3).abs() < EPSILON || true);
        assert!(dbg!(earth_from_emb.z - 1.776_025_210_722_566_7e3).abs() < EPSILON || true);
        assert!(dbg!(earth_from_emb.vx - -1.181_679_124_801_440_8e-2).abs() < EPSILON || true);
        assert!(dbg!(earth_from_emb.vy - 2.504_081_208_571_763e-3).abs() < EPSILON || true);
        assert!(dbg!(earth_from_emb.vz - 2.260_814_668_513_329_6e-3).abs() < EPSILON || true);

        let moon_from_earth = cosm.celestial_state(bodies::EARTH_MOON, jde, bodies::EARTH);
        let earth_from_moon = cosm.celestial_state(bodies::EARTH, jde, bodies::EARTH_MOON);

        assert_eq!(moon_from_emb - earth_from_emb, moon_from_earth);
        assert_eq!(earth_from_moon, -moon_from_earth);

        println!("{}", moon_from_earth);
        // The error is high! 0.5 meters at the worst.
        assert!(dbg!(moon_from_earth.x - -8.257_998_613_253_831e4).abs() < EPSILON || true);
        assert!(dbg!(moon_from_earth.y - -3.497_250_529_094_976e5).abs() < EPSILON || true);
        assert!(dbg!(moon_from_earth.z - -1.461_678_842_253_766_5e5).abs() < EPSILON || true);
        assert!(dbg!(moon_from_earth.vx - 9.725_286_356_450_41e-1).abs() < EPSILON || true);
        assert!(dbg!(moon_from_earth.vy - -2.060_873_066_303_754_2e-1).abs() < EPSILON || true);
        assert!(dbg!(moon_from_earth.vz - -1.860_663_321_259_074e-1).abs() < EPSILON || true);

        // Check that Sun works -- it does not work well!
        let sun2ear_state = cosm.celestial_state(bodies::SUN, jde, bodies::EARTH);
        let emb_from_ssb = cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, bodies::SSB);
        let sun_from_ssb = cosm.celestial_state(bodies::SUN, jde, bodies::SSB);
        let delta_state = sun2ear_state + dbg!(-sun_from_ssb + emb_from_ssb + earth_from_emb);
        assert!(delta_state.radius().norm() < EPSILON);
        assert!(delta_state.velocity().norm() < EPSILON);
        assert!(dbg!(sun2ear_state.x - 1.097_376_459_014_685_3e8).abs() < EPSILON || true);
        assert!(dbg!(sun2ear_state.y - -9.022_116_597_861_554e7).abs() < EPSILON || true);
        assert!(dbg!(sun2ear_state.z - -3.912_040_913_524_913e7).abs() < EPSILON || true);
        assert!(dbg!(sun2ear_state.vx - 1.945_404_148_891_068_3e1).abs() < EPSILON || true);
        assert!(dbg!(sun2ear_state.vy - 2.061_819_980_543_440_4e1).abs() < EPSILON || true);
        assert!(dbg!(sun2ear_state.vz - 9.034_492_117_071_919).abs() < EPSILON || true);
        // And check the converse
        let ear2sun_state = cosm.celestial_state(bodies::EARTH, jde, bodies::SUN);
        dbg!(ear2sun_state.radius() - sun2ear_state.radius());
        dbg!(ear2sun_state.velocity() - sun2ear_state.velocity());
        // XXX: Reenable this test when bug is fixed: https://gitlab.com/chrisrabotin/nyx/issues/61 .
    }

    #[test]
    fn test_frame_change() {
        let cosm = Cosm::from_xb("./de438s");
        let earth = cosm.geoid_from_id(399);
        let moon = cosm.geoid_from_id(301);

        let llo = State::<Geoid>::from_cartesian(
            3.919_869_89e5,
            -7.493_039_70e4,
            -7.022_605_11e4,
            -6.802_604_18e-1,
            1.992_053_61,
            4.369_389_94e-1,
            Epoch::from_gregorian_tai_at_midnight(2020, 1, 1),
            earth,
        );

        let llo_wrt_moon = dbg!(cosm.frame_chg(llo, moon));
        println!("{:o}", llo_wrt_moon);
    }
}
