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
                        None => 0.0,
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
                        399 => {
                            // Compute the Earth GM by subtracting the Moon GM from the Earth System GM
                            let earth_sys = cosm
                                .try_geoid_from_id(3)
                                .expect("Earth Barycenter must be in EXB prior to Earth itself");
                            let moon_gm = cosm
                                .try_geoid_from_id(301)
                                .expect("Earth Moon must be in EXB prior to Earth itself")
                                .gm;
                            let earth = Geoid {
                                id: id.number,
                                center_id: exb_id,
                                orientation_id: axb_id,
                                gm: earth_sys.gm - moon_gm,
                                flattening: earth_sys.flattening,
                                equatorial_radius: earth_sys.equatorial_radius,
                                semi_major_radius: earth_sys.semi_major_radius,
                            };
                            cosm.geoids.insert(exb_tpl, earth);
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

        let exb_states = match interp.state_data.as_ref().ok_or(CosmError::NoStateData(exb.number))? {
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
            interp_dt[i] = (2.0 * t1) * interp_dt[i - 1] - interp_dt[i - 2] + interp_t[i - 1] + interp_t[i - 1];
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

    /// Returns the state of the celestial object of EXB ID `exb_id` (the target) at time `jde` `as_seen_from`
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
        let mut state = State::<Geoid>::from_cartesian(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, Epoch::from_jde_tai(jde), as_seen_from);
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

    pub fn celestial_state(&self, target_exb_id: i32, jde: f64, as_seen_from_exb_id: i32) -> State<Geoid> {
        self.try_celestial_state(target_exb_id, jde, as_seen_from_exb_id).unwrap()
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
        match astar(&self.exb_map, start_idx, |finish| finish == end_idx, |e| *e.weight(), |_| 0) {
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
            None => Err(CosmError::DisjointFrameCenters(from.center_id(), to.center_id())),
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

    let ephcnt = EphemerisContainer::decode(input_exb_buf.into_buf()).expect("could not decode EXB");

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
    info!("Loaded {} frames in {} seconds.", num_eph, decode_start.elapsed().as_secs());
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
        assert!((out_state.x - -81_638.253_069_843_03).abs() < 1e-9);
        assert!((out_state.y - -345_462.617_249_631_9).abs() < 1e-9);
        assert!((out_state.z - -144_380.059_413_586_45).abs() < 1e-9);
        assert!((out_state.vx - 0.960_674_300_894_127_2).abs() < 1e-12);
        assert!((out_state.vy - -0.203_736_475_764_411_6).abs() < 1e-12);
        assert!((out_state.vz - -0.183_869_552_742_917_6).abs() < 1e-12);
    }

    #[test]
    fn test_cosm_indirect() {
        let jde = 2_452_312.5;

        let cosm = Cosm::from_xb("./de438s");

        let ven2ear = cosm
            .intermediate_geoid(
                &cosm.geoid_from_id(bodies::VENUS_BARYCENTER),
                &cosm.geoid_from_id(bodies::EARTH_MOON),
            )
            .unwrap();
        assert_eq!(ven2ear.len(), 3, "Venus -> (SSB) -> Earth Barycenter -> Earth Moon");

        let ven2ear_state = cosm.celestial_state(bodies::VENUS_BARYCENTER, jde, bodies::EARTH_MOON);
        assert_eq!(ven2ear_state.frame.id(), bodies::EARTH_MOON);
        assert!((ven2ear_state.x - 205_123_905.586_005_96).abs() < 1e-7);
        assert!((ven2ear_state.y - -135_615_685.849_769_38).abs() < 1e-7);
        assert!((ven2ear_state.z - -65_579_728.485_815_55).abs() < 1e-7);
        assert!((ven2ear_state.vx - 36.052_331_787_059_61).abs() < 1e-12);
        assert!((ven2ear_state.vy - 48.888_638_291_647_33).abs() < 1e-12);
        assert!((ven2ear_state.vz - 20.702_719_193_464_77).abs() < 1e-12);

        // Check that conversion via a center frame works
        let moon_from_emb = cosm.celestial_state(bodies::EARTH_MOON, jde, bodies::EARTH_BARYCENTER);

        let earth_from_emb = cosm.celestial_state(bodies::EARTH, jde, bodies::EARTH_BARYCENTER);

        let moon_from_earth = cosm.celestial_state(bodies::EARTH_MOON, jde, bodies::EARTH);
        let earth_from_moon = cosm.celestial_state(bodies::EARTH, jde, bodies::EARTH_MOON);

        assert_eq!(moon_from_emb - earth_from_emb, moon_from_earth);
        assert_eq!(earth_from_moon, -moon_from_earth);

        // Check that Moon <-> Earth state is correct from JPL Horizons
        println!("{}", moon_from_earth);
        assert!(dbg!(moon_from_earth.x - -8.264240671516322E4).abs() < 1e-1);
        assert!(dbg!(moon_from_earth.y - -3.497118204014440E5).abs() < 1e-1);
        assert!(dbg!(moon_from_earth.z - -1.461559389839604E5).abs() < 1e-1);
        assert!(dbg!(moon_from_earth.vx - 9.724906303078747E-1).abs() < 1e-7);
        assert!(dbg!(moon_from_earth.vy - -2.062424425685089E-1).abs() < 1e-7);
        assert!(dbg!(moon_from_earth.vz - -1.861311547467976E-1).abs() < 1e-7);

        // Check that Sun works
        let sun2ear_state = cosm.celestial_state(bodies::SUN, jde, bodies::EARTH);
        println!("{}", sun2ear_state);
        assert!((sun2ear_state.x - 1.096_537_548_791_171E8).abs() < 1e-1);
        assert!((sun2ear_state.y - -9.057_220_115_376_98E7).abs() < 1e-1);
        assert!((sun2ear_state.z - -3.926_714_485_809_306E7).abs() < 1e-1);
        assert!((sun2ear_state.vx - 2.042_686_047_796_898E1).abs() < 1e-7);
        assert!((sun2ear_state.vy - 2.041_187_043_287_114E1).abs() < 1e-7);
        assert!((sun2ear_state.vz - 8.848_320_853_924_715).abs() < 1e-7);
        // And check the converse
        let ear2sun_state = cosm.celestial_state(bodies::EARTH, jde, bodies::SUN);
        dbg!(ear2sun_state.radius() - sun2ear_state.radius());
        dbg!(ear2sun_state.velocity() - sun2ear_state.velocity());
        // XXX: Reenable this test when bug is fixed: https://gitlab.com/chrisrabotin/nyx/issues/61 .
    }
}
