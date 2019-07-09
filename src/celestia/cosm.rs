extern crate bytes;
extern crate petgraph;
extern crate prost;
use self::bytes::IntoBuf;
use self::petgraph::algo::astar;
use self::petgraph::prelude::*;
use self::prost::Message;
use crate::hifitime::julian::ModifiedJulian;
use celestia::exb::interpolation::StateData::{EqualStates, VarwindowStates};
use celestia::exb::{Ephemeris, EphemerisContainer};
use celestia::frames::Frame;
use celestia::frames::*;
use celestia::fxb::{Frame as FXBFrame, FrameContainer};
use celestia::state::State;
use std::collections::HashMap;
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
        let sun_gm = 1.327_124_400_18e20;
        cosm.geoids.insert(
            (0, "Solar System Barycenter".to_string()),
            Geoid::perfect_sphere(0, 0, 1, ss_mass * sun_gm),
        );

        cosm.exb_map.add_node(0); // Add the SSB

        let ephemerides = load_ephemeris(&(filename.to_string() + ".exb"));
        for ephem in &ephemerides {
            let id = ephem.id.clone().unwrap();
            let exb_tpl = (id.number, id.name.clone());
            println!("{:?}", exb_tpl);

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
                                .geoid_from_id(3)
                                .expect("Earth Barycenter must be in EXB prior to Earth itself");
                            let moon_gm = cosm
                                .geoid_from_id(301)
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

    pub fn geoid_from_id(&self, id: i32) -> Result<Geoid, CosmError> {
        for ((geoid_id, _), geoid) in &self.geoids {
            if *geoid_id == id {
                return Ok(*geoid);
            }
        }
        Err(CosmError::ObjectIDNotFound(id))
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
        let interp = ephem.clone().interpolator.ok_or(CosmError::NoInterpolationData(exb.number))?;

        let start_mod_julian: f64 = interp.start_mod_julian;
        let coefficient_count: usize = interp.position_degree as usize;

        let exb_states = match interp.state_data.ok_or(CosmError::NoStateData(exb.number))? {
            EqualStates(states) => states.clone(),
            VarwindowStates(_) => panic!("var window not yet supported by Cosm"),
        };

        let interval_length: f64 = exb_states.window_duration;

        let delta_jde = jde - start_mod_julian;
        let index_f = (delta_jde / interval_length).round();
        let offset = delta_jde - index_f * interval_length;
        let index = index_f as usize;

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
        let storage_geoid = self.geoid_from_id(ephem.id.clone().unwrap().number).unwrap();
        let dt = ModifiedJulian { days: jde - 2_400_000.5 };
        Ok(State::<Geoid>::from_cartesian(x, y, z, vx, vy, vz, dt, storage_geoid))
    }

    /// Returns the state of the celestial object of EXB ID `exb_id` (the target) at time `jde` `as_seen_from`
    pub fn celestial_state(&self, target_exb_id: i32, jde: f64, as_seen_from_exb_id: i32) -> Result<State<Geoid>, CosmError> {
        let target_geoid = self.geoid_from_id(target_exb_id)?;
        let as_seen_from = self.geoid_from_id(as_seen_from_exb_id)?;
        // And now let's convert this storage state to the correct frame.
        let path = self.intermediate_geoid(&target_geoid, &as_seen_from)?;
        // Maybe make the path mutable and pop each item as they come?
        if path.len() == 1 {
            // This means the target or the origin is exactly this path.
            if path[0].center_id() == as_seen_from_exb_id {
                Ok(self.raw_celestial_state(target_exb_id, jde)?)
            } else {
                // Let's invert the state (sicne it's in the wrong frame), and fix the frame.
                let mut state = -self.raw_celestial_state(target_exb_id, jde)?;
                state.frame = as_seen_from;
                Ok(state)
            }
        } else {
            unimplemented!("convert celestial state to a different geoid");
        }
        /*
        for fno in 1..path.len() {
            let int_st = self.celestial_state(*path[fno - 1].id(), jde, &path[fno])?;
            state = state + int_st;
            println!("OK");
        }

        // And finally update the frame to the requested frame.
        state.frame = *origin;

        Ok(state)*/
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

        let start_idx = self.exbid_to_map_idx(from.id()).unwrap();
        let end_idx = self.exbid_to_map_idx(to.center_id()).unwrap();
        match astar(&self.exb_map, start_idx, |finish| finish == end_idx, |e| *e.weight(), |_| 0) {
            Some((weight, path)) => {
                // Build the path with the frames
                let mut f_path = Vec::new();
                for idx in path {
                    f_path.push(self.geoid_from_id(self.exb_map[idx]).unwrap());
                }
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

    println!("Decoding EXB (this may take a while for large files).");

    let decode_start = Instant::now();

    let ephcnt = EphemerisContainer::decode(input_exb_buf.into_buf()).expect("could not decode EXB");

    let ephemerides = ephcnt.ephemerides;
    let num_eph = ephemerides.len();
    if num_eph == 0 {
        panic!("no ephemerides found in EXB");
    }
    println!(
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

    println!("Decoding FXB (this may take a while for large files).");

    let decode_start = Instant::now();

    let cnt = FrameContainer::decode(input_fxb_buf.into_buf()).expect("could not decode FXB");

    let frames = cnt.frames;
    let num_eph = frames.len();
    if num_eph == 0 {
        panic!("no frames found in FXB");
    }
    println!("Loaded {} frames in {} seconds.", num_eph, decode_start.elapsed().as_secs());
    frames
}

#[cfg(test)]
mod tests {
    use super::*;
    use celestia::bodies;

    #[test]
    fn test_cosm_direct() {
        let cosm = Cosm::from_xb("./de438s");

        assert_eq!(
            cosm.intermediate_geoid(
                &cosm.geoid_from_id(bodies::EARTH_BARYCENTER).unwrap(),
                &cosm.geoid_from_id(bodies::EARTH_BARYCENTER).unwrap(),
            )
            .unwrap()
            .len(),
            0,
            "Conversions within Earth does not require any transformation"
        );

        // Using an odd time to check that the computation is correct for small differences.
        let out_state = cosm
            .celestial_state(bodies::EARTH_BARYCENTER, 2_474_160.159_786, bodies::SSB)
            .unwrap();

        /*
        Expected data from jplephem
        (array([5.29841561e+07, 1.25367735e+08, 5.43395074e+07]), array([-2445159.18641162,   833442.50611388,   361194.96441988]))
        */
        assert!((out_state.x - 5.298_415_61e+7).abs() < 1e-1);
        assert!((out_state.y - 1.253_677_35e+8).abs() < 1e0);
        assert!((out_state.z - 5.433_950_74e+7).abs() < 1e-1);
        assert!((out_state.vx - -2_445_159.186_411_62).abs() < 1e-3);
        assert!((out_state.vy - 833_442.506_113_88).abs() < 1e-3);
        assert!((out_state.vz - 361_194.964_419_88).abs() < 1e-3);

        let jde = 2_474_160.0;

        let out_state = cosm.celestial_state(bodies::EARTH_BARYCENTER, jde, bodies::SSB).unwrap();

        /*
        Expected data from jplephem
        (array([5.33746506e+07, 1.25234065e+08, 5.42815776e+07]), array([-2442556.02759262,   839674.57090544,   363895.83296089]))
        */
        assert!((out_state.x - 5.337_465_06e+7).abs() < 1e-1);
        assert!((out_state.y - 1.252_340_65e+8).abs() < 1e0);
        assert!((out_state.z - 5.428_157_76e+7).abs() < 1e-1);
        assert!((out_state.vx - -2_442_556.027_592_62).abs() < 1e-3);
        assert!((out_state.vy - 839_674.570_905_44).abs() < 1e-3);
        assert!((out_state.vz - 363_895.832_960_89).abs() < 1e-3);

        /*
        Expected data from jplephem on de438s.bsp
        >>> s = SPK.open('bsp/de438s.bsp')
        >>> moon_sg = s.pairs[(3, 301)]
        >>> moon_sg.compute_and_differentiate(2474160.0)
        (array([-223912.84465084,  251062.86640476,  139189.45802101]), array([-74764.97976908, -45782.16541702, -23779.8893784 ]))
        */

        // BUG: Why is this frame seemingly inverted?!
        let out_state = cosm
            .celestial_state(bodies::EARTH_BARYCENTER, jde, bodies::EARTH_MOON)
            .unwrap();
        assert_eq!(out_state.frame.id(), bodies::EARTH_MOON);
        assert!((out_state.x - -223_912.844_650_84).abs() < 1e-3);
        assert!((out_state.y - 251_062.866_404_76).abs() < 1e-3);
        assert!((out_state.z - 139_189.458_021_01).abs() < 1e-3);
        assert!((out_state.vx - -74_764.979_769_08).abs() < 1e-3);
        assert!((out_state.vy - -45_782.165_417_02).abs() < 1e-3);
        assert!((out_state.vz - -23_779.889_378_4).abs() < 1e-3);

        // Add the reverse test too
        let out_state = cosm
            .celestial_state(bodies::EARTH_MOON, jde, bodies::EARTH_BARYCENTER)
            .unwrap();
        assert_eq!(out_state.frame.id(), 3);
        assert!((out_state.x - 223_912.844_650_84).abs() < 1e-3);
        assert!((out_state.y - -251_062.866_404_76).abs() < 1e-3);
        assert!((out_state.z - -139_189.458_021_01).abs() < 1e-3);
        assert!((out_state.vx - 74_764.979_769_08).abs() < 1e-3);
        assert!((out_state.vy - 45_782.165_417_02).abs() < 1e-3);
        assert!((out_state.vz - 23_779.889_378_4).abs() < 1e-3);
    }

    #[test]
    fn test_cosm_indirect() {
        let cosm = Cosm::from_xb("./de438s");

        let ven2ear = cosm
            .intermediate_geoid(
                &cosm.geoid_from_id(bodies::VENUS_BARYCENTER).unwrap(),
                &cosm.geoid_from_id(bodies::EARTH_BARYCENTER).unwrap(),
            )
            .unwrap();
        assert_eq!(ven2ear.len(), 2, "Venus and Earth are in the same frame via the Sun");

        let path = cosm
            .intermediate_geoid(
                &cosm.geoid_from_id(bodies::VENUS).unwrap(),
                &cosm.geoid_from_id(bodies::EARTH_MOON).unwrap(),
            )
            .unwrap();
        assert_eq!(path.len(), 3, "Venus and Moon should convert via Sun and Earth");

        /*

        let out_state = cosm.celestial_state(bodies::VENUS, 2474160.0, &moon).unwrap();
        println!("{}", out_state);

        /*
        Expected data from jplephem (continued from above)
        >>> moon_state[0] + earth_state[0] - venus_state[0]
        array([-11582342.64316903,  45964259.67546433,  22672732.01650738])
        >>> moon_state[1] + earth_state[1] - venus_state[1]
        array([ -78659.10187959, -782169.48874783, -523505.15611691])
        */

        dbg!((out_state.x - -11582342.64316903).abs());
        dbg!((out_state.y - 45964259.67546433).abs());
        dbg!((out_state.z - 22672732.01650738).abs());
        dbg!((out_state.vx - -78659.10187959).abs());
        dbg!((out_state.vy - -782169.48874783).abs());
        dbg!((out_state.vz - -523505.15611691).abs());

        assert!((out_state.x - 1.158234274236687E+07).abs() < 1e-1);
        assert!((out_state.y - -5.119007862873630E+07).abs() < 1e-0);
        assert!((out_state.z - -2.518292393506091E+06).abs() < 1e-1);
        assert!((out_state.vx - 9.104062558220569E-01).abs() < 1e-5);
        assert!((out_state.vy - 1.071602860104991E+01).abs() < 1e-5);
        assert!((out_state.vz - 1.958072164387888E+00).abs() < 1e-5);*/

        // Let's simply check that the Venus state also works. No easy way to get that data, so no verification.
    }
}
