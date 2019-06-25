extern crate bytes;
extern crate prost;
use self::bytes::IntoBuf;
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
}

#[derive(Debug)]
pub enum CosmError {
    ObjectIDNotFound(i32),
    ObjectNameNotFound(String),
    NoInterpolationData(i32, String),
    NoStateData(i32, String),
}

impl Cosm {
    /// Builds a Cosm from the *XB files. Path should _not_ contain file extension.
    pub fn from_xb(filename: &str) -> Cosm {
        let mut cosm = Cosm {
            ephemerides: HashMap::new(),
            frames: HashMap::new(),
            geoids: HashMap::new(),
        };

        let ephemerides = load_ephemeris(&(filename.to_string() + ".exb"));
        for ephem in &ephemerides {
            let id = ephem.id.clone().unwrap();
            let exb_tpl = (id.number, id.name.clone());

            cosm.ephemerides.insert(exb_tpl.clone(), ephem.clone());

            // Compute the exb_id and axb_id from the ref frame.
            let ref_frame_id = ephem.ref_frame.clone().unwrap().number;
            let exb_id = ref_frame_id % 100000;
            let axb_id = ref_frame_id / 100000;

            // Build the Geoid frames -- assume all frames are geoids if they have a GM parameter

            // Ephemeris exists
            if let Some(gm) = ephem.ephem_parameters.get("GM") {
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
                    center_id: exb_id,
                    orientation_id: axb_id,
                    gm: gm.value,
                    flattening,
                    equatorial_radius,
                    semi_major_radius,
                };
                cosm.geoids.insert(exb_tpl, geoid);
            }
        }

        for frame in load_frames(&(filename.to_string() + ".fxb")) {
            let id = frame.id.clone().unwrap();
            cosm.frames.insert((id.number, id.name), frame.clone());
        }

        // Solar System Barycenter
        cosm.geoids.insert(
            (0, "Solar System Barycenter".to_string()),
            Geoid::perfect_sphere(1, 1, 1.327_124_400_18e20),
        );

        cosm
    }

    pub fn geoid_from_id(&self, id: i32) -> Result<Geoid, CosmError> {
        for ((geoid_id, _), geoid) in &self.geoids {
            if *geoid_id == id {
                return Ok(geoid.clone());
            }
        }
        Err(CosmError::ObjectIDNotFound(id))
    }

    pub fn geoid_from_name(&self, name: String) -> Result<Geoid, CosmError> {
        for ((_, geoid_name), geoid) in &self.geoids {
            if *geoid_name == name {
                return Ok(geoid.clone());
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

    pub fn celestial_state<B: Frame>(&self, exb_id: i32, jde: f64, frame: B) -> Result<State<B>, CosmError> {
        let exb = self.exbid_from_id(exb_id)?;
        self.state(exb, jde, frame)
    }

    pub fn state<B: Frame>(&self, exb: EXBID, jde: f64, frame: B) -> Result<State<B>, CosmError> {
        let ephem = self
            .ephemerides
            .get(&(exb.number, exb.name))
            .ok_or(CosmError::ObjectIDNotFound(exb.number))?;

        // Compute the position
        // TODO: Maybe should this cache the previous ephemeris retrieved?
        let interp = ephem
            .clone()
            .interpolator
            .ok_or(CosmError::NoInterpolationData(exb.number, "exb.name".to_string()))?;

        let start_mod_julian: f64 = interp.start_mod_julian;
        let coefficient_count: usize = interp.position_degree as usize;

        let exb_states = match interp
            .state_data
            .ok_or(CosmError::NoStateData(exb.number, "exb.name".to_string()))?
        {
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

        // XXX: This uses the positions to compute the velocity.
        let mut interp_dt = vec![0.0; coefficient_count];
        interp_dt[0] = 0.0;
        interp_dt[1] = 1.0;
        interp_dt[2] = 2.0 * (2.0 * t1);
        for i in 3..coefficient_count {
            interp_dt[i] = (2.0 * t1) * interp_dt[i - 1] - interp_dt[i - 2] + interp_t[i - 1] + interp_t[i - 1];
        }
        for i in 0..interp_dt.len() {
            interp_dt[i] *= 2.0 / interval_length;
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

        let ref_frame = ephem.ref_frame.clone().unwrap();
        // Get the EXB of that frame
        let frame_center = self
            .frames
            .get(&(ref_frame.number, ref_frame.name))
            .unwrap()
            .exb_id
            .clone()
            .unwrap();
        // Get the Geoid associated with the ephemeris frame
        let storage_geoid = self.geoids.get(&(frame_center.number, frame_center.name)).unwrap();

        let dt = ModifiedJulian { days: jde - 2_400_000.5 };

        // We now have the state of the body in its storage frame.
        let storage_state = State::<Geoid>::from_cartesian(x, y, z, vx, vy, vz, dt, storage_geoid.clone());
        println!("{} {}", storage_state.rmag(), storage_state.vmag());
        if storage_geoid.orientation_id() != frame.orientation_id() {
            unimplemented!("orientation changes are not yet implemented");
        }
        if storage_geoid.center_id() != frame.center_id() {
            // Must convert the storage
            // TODO: Add a graph and find the path from storage frame to destination frame
        }

        // BUG: This does not perform any frame transformation
        Ok(State::<B>::from_cartesian(x, y, z, vx, vy, vz, dt, frame))
    }
}

/// Loads the provided input_filename as an EXB
///
/// This function may panic!
pub fn load_ephemeris(input_filename: &str) -> Vec<Ephemeris> {
    let mut input_exb_buf = Vec::new();

    File::open(input_filename)
        .expect(&format!("could not open EXB file {}", input_filename))
        .read_to_end(&mut input_exb_buf)
        .expect("something went wrong reading the file");

    if input_exb_buf.len() == 0 {
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
        .expect(&format!("could not open FXB file {}", input_filename))
        .read_to_end(&mut input_fxb_buf)
        .expect("something went wrong reading the file");

    if input_fxb_buf.len() == 0 {
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
    fn test_cosm() {
        let cosm = Cosm::from_xb("./de438s");
        for key in cosm.frames.keys() {
            println!("frame -- {:?} => {:?}", key, cosm.frames[&key]);
        }
        for ek in cosm.ephemerides.keys() {
            println!("ephem -- {:?} {:?}", ek, cosm.ephemerides[&ek].ephem_parameters);
        }
        for ek in cosm.geoids.keys() {
            println!("geoid -- {:?} {:?}", ek, cosm.geoids[&ek]);
        }

        let exb_id = EXBID {
            number: 3,
            name: "Earth Barycenter".to_string(),
        };

        let out_body = cosm.geoid_from_id(bodies::SUN).unwrap();

        let out_state = cosm.state(exb_id, 2474160.13175, out_body).unwrap();
        println!("{:?}", out_state);

        /*
        Expected data from jplephem
        (array([5.30527022e+07, 1.25344353e+08, 5.43293743e+07]), array([-2444703.8160139 ,   834536.49356688,   361669.07958066]))
        */
        // XXX: Why is the position that far off?!
        assert!((out_state.x - 5.30527022e+07).abs() < 1e-1);
        assert!((out_state.y - 1.25344353e+08).abs() < 1e-0);
        assert!((out_state.z - 5.43293743e+07).abs() < 1e-1);
        assert!((out_state.vx - -2444703.8160139).abs() < 1e-5);
        assert!((out_state.vy - 834536.49356688).abs() < 1e-5);
        assert!((out_state.vz - 361669.07958066).abs() < 1e-5);
    }

    #[test]
    fn test_cosm_transform() {
        let cosm = Cosm::from_xb("./de438s");

        let exb_id = EXBID {
            number: 3,
            name: "Earth Barycenter".to_string(),
        };

        let out_body = cosm.geoid_from_id(bodies::VENUS).unwrap();

        let out_state = cosm.state(exb_id, 2474160.13175, out_body).unwrap();
        println!("{:?}", out_state);

        /*
        Expected data from jplephem
        (array([5.30527022e+07, 1.25344353e+08, 5.43293743e+07]), array([-2444703.8160139 ,   834536.49356688,   361669.07958066]))
        */
        // XXX: Why is the position that far off?!
        assert!((out_state.x - 5.30527022e+07).abs() < 1e-1);
        assert!((out_state.y - 1.25344353e+08).abs() < 1e-0);
        assert!((out_state.z - 5.43293743e+07).abs() < 1e-1);
        assert!((out_state.vx - -2444703.8160139).abs() < 1e-5);
        assert!((out_state.vy - 834536.49356688).abs() < 1e-5);
        assert!((out_state.vz - 361669.07958066).abs() < 1e-5);
    }
}
