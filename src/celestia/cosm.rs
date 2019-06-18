extern crate bytes;
extern crate prost;
use self::bytes::IntoBuf;
use self::prost::Message;
use celestia::exb::interpolation::StateData::{EqualStates, VarwindowStates};
use celestia::exb::{Ephemeris, EphemerisContainer};
use celestia::frames::*;
use celestia::fxb::{Frame, FrameContainer};
use celestia::state::State;
use std::collections::HashMap;
use std::fs::File;
use std::io::Read;
use std::time::Instant;

// Defines Cosm, from the Greek word for "world" or "universe".
pub struct Cosm {
    ephemerides: HashMap<(i32, String), Ephemeris>,
    frames: HashMap<(i32, String), Frame>,
    geoids: HashMap<(i32, String), Geoid>, // TODO: Change to Graph (must confirm traverse)
}

#[derive(Debug)]
pub enum CosmError {
    ObjectNotFound(i32, String),
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
        for ephem in ephemerides {
            let id = ephem.id.clone().unwrap();
            cosm.ephemerides.insert((id.number, id.name), ephem);
        }

        for frame in load_frames(&(filename.to_string() + ".fxb")) {
            let id = frame.id.clone().unwrap();
            cosm.frames.insert((id.number, id.name), frame.clone());

            // Build the Geoid frames -- assume all frames are geoids if they have a GM parameter
            let exb_id = frame.exb_id.clone().unwrap();
            let exb_tpl = (exb_id.number, exb_id.name.clone());
            if let Some(ephem) = cosm.ephemerides.get(&exb_tpl) {
                // Ephemeris exists
                if let Some(gm) = ephem.ephem_parameters.get("GM") {
                    // It's a geoid, and we assume everything else is there
                    let geoid = Geoid {
                        frame_id: exb_id.clone(),
                        gm: gm.value,
                        flattening: ephem.ephem_parameters.get("Flattening").unwrap().value,
                        equatorial_radius: ephem.ephem_parameters.get("Equatorial radius").unwrap().value,
                        semi_major_radius: if exb_id.name == "EARTH BARYCENTER" {
                            println!("FOUND");
                            6378.1370
                        } else {
                            ephem.ephem_parameters.get("Equatorial radius").unwrap().value
                        },
                    };

                    cosm.geoids.insert(exb_tpl, geoid);
                }
            } else if exb_id.number == 0 {
                // Solar System Barycenter
                cosm.geoids
                    .insert(exb_tpl, Geoid::perfect_sphere(exb_id, 1.327_124_400_18e20));
            }
        }

        cosm
    }

    pub fn state<B: Body>(&self, exb: EXBID, jde: f64, frame: B) -> Result<State<B>, CosmError> {
        let ephem = self
            .ephemerides
            .get(&(exb.number, exb.name))
            .ok_or(CosmError::ObjectNotFound(exb.number, "exb.name".to_string()))?;

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

        // We now have the state of the body in its storage frame.
        let storage_state = State::<Geoid>::from_position_velocity(x, y, z, vx, vy, vz, storage_geoid.clone());
        println!("{} {}", storage_state.rmag(), storage_state.vmag());

        // BUG: This does not perform any frame transformation
        Ok(State::<B>::from_position_velocity(x, y, z, vx, vy, vz, frame))
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
pub fn load_frames(input_filename: &str) -> Vec<Frame> {
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

        let out_body = cosm.geoids[&(0, "Solar System Barycenter".to_string())].clone();

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
