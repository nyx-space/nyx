/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

extern crate bytes;
extern crate meval;
extern crate prost;
extern crate rust_embed;
extern crate toml;

use self::meval::Expr;
use self::rust_embed::RustEmbed;
use super::frames::*;
use super::orbit::Orbit;
use super::rotations::*;
use super::xb::ephem_interp::StateData::{EqualStates, VarwindowStates};
use super::xb::{Ephemeris, Xb};
use super::SPEED_OF_LIGHT_KMS;
use crate::errors::NyxError;
use crate::hifitime::{Epoch, Unit, SECONDS_PER_DAY};
use crate::io::frame_serde;
use crate::na::{Matrix3, Matrix6};
use crate::utils::{capitalize, dcm_finite_differencing, rotv};
use std::collections::HashMap;
use std::fmt;
pub use std::io::{Error as IoError, ErrorKind as IoErrorKind};
use std::str::FromStr;
use std::sync::Arc;

#[derive(RustEmbed)]
#[folder = "data/embed/"]
struct EmbeddedAsset;

/// Mass of the solar system from https://en.wikipedia.org/w/index.php?title=Special:CiteThisPage&page=Solar_System&id=905437334
pub const SS_MASS: f64 = 1.0014;
/// GM of the Sun in km^3/s^2
pub const SUN_GM: f64 = 132_712_440_041.939_38;

/// Enable or not light time correction for the computation of the celestial states
#[derive(Copy, Clone, Debug, PartialEq)]
#[allow(clippy::upper_case_acronyms)]
pub enum LightTimeCalc {
    /// No correction, i.e. assumes instantaneous propagation of photons
    None,
    /// Accounts for light-time correction. This is corresponds to CN in SPICE.
    LightTime,
    /// Accounts for light-time and stellar abberation where the solar system barycenter is the inertial frame. Corresponds to CN+S in SPICE.
    Abberation,
}

#[derive(Debug)]
pub struct FrameTree {
    name: String,
    frame: Frame,
    // If None, and has a parent (check frame.frame_path()), then rotation is I33 (and therefore no need to compute it)
    parent_rotation: Option<Box<dyn ParentRotation>>,
    children: Vec<FrameTree>,
}

impl FrameTree {
    /// Seek an ephemeris from its celestial name (e.g. Earth Moon Barycenter)
    fn frame_seek_by_name(
        name: &str,
        cur_path: &mut Vec<usize>,
        f: &FrameTree,
    ) -> Result<Vec<usize>, NyxError> {
        if f.name == name {
            Ok(cur_path.to_vec())
        } else if f.children.is_empty() {
            Err(NyxError::ObjectNotFound(name.to_string()))
        } else {
            for (cno, child) in f.children.iter().enumerate() {
                let mut this_path = cur_path.clone();
                this_path.push(cno);
                let child_attempt = Self::frame_seek_by_name(name, &mut this_path, child);
                if let Ok(found_path) = child_attempt {
                    return Ok(found_path);
                }
            }
            // Could not find name in iteration, fail
            Err(NyxError::ObjectNotFound(name.to_string()))
        }
    }
}

// Defines Cosm, from the Greek word for "world" or "universe".
pub struct Cosm {
    pub xb: Xb,
    pub frame_root: FrameTree,
    // Maps the ephemeris path to the frame root path (remove this with the upcoming xb file)
    ephem2frame_map: HashMap<Vec<usize>, Vec<usize>>,
}

impl fmt::Debug for Cosm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(
            f,
            "Cosm with `{}` as ephemerides root",
            match &self.xb.ephemeris_root {
                Some(r) => r.name.clone(),
                None => "NONE".to_string(),
            }
        )?;
        for frame in self.frames_get() {
            writeln!(f, "\t{:?}", frame)?;
        }
        write!(f, "")
    }
}

impl fmt::Display for Cosm {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?}", self)
    }
}

impl Cosm {
    /// Builds a Cosm from the *XB files. Path should _not_ contain file extension. Panics if the files could not be loaded.
    pub fn from_xb(filename: &str) -> Result<Self, NyxError> {
        Self::try_from_xb(Xb::from_file(filename)?)
    }

    /// Tries to load a subset of the DE438 XB from the embedded files, bounded between 01 Jan 2000 and 31 Dec 2050 TDB.
    pub fn try_de438() -> Result<Self, NyxError> {
        let de438_buf =
            EmbeddedAsset::get("de438s-00-50.xb").expect("Could not find asset de438s-00-550.xb");
        Self::try_from_xb(Xb::from_buffer(&de438_buf.data)?)
    }

    /// Load a subset of the DE438 XB from the embedded files, bounded between 01 Jan 2000 and 31 Dec 2050 TAI.
    pub fn de438() -> Arc<Self> {
        Arc::new(Self::try_de438().expect("could not load embedded de438s XB file"))
    }

    /// Load a subset of the DE438 XB from the embedded files, bounded between 01 Jan 2000 and 31 Dec 2050 TAI.
    pub fn de438_raw() -> Self {
        Self::try_de438().expect("could not load embedded de438s XB file")
    }

    /// Load a subset of the DE438 XB from the embedded files, bounded between 01 Jan 2000 and 31 Dec 2050 TAI.
    pub fn de438_gmat() -> Arc<Self> {
        let mut cosm = Self::try_de438().expect("could not load embedded de438s XB file");
        cosm.use_gmat_gm();
        Arc::new(cosm)
    }

    /// Attempts to build a Cosm from the XB files and the embedded IAU frames
    pub fn try_from_xb(xb: Xb) -> Result<Self, NyxError> {
        let mut cosm = Cosm {
            xb,
            frame_root: FrameTree {
                name: "SSB J2000".to_string(),
                frame: Frame::Celestial {
                    gm: SS_MASS * SUN_GM,
                    ephem_path: [None, None, None],
                    frame_path: [None, None, None],
                },
                parent_rotation: None,
                children: Vec::new(),
            },
            ephem2frame_map: HashMap::new(),
        };
        cosm.append_xb();
        cosm.load_iau_frames()?;
        Ok(cosm)
    }

    /// Switch the GM values to those from GMAT
    pub fn use_gmat_gm(&mut self) {
        // Set all of the GMs and their body fixed frames too
        self.frame_mut_gm("Sun J2000", 132_712_440_017.99);
        self.frame_mut_gm("IAU Sun", 132_712_440_017.99);
        self.frame_mut_gm("Mercury Barycenter J2000", 22_032.080_486_418);
        // No IAU mercury
        self.frame_mut_gm("Venus Barycenter J2000", 324_858.598_826_46);
        self.frame_mut_gm("IAU Venus", 324_858.598_826_46);
        self.frame_mut_gm("EME2000", 398_600.441_5);
        self.frame_mut_gm("IAU Earth", 398_600.441_5);
        self.frame_mut_gm("Luna", 4_902.800_582_147_8);
        self.frame_mut_gm("IAU Moon", 4_902.800_582_147_8);
        self.frame_mut_gm("Mars Barycenter J2000", 42_828.314258067);
        self.frame_mut_gm("IAU Mars", 42_828.314258067);
        self.frame_mut_gm("Jupiter Barycenter J2000", 126_712_767.857_80);
        self.frame_mut_gm("IAU Jupiter", 126_712_767.857_80);
        self.frame_mut_gm("Saturn Barycenter J2000", 37_940_626.061_137);
        self.frame_mut_gm("IAU Saturn", 37_940_626.061_137);
        self.frame_mut_gm("Uranus Barycenter J2000", 5_794_549.007_071_9);
        self.frame_mut_gm("IAU Uranus", 5_794_549.007_071_9);
        self.frame_mut_gm("Neptune Barycenter J2000", 6_836_534.063_879_3);
        self.frame_mut_gm("IAU Neptune", 6_836_534.063_879_3);
    }

    /// Load the IAU Frames as defined in Celest Mech Dyn Astr (2018) 130:22 (https://doi.org/10.1007/s10569-017-9805-5)
    pub fn load_iau_frames(&mut self) -> Result<(), NyxError> {
        // Load the IAU frames from the embedded TOML
        let iau_toml_str =
            EmbeddedAsset::get("iau_frames.toml").expect("Could not find iau_frames.toml as asset");
        self.append_frames(
            std::str::from_utf8(&iau_toml_str.data)
                .expect("Could not deserialize iau_frames.toml as string"),
        )
    }

    /// Returns the machine path of the ephemeris whose orientation is requested
    pub fn frame_find_path_for_orientation(&self, name: &str) -> Result<Vec<usize>, NyxError> {
        if self.frame_root.name == name {
            // Return an empty vector (but OK because we're asking for the root)
            Ok(Vec::new())
        } else {
            let mut path = Vec::new();
            Self::frame_find_path(name, &mut path, &self.frame_root)
        }
    }

    /// Seek a frame from its orientation name
    fn frame_find_path(
        frame_name: &str,
        cur_path: &mut Vec<usize>,
        f: &FrameTree,
    ) -> Result<Vec<usize>, NyxError> {
        if f.name == frame_name {
            Ok(cur_path.to_vec())
        } else if f.children.is_empty() {
            Err(NyxError::ObjectNotFound(frame_name.to_string()))
        } else {
            for child in &f.children {
                let mut this_path = cur_path.clone();
                let child_attempt = Self::frame_find_path(frame_name, &mut this_path, child);
                if let Ok(found_path) = child_attempt {
                    return Ok(found_path);
                }
            }
            // Could not find name in iteration, fail
            Err(NyxError::ObjectNotFound(frame_name.to_string()))
        }
    }

    /// Returns the correct frame for this ephemeris
    fn default_frame_value(
        e: &Ephemeris,
        ephem_path: [Option<usize>; 3],
        pos: usize,
    ) -> Option<FrameTree> {
        match e.constants.get("GM") {
            Some(gm) => {
                // It's a geoid, and we assume everything else is there
                let flattening = match e.constants.get("Flattening") {
                    Some(param) => param.value,
                    None => {
                        if e.name == "Moon" {
                            0.0012
                        } else {
                            0.0
                        }
                    }
                };
                let equatorial_radius = match e.constants.get("Equatorial radius") {
                    Some(param) => param.value,
                    None => {
                        if e.name == "Moon" {
                            1738.1
                        } else {
                            0.0
                        }
                    }
                };
                let semi_major_radius = match e.constants.get("Equatorial radius") {
                    Some(param) => {
                        if e.name == "Earth Barycenter" {
                            6378.1370
                        } else {
                            param.value
                        }
                    }
                    None => equatorial_radius, // assume spherical if unspecified
                };

                // Let's now build the J2000 version of this body
                Some(FrameTree {
                    name: format!("{} J2000", e.name.clone()),
                    frame: Frame::Geoid {
                        gm: gm.value,
                        flattening,
                        equatorial_radius,
                        semi_major_radius,
                        ephem_path,
                        frame_path: [Some(pos), None, None],
                    },
                    parent_rotation: None,
                    children: Vec::new(),
                })
            }
            None => {
                if e.name.to_lowercase() == *"sun" {
                    // Build the Sun frame in J2000
                    Some(FrameTree {
                        name: "Sun J2000".to_string(),
                        frame: Frame::Geoid {
                            gm: SUN_GM,
                            flattening: 0.0,
                            // From https://iopscience.iop.org/article/10.1088/0004-637X/750/2/135
                            equatorial_radius: 696_342.0,
                            semi_major_radius: 696_342.0,
                            ephem_path,
                            frame_path: [Some(pos), None, None],
                        },
                        parent_rotation: None,
                        children: Vec::new(),
                    })
                } else {
                    warn!("no GM value for XB {}", e.name);
                    None
                }
            }
        }
    }

    pub fn append_xb(&mut self) {
        // Insert the links between the SSB ephem and the J2000 frame (data stored in self.frame_root!)
        self.ephem2frame_map.insert(Vec::new(), Vec::new());

        // Build the frames
        for i in 0..self.xb.ephemeris_root.as_ref().unwrap().children.len() {
            if let Ok(child) = self.xb.ephemeris_from_path(&[i]) {
                // Add base J2000 frame for all barycenters
                if let Some(frame) = Self::default_frame_value(
                    child,
                    [Some(i), None, None],
                    self.frame_root.children.len(),
                ) {
                    self.frame_root.children.push(frame);
                    let frame_path = vec![self.frame_root.children.len() - 1];
                    self.ephem2frame_map.insert(vec![i], frame_path);
                }

                // Try to go one level deeper
                for j in 0..child.children.len() {
                    if let Ok(next_child) = self.xb.ephemeris_from_path(&[i, j]) {
                        // Create the frame
                        if let Some(frame) = Self::default_frame_value(
                            next_child,
                            [Some(i), Some(j), None],
                            self.frame_root.children.len(),
                        ) {
                            // At this stage, they are all children of the J2000 frame
                            // Bug: This should eventually use the orientation of the XB or it'll fail if it isn't J2000 based
                            self.frame_root.children.push(frame);
                            let frame_path = vec![self.frame_root.children.len() - 1];
                            self.ephem2frame_map.insert(vec![i, j], frame_path);
                        }
                    }
                }
            }
        }
    }

    /// Append Cosm with the contents of this TOML (must _not_ be the filename)
    pub fn append_frames(&mut self, toml_content: &str) -> Result<(), NyxError> {
        let maybe_frames: Result<frame_serde::FramesSerde, _> = toml::from_str(toml_content);
        match maybe_frames {
            Ok(mut frames) => {
                for (ref name, ref mut definition) in frames.frames.drain() {
                    if self.try_frame(name.as_str()).is_ok() {
                        warn!("overwriting frame `{}`", name);
                    }
                    if let Some(src_frame_name) = &definition.inherit {
                        match self.try_frame(src_frame_name.as_str()) {
                            Ok(src_frame) => {
                                definition.update_from(&src_frame);
                            }
                            Err(_) => error!(
                                "frame `{}` is derived from unknown frame `{}`, skipping!",
                                name, src_frame_name
                            ),
                        }
                    }
                    let rot = &definition.rotation;
                    let right_asc: Expr = match rot.right_asc.parse() {
                        Ok(expr) => expr,
                        Err(e) => {
                            let msg = format!("[frame.{}] - could not parse right_asc `{}` - are there any special characters? {}",
                            &name, &rot.right_asc, e);
                            error!("{}", msg);
                            return Err(NyxError::LoadingError(msg));
                        }
                    };
                    let declin: Expr = match rot.declin.parse() {
                        Ok(expr) => expr,
                        Err(e) => {
                            let msg = format!("[frame.{}] - could not parse declin `{}` - are there any special characters? {}",
                            &name, &rot.declin, e);
                            error!("{}", msg);
                            return Err(NyxError::LoadingError(msg));
                        }
                    };
                    let w_expr: Expr = match rot.w.parse() {
                        Ok(expr) => expr,
                        Err(e) => {
                            let msg = format!("[frame.{}] - could not parse w `{}` - are there any special characters? {}",
                            &name, &rot.w, e);
                            error!("{}", msg);
                            return Err(NyxError::LoadingError(msg));
                        }
                    };

                    let frame_rot = Euler3AxisDt::from_ra_dec_w(
                        right_asc,
                        declin,
                        w_expr,
                        match &rot.context {
                            Some(ctx) => ctx.clone(),
                            None => HashMap::new(),
                        },
                        match &rot.angle_unit {
                            Some(val) => AngleUnit::from_str(val.as_str()).unwrap(),
                            None => AngleUnit::Degrees,
                        },
                    );

                    // Let's now create the Frame, we'll add the ephem path and frame path just after
                    let mut new_frame = definition.as_frame();
                    let frame_name = name.replace("_", " ").trim().to_string();

                    // Grab the inherited frame again so we know how to place it in the frame tree
                    if let Some(src_frame_name) = &definition.inherit {
                        debug!("Loaded frame {}", frame_name);
                        let src_frame = self.try_frame(src_frame_name.as_str())?;
                        let mut fpath = src_frame.frame_path();
                        // And find the correct children
                        let children = match fpath.len() {
                            2 => {
                                &mut self.frame_root.children[fpath[0]].children[fpath[1]].children
                            }
                            1 => &mut self.frame_root.children[fpath[0]].children,
                            0 => &mut self.frame_root.children,
                            _ => unimplemented!("Too many children for now"),
                        };
                        fpath.push(children.len());
                        // Set the frame path and ephem path for this new frame
                        match new_frame {
                            Frame::Celestial {
                                ref mut ephem_path,
                                ref mut frame_path,
                                ..
                            }
                            | Frame::Geoid {
                                ref mut ephem_path,
                                ref mut frame_path,
                                ..
                            } => {
                                match fpath.len() {
                                    3 => {
                                        *frame_path =
                                            [Some(fpath[0]), Some(fpath[1]), Some(fpath[2])]
                                    }
                                    2 => *frame_path = [Some(fpath[0]), Some(fpath[1]), None],
                                    1 => *frame_path = [Some(fpath[0]), None, None],
                                    _ => unimplemented!(),
                                };

                                let epath = src_frame.ephem_path();
                                match epath.len() {
                                    3 => {
                                        *ephem_path =
                                            [Some(epath[0]), Some(epath[1]), Some(epath[2])]
                                    }
                                    2 => *ephem_path = [Some(epath[0]), Some(epath[1]), None],
                                    1 => *ephem_path = [Some(epath[0]), None, None],
                                    _ => unimplemented!(),
                                };
                            }
                            _ => unimplemented!(),
                        }

                        // And create and insert
                        // Create the new FrameTree node, and insert it as a child of the current path
                        let fnode = FrameTree {
                            name: frame_name,
                            frame: new_frame,
                            parent_rotation: Some(Box::new(frame_rot)),
                            children: Vec::new(),
                        };

                        children.push(fnode);
                    } else {
                        warn!(
                            "Frame `{}` does not inherit from anyone, cannot organize tree",
                            frame_name
                        );
                    }
                }
                Ok(())
            }
            Err(e) => {
                error!("{}", e);
                Err(NyxError::LoadingError(format!("{}", e)))
            }
        }
    }

    /// Returns the expected frame name with its ephemeris name for querying
    fn fix_frame_name(name: &str) -> String {
        let name = name.to_lowercase().trim().replace("_", " ");
        // Handle the specific frames
        if name == "eme2000" {
            String::from("Earth J2000")
        } else if name == "luna" {
            String::from("Moon J2000")
        } else if name == "earth moon barycenter" {
            String::from("Earth Barycenter J2000")
        } else if name == "ssb" || name == "ssb j2000" {
            String::from("SSB J2000")
        } else {
            let splt: Vec<_> = name.split(' ').collect();
            if splt[0] == "iau" {
                // This is an IAU frame, so the orientation is specified first, and we don't capitalize the ephemeris name
                vec![splt[0].to_string(), splt[1..splt.len()].join(" ")].join(" ")
            } else {
                // Likely a default center and frame, so let's do some clever guessing and capitalize the words
                let frame_name = capitalize(&splt[splt.len() - 1].to_string());
                let ephem_name = splt[0..splt.len() - 1]
                    .iter()
                    .map(|word| capitalize(word))
                    .collect::<Vec<_>>()
                    .join(" ");

                vec![ephem_name, frame_name].join(" ")
            }
        }
    }

    /// Fetch the frame associated with this ephemeris name
    /// This is slow, so avoid using it.
    pub fn try_frame(&self, name: &str) -> Result<Frame, NyxError> {
        let name = Self::fix_frame_name(name);
        if self.frame_root.name == name {
            // Return an empty vector (but OK because we're asking for the root)
            Ok(self.frame_root.frame)
        } else {
            let mut path = Vec::new();
            Ok(self.frame_from_frame_path(&FrameTree::frame_seek_by_name(
                &name,
                &mut path,
                &self.frame_root,
            )?))
        }
    }

    /// Provided an ephemeris path and an optional frame name, returns the Frame of that ephemeris.
    /// For example, if [3, 1] is provided (Moon in J2000 in the DE file), return Moon J2000
    /// If no frame name is provided, then the storage frame is returned. Otherwise, the correct frame is returned.
    pub fn frame_from_ephem_path(&self, ephem_path: &[usize]) -> Frame {
        self.frame_from_frame_path(self.ephem2frame_map.get(&ephem_path.to_vec()).unwrap())
    }

    /// Provided a frame path returns the Frame.
    pub fn frame_from_frame_path(&self, frame_path: &[usize]) -> Frame {
        match frame_path.len() {
            2 => self.frame_root.children[frame_path[0]].children[frame_path[1]].frame,
            1 => self.frame_root.children[frame_path[0]].frame,
            0 => self.frame_root.frame,
            _ => unimplemented!("Not expecting three layers of attitude frames"),
        }
    }

    fn frame_names(mut names: &mut Vec<String>, f: &FrameTree) {
        names.push(f.name.clone());
        for child in &f.children {
            Self::frame_names(&mut names, child);
        }
    }

    pub fn frames_get_names(&self) -> Vec<String> {
        let mut names = Vec::new();
        Self::frame_names(&mut names, &self.frame_root);
        names
    }

    fn frames(mut frames: &mut Vec<Frame>, f: &FrameTree) {
        frames.push(f.frame);
        for child in &f.children {
            Self::frames(&mut frames, child);
        }
    }

    /// Returns all of the frames defined in this Cosm
    pub fn frames_get(&self) -> Vec<Frame> {
        let mut frames = Vec::new();
        Self::frames(&mut frames, &self.frame_root);
        frames
    }

    /// Returns the geoid from the loaded XB, if it is in there, else panics!
    pub fn frame(&self, name: &str) -> Frame {
        self.try_frame(name).unwrap()
    }

    /// Mutates the GM value for the provided geoid id. Panics if ID not found.
    pub fn frame_mut_gm(&mut self, name: &str, new_gm: f64) {
        // Grab the frame -- this may panic!
        let frame_path = self.frame(name).frame_path();

        match frame_path.len() {
            2 => self.frame_root.children[frame_path[0]].children[frame_path[1]]
                .frame
                .gm_mut(new_gm),
            1 => self.frame_root.children[frame_path[0]].frame.gm_mut(new_gm),
            0 => self.frame_root.frame.gm_mut(new_gm),
            _ => unimplemented!("Not expecting three layers of attitude frames"),
        }
    }

    /// Returns the celestial state as computed from a de4xx.{FXB,XB} file in the original frame
    #[allow(clippy::comparison_chain)]
    pub fn raw_celestial_state(&self, path: &[usize], epoch: Epoch) -> Result<Orbit, NyxError> {
        if path.is_empty() {
            // This is the solar system barycenter, so we just return a state of zeros
            return Ok(Orbit::cartesian(
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                0.0,
                epoch,
                self.frame_root.frame,
            ));
        }
        let ephem = self.xb.ephemeris_from_path(path)?;

        // Compute the position as per the algorithm from jplephem
        let interp = ephem
            .interpolator
            .as_ref()
            .ok_or_else(|| NyxError::NoInterpolationData(ephem.name.clone()))?;

        // the DE file epochs are all in ET mod julian
        let start_mod_julian_f64 = ephem.start_epoch.as_ref().unwrap().as_raw();
        let coefficient_count: usize = interp.position_degree as usize;
        if coefficient_count <= 2 {
            // Cf. https://gitlab.com/chrisrabotin/nyx/-/issues/131
            return Err(NyxError::InvalidInterpolationData(format!(
                "position_degree is less than 3 for XB {}",
                ephem.name
            )));
        }

        let exb_states = match interp
            .state_data
            .as_ref()
            .ok_or_else(|| NyxError::NoStateData(ephem.name.clone()))?
        {
            EqualStates(states) => states,
            VarwindowStates(_) => panic!("variable window not yet supported by Cosm"),
        };

        let interval_length: f64 = exb_states.window_duration;

        let epoch_jde = epoch.as_jde_tdb_days();
        let delta_jde = epoch_jde - start_mod_julian_f64;

        let index_f = (delta_jde / interval_length).floor();
        let mut offset = delta_jde - index_f * interval_length;
        let mut index = index_f as usize;
        if index == exb_states.position.len() {
            index -= 1;
            offset = interval_length;
        } else if index > exb_states.position.len() {
            return Err(NyxError::NoInterpolationData(format!(
                "No interpolation data for date {}",
                epoch
            )));
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
        let storage_geoid = self.frame_from_ephem_path(path);
        Ok(Orbit::cartesian(
            x,
            y,
            z,
            vx / SECONDS_PER_DAY,
            vy / SECONDS_PER_DAY,
            vz / SECONDS_PER_DAY,
            epoch,
            storage_geoid,
        ))
    }

    /// Attempts to return the state of the celestial object at the provided time
    ///
    /// The light time correction is based on SPICE's implementation: https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/cspice/spkezr_c.html .
    /// Aberration computation is a conversion of the stelab function in SPICE, available here
    /// https://github.com/ChristopherRabotin/cspice/blob/26c72936fb7ff6f366803a1419b7cc3c61e0b6e5/src/cspice/stelab.c#L255
    pub fn try_celestial_state(
        &self,
        target_ephem: &[usize],
        datetime: Epoch,
        frame: Frame,
        correction: LightTimeCalc,
    ) -> Result<Orbit, NyxError> {
        let target_frame = self.frame_from_ephem_path(target_ephem);
        match correction {
            LightTimeCalc::None => {
                let state = Orbit::cartesian(0.0, 0.0, 0.0, 0.0, 0.0, 0.0, datetime, target_frame);
                Ok(-self.try_frame_chg(&state, frame)?)
            }
            LightTimeCalc::LightTime | LightTimeCalc::Abberation => {
                // Get the geometric states as seen from SSB
                let ssb2k = self.frame_root.frame;

                let obs = self.try_celestial_state(
                    &frame.ephem_path(),
                    datetime,
                    ssb2k,
                    LightTimeCalc::None,
                )?;
                let mut tgt =
                    self.try_celestial_state(target_ephem, datetime, ssb2k, LightTimeCalc::None)?;
                // It will take less than three iterations to converge
                for _ in 0..3 {
                    // Compute the light time
                    let lt = (tgt - obs).rmag() / SPEED_OF_LIGHT_KMS;
                    // Compute the new target state
                    let lt_dt = datetime - lt * Unit::Second;
                    tgt =
                        self.try_celestial_state(target_ephem, lt_dt, ssb2k, LightTimeCalc::None)?;
                }
                // Compute the correct state
                let mut state = Orbit::cartesian(
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

                if correction == LightTimeCalc::Abberation {
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

    /// Returns the state of the celestial object (target ephem) as seen in the requested frame at the provided time
    pub fn celestial_state(
        &self,
        target_ephem: &[usize],
        datetime: Epoch,
        frame: Frame,
        correction: LightTimeCalc,
    ) -> Orbit {
        self.try_celestial_state(target_ephem, datetime, frame, correction)
            .unwrap()
    }

    /// Return the position DCM (3x3) to go from the `from` frame to the `to` frame
    pub fn try_position_dcm_from_to(
        &self,
        from: &Frame,
        to: &Frame,
        dt: Epoch,
    ) -> Result<Matrix3<f64>, NyxError> {
        // And now let's compute the rotation path
        let mut dcm = Matrix3::<f64>::identity();

        if to.frame_path() == from.frame_path() {
            // No need to go any further
            return Ok(dcm);
        }

        let state_frame_path = from.frame_path();
        let new_frame_path = to.frame_path();

        // Let's get the translation path between both both states.
        let f_common_path = self.find_common_root(&new_frame_path, &state_frame_path);

        let get_dcm = |path: &[usize]| -> &FrameTree {
            // This is absolutely terrible, and there must be a better way to do it, but it's late.
            match path.len() {
                1 => &self.frame_root.children[path[0]],
                2 => &self.frame_root.children[path[0]].children[path[1]],
                3 => &self.frame_root.children[path[0]].children[path[1]].children[path[2]],
                _ => unimplemented!(),
            }
        };

        // Walk forward from the destination state

        for i in (f_common_path.len()..new_frame_path.len()).rev() {
            if let Some(parent_rot) = &get_dcm(&new_frame_path[0..=i]).parent_rotation {
                if let Some(next_dcm) = parent_rot.dcm_to_parent(dt) {
                    dcm *= next_dcm;
                }
            }
        }
        // Walk backward from current state up to common node (we transpose all backward rotations)
        for i in (f_common_path.len()..state_frame_path.len()).rev() {
            if let Some(parent_rot) = &get_dcm(&state_frame_path[0..=i]).parent_rotation {
                if let Some(next_dcm) = parent_rot.dcm_to_parent(dt) {
                    dcm *= next_dcm.transpose();
                }
            }
        }

        Ok(dcm)
    }

    /// Return the position and velocity DCM (6x6) to go from the `from` frame to the `to` frame
    #[allow(clippy::identity_op)]
    pub fn try_dcm_from_to(
        &self,
        from: &Frame,
        to: &Frame,
        dt: Epoch,
    ) -> Result<Matrix6<f64>, NyxError> {
        let r_dcm = self.try_position_dcm_from_to(from, to, dt)?;
        // Compute the dRdt DCM with finite differencing
        let pre_r_dcm = self.try_position_dcm_from_to(from, to, dt - 1 * Unit::Second)?;
        let post_r_dcm = self.try_position_dcm_from_to(from, to, dt + 1 * Unit::Second)?;

        Ok(dcm_finite_differencing(pre_r_dcm, r_dcm, post_r_dcm))
    }

    /// Return the position and velocity DCM (two 3x3 matrices) to go from the `from` frame to the `to` frame
    #[allow(clippy::identity_op)]
    pub fn try_dcm_from_to_in_parts(
        &self,
        from: &Frame,
        to: &Frame,
        dt: Epoch,
    ) -> Result<(Matrix3<f64>, Matrix3<f64>), NyxError> {
        let r_dcm = self.try_position_dcm_from_to(from, to, dt)?;
        // Compute the dRdt DCM with finite differencing
        let pre_r_dcm = self.try_position_dcm_from_to(from, to, dt - 1 * Unit::Second)?;
        let post_r_dcm = self.try_position_dcm_from_to(from, to, dt + 1 * Unit::Second)?;

        let drdt = 0.5 * post_r_dcm - 0.5 * pre_r_dcm;

        Ok((r_dcm, drdt))
    }

    /// Attempts to only perform a translation without rotation between two frames.
    /// You really shouldn't be using this unless you know exactly what you're doing.
    /// Typically, you want to use `try_frame_chg`.
    /// **WARNING:** This will update the Frame of the Orbit to the requested one EVEN IF it doesn't rotate it.
    pub fn try_frame_translation(
        &self,
        state: &Orbit,
        new_frame: Frame,
    ) -> Result<Orbit, NyxError> {
        let new_ephem_path = new_frame.ephem_path();
        let state_ephem_path = state.frame.ephem_path();

        // This doesn't make sense, but somehow the following algorithm only works when converting spacecraft states
        let mut new_state = if state.rmag() > 0.0 {
            *state
        } else if state_ephem_path.is_empty() {
            // SSB, let's invert this
            -*state
        } else {
            *state
        };

        // If we only need a rotation, let's skip trying to find the translation
        if new_ephem_path != state_ephem_path {
            // Let's get the translation path between both both states.
            let e_common_path = self.find_common_root(&new_ephem_path, &state_ephem_path);

            // This doesn't make sense, but somehow the following algorithm only works when converting spacecraft states
            new_state = if state.rmag() > 0.0 {
                // Walk backward from current state up to common node
                for i in (e_common_path.len()..state_ephem_path.len()).rev() {
                    let next_state =
                        self.raw_celestial_state(&state_ephem_path[0..=i], state.dt)?;
                    new_state += next_state;
                }

                // Walk forward from the destination state
                for i in (e_common_path.len()..new_ephem_path.len()).rev() {
                    let next_state = self.raw_celestial_state(&new_ephem_path[0..=i], state.dt)?;
                    new_state -= next_state;
                }

                new_state
            } else {
                let mut negated_fwd = false;

                // Walk forward from the destination state
                for i in (e_common_path.len()..new_ephem_path.len()).rev() {
                    let next_state = self.raw_celestial_state(&new_ephem_path[0..=i], state.dt)?;
                    if new_ephem_path.len() < state_ephem_path.len() && i == e_common_path.len() {
                        // We just crossed the common point going forward, so let's add the opposite of this state
                        new_state -= next_state;
                        negated_fwd = true;
                    } else {
                        new_state += next_state;
                    }
                }
                // Walk backward from current state up to common node
                for i in (e_common_path.len()..state_ephem_path.len()).rev() {
                    let next_state =
                        self.raw_celestial_state(&state_ephem_path[0..=i], state.dt)?;
                    if !negated_fwd && i == e_common_path.len() {
                        // We just crossed the common point (and haven't passed it going forward), so let's negate this state
                        new_state -= next_state;
                    } else {
                        new_state += next_state;
                    }
                }

                if negated_fwd {
                    // Because we negated the state going forward, let's flip it back to its correct orientation now.
                    -new_state
                } else {
                    new_state
                }
            };
        }
        new_state.frame = new_frame;
        Ok(new_state)
    }

    /// Attempts to return the provided state in the provided frame.
    pub fn try_frame_chg(&self, state: &Orbit, new_frame: Frame) -> Result<Orbit, NyxError> {
        if state.frame == new_frame {
            return Ok(*state);
        }
        // Let's perform the translation
        let mut new_state = self.try_frame_translation(state, new_frame)?;
        // And now let's compute the rotation path
        new_state.rotate_by(self.try_dcm_from_to(&state.frame, &new_frame, state.dt)?);
        Ok(new_state)
    }

    /// Return the provided state in the provided frame, or panics
    pub fn frame_chg(&self, state: &Orbit, new_frame: Frame) -> Orbit {
        self.try_frame_chg(state, new_frame).unwrap()
    }

    /// Returns the conversion path from the target ephemeris or frame `from` as seen from `to`.
    fn find_common_root(&self, from: &[usize], to: &[usize]) -> std::vec::Vec<usize> {
        let mut common_root = Vec::with_capacity(3); // Unlikely to be more than 3 items
        if !from.is_empty() && !to.is_empty() {
            // This is NOT the root of the ephemeris
            if from.len() < to.len() {
                // Iterate through the items in from
                for (n, obj) in from.iter().enumerate() {
                    if &to[n] == obj {
                        common_root.push(*obj);
                    } else {
                        // Found the end of the matching objects
                        break;
                    }
                }
            } else {
                // Iterate through the items in from
                for (n, obj) in to.iter().enumerate() {
                    if &from[n] == obj {
                        common_root.push(*obj);
                    } else {
                        // Found the end of the matching objects
                        break;
                    }
                }
            }
        }
        common_root
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cosmic::Bodies;

    #[test]
    fn test_cosm_indirect() {
        use crate::utils::is_diagonal;
        use std::f64::EPSILON;

        let jde = Epoch::from_gregorian_utc_at_midnight(2002, 2, 7);

        let cosm = Cosm::de438();

        let ven2ear = cosm.find_common_root(Bodies::Venus.ephem_path(), Bodies::Luna.ephem_path());
        assert_eq!(
            ven2ear.len(),
            0,
            "Venus -> (SSB) -> Earth Barycenter -> Earth Moon, therefore common root is zero lengthed"
        );

        // In this first part of the tests, we want to check that the DCM corresponds to no rotation, or none.
        // If the DCM is diagonal, then it does not have a rotation.

        assert!(
            is_diagonal(
                &cosm
                    .try_position_dcm_from_to(
                        &cosm.frame("Earth Barycenter J2000"),
                        &cosm.frame("Earth Barycenter J2000"),
                        jde
                    )
                    .unwrap()
            ),
            "Conversion does not require any rotation"
        );

        assert!(
            is_diagonal(
                &cosm
                    .try_position_dcm_from_to(
                        &cosm.frame("Venus Barycenter J2000"),
                        &cosm.frame("EME2000"),
                        jde
                    )
                    .unwrap()
            ),
            "Conversion does not require any rotation"
        );

        assert!(
            !is_diagonal(
                &cosm
                    .try_position_dcm_from_to(
                        &cosm.frame("Venus Barycenter J2000"),
                        &cosm.frame("IAU Sun"),
                        jde
                    )
                    .unwrap()
            ),
            "Conversion to Sun IAU from Venus J2k requires one rotation"
        );

        assert!(
            !is_diagonal(
                &cosm
                    .try_position_dcm_from_to(
                        &cosm.frame("IAU Sun"),
                        &cosm.frame("Venus Barycenter J2000"),
                        jde
                    )
                    .unwrap()
            ),
            "Conversion to Sun IAU from Venus J2k requires one rotation"
        );

        let c = LightTimeCalc::None;

        /*
        # Preceed all of the following python examples with
        >>> import spiceypy as sp
        >>> sp.furnsh('bsp/de438s.bsp')
        >>> et = 66312064.18493939
        */

        let earth_moon = cosm.frame("Luna");
        let ven2ear_state =
            cosm.celestial_state(Bodies::VenusBarycenter.ephem_path(), jde, earth_moon, c);
        assert_eq!(ven2ear_state.frame.ephem_path(), Bodies::Luna.ephem_path());
        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(1, et, "J2000", "NONE", 301)[0]]
        ['2.0512621957200775e+08', '-1.3561254792308527e+08', '-6.5578399676151529e+07', '3.6051374278177832e+01', '4.8889024622170766e+01', '2.0702933800843084e+01']
        */
        // NOTE: Venus position is quite off, not sure why.
        assert!(dbg!(ven2ear_state.x - 2.051_262_195_720_077_5e8).abs() < 5e-4);
        assert!(dbg!(ven2ear_state.y - -1.356_125_479_230_852_7e8).abs() < 7e-4);
        assert!(dbg!(ven2ear_state.z - -6.557_839_967_615_153e7).abs() < 4e-4);
        assert!(dbg!(ven2ear_state.vx - 3.605_137_427_817_783e1).abs() < 1e-8);
        assert!(dbg!(ven2ear_state.vy - 4.888_902_462_217_076_6e1).abs() < 1e-8);
        assert!(dbg!(ven2ear_state.vz - 2.070_293_380_084_308_4e1).abs() < 1e-8);

        // Check that conversion via a center frame works
        let earth_bary = cosm.frame("Earth Barycenter J2000");
        let moon_from_emb = cosm.celestial_state(Bodies::Luna.ephem_path(), jde, earth_bary, c);
        // Check this state again, as in the direct test
        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(301, et, "J2000", "NONE", 3)[0]]
        ['-8.1576591043050896e+04', '-3.4547568914480874e+05', '-1.4439185901465410e+05', '9.6071184439702662e-01', '-2.0358322542180365e-01', '-1.8380551745739407e-01']
        */
        assert_eq!(moon_from_emb.frame, earth_bary);
        assert!(dbg!(moon_from_emb.x - -8.157_659_104_305_09e4).abs() < 1e-4);
        assert!(dbg!(moon_from_emb.y - -3.454_756_891_448_087_4e5).abs() < 2e-5);
        assert!(dbg!(moon_from_emb.z - -1.443_918_590_146_541e5).abs() < 1e-5);
        assert!(dbg!(moon_from_emb.vx - 9.607_118_443_970_266e-1).abs() < 1e-8);
        assert!(dbg!(moon_from_emb.vy - -2.035_832_254_218_036_5e-1).abs() < 1e-8);
        assert!(dbg!(moon_from_emb.vz - -1.838_055_174_573_940_7e-1).abs() < 1e-8);

        let earth_from_emb = cosm.celestial_state(Bodies::Earth.ephem_path(), jde, earth_bary, c);
        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(399, et, "J2000", "NONE", 3)[0]]
        ['1.0033950894874154e+03', '4.2493637646888546e+03', '1.7760252107225667e+03', '-1.1816791248014408e-02', '2.5040812085717632e-03', '2.2608146685133296e-03']
        */
        assert!(dbg!(earth_from_emb.x - 1.003_395_089_487_415_4e3).abs() < 1e-6);
        assert!(dbg!(earth_from_emb.y - 4.249_363_764_688_855e3).abs() < 1e-6);
        assert!(dbg!(earth_from_emb.z - 1.776_025_210_722_566_7e3).abs() < 1e-6);
        assert!(dbg!(earth_from_emb.vx - -1.181_679_124_801_440_8e-2).abs() < 1e-9);
        assert!(dbg!(earth_from_emb.vy - 2.504_081_208_571_763e-3).abs() < 1e-9);
        assert!(dbg!(earth_from_emb.vz - 2.260_814_668_513_329_6e-3).abs() < 1e-9);

        let eme2k = cosm.frame("EME2000");
        let moon_from_earth = cosm.celestial_state(Bodies::Luna.ephem_path(), jde, eme2k, c);
        let earth_from_moon = cosm.celestial_state(Bodies::Earth.ephem_path(), jde, earth_moon, c);

        assert_eq!(earth_from_moon.radius(), -moon_from_earth.radius());
        assert_eq!(earth_from_moon.velocity(), -moon_from_earth.velocity());

        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(301, et, "J2000", "NONE", 399)[0]]
        ['-8.2579986132538310e+04', '-3.4972505290949758e+05', '-1.4616788422537665e+05', '9.7252863564504100e-01', '-2.0608730663037542e-01', '-1.8606633212590740e-01']
        */
        assert!(dbg!(moon_from_earth.x - -8.257_998_613_253_831e4).abs() < 1e-4);
        assert!(dbg!(moon_from_earth.y - -3.497_250_529_094_976e5).abs() < 1e-4);
        assert!(dbg!(moon_from_earth.z - -1.461_678_842_253_766_5e5).abs() < 1e-4);
        assert!(dbg!(moon_from_earth.vx - 9.725_286_356_450_41e-1).abs() < 1e-9);
        assert!(dbg!(moon_from_earth.vy - -2.060_873_066_303_754_2e-1).abs() < 1e-9);
        assert!(dbg!(moon_from_earth.vz - -1.860_663_321_259_074e-1).abs() < 1e-9);

        /*
        >>> ['{:.16e}'.format(x) for x in sp.spkez(10, et, "J2000", "NONE", 399)[0]]
        ['1.0965506591533598e+08', '-9.0570891031525031e+07', '-3.9266577019474506e+07', '2.0426570124555724e+01', '2.0412112498804031e+01', '8.8484257849460111e+00']
        */
        let sun2ear_state = cosm.celestial_state(Bodies::Sun.ephem_path(), jde, eme2k, c);
        let ssb_frame = cosm.frame("SSB");
        let emb_from_ssb =
            cosm.celestial_state(Bodies::EarthBarycenter.ephem_path(), jde, ssb_frame, c);
        let sun_from_ssb = cosm.celestial_state(Bodies::Sun.ephem_path(), jde, ssb_frame, c);
        let delta_state = sun2ear_state + (-sun_from_ssb + emb_from_ssb + earth_from_emb);

        assert!(delta_state.radius().norm() < EPSILON);
        assert!(delta_state.velocity().norm() < EPSILON);

        assert!(dbg!(sun2ear_state.x - 1.096_550_659_153_359_8e8).abs() < 1e-3);
        assert!(dbg!(sun2ear_state.y - -9.057_089_103_152_503e7).abs() < 1e-3);
        assert!(dbg!(sun2ear_state.z - -3.926_657_701_947_451e7).abs() < 1e-3);
        assert!(dbg!(sun2ear_state.vx - 2.042_657_012_455_572_4e1).abs() < 1e-9);
        assert!(dbg!(sun2ear_state.vy - 2.041_211_249_880_403e1).abs() < 1e-9);
        assert!(dbg!(sun2ear_state.vz - 8.848_425_784_946_011).abs() < 1e-9);
        // And check the converse
        let sun2k = cosm.frame("Sun J2000");
        let sun2ear_state = cosm.celestial_state(&sun2k.ephem_path(), jde, eme2k, c);
        let ear2sun_state = cosm.celestial_state(&eme2k.ephem_path(), jde, sun2k, c);
        let state_sum = ear2sun_state + sun2ear_state;
        assert!(state_sum.rmag() < 1e-8);
        assert!(state_sum.vmag() < 1e-11);
    }

    #[test]
    fn test_cosm_frame_change_earth2luna() {
        let cosm = Cosm::de438();
        let eme2k = cosm.frame("EME2000");
        let luna = cosm.frame("Luna");

        let jde = Epoch::from_jde_et(2_458_823.5);
        // From JPL HORIZONS
        let lro = Orbit::cartesian(
            4.017_685_334_718_784E5,
            2.642_441_356_763_487E4,
            -3.024_209_691_251_325E4,
            -6.168_920_999_978_097E-1,
            -6.678_258_076_726_339E-1,
            4.208_264_479_358_517E-1,
            jde,
            eme2k,
        );

        let lro_jpl = Orbit::cartesian(
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
        let lro_wrt_earth = cosm.frame_chg(&lro_wrt_moon, eme2k);
        assert!((lro_wrt_earth - lro).rmag() < std::f64::EPSILON);
        assert!((lro_wrt_earth - lro).vmag() < std::f64::EPSILON);
    }

    #[test]
    fn test_cosm_frame_change_ven2luna() {
        let cosm = Cosm::de438();
        let luna = cosm.frame("Luna");
        let venus = cosm.frame("Venus Barycenter J2000");

        let jde = Epoch::from_jde_et(2_458_823.5);
        // From JPL HORIZONS
        let lro = Orbit::cartesian(
            -4.393_308_217_174_602E7,
            1.874_075_194_166_327E8,
            8.763_986_396_329_135E7,
            -5.054_051_490_556_286E1,
            -1.874_720_232_671_061E1,
            -6.518_342_268_306_54,
            jde,
            venus,
        );

        let lro_jpl = Orbit::cartesian(
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
        assert!(lro_moon_earth_delta.rmag() < 0.3);
        assert!(lro_moon_earth_delta.vmag() < 1e-5);
        // And the converse
        let lro_wrt_venus = cosm.frame_chg(&lro_wrt_moon, venus);
        assert!((lro_wrt_venus - lro).rmag() < std::f64::EPSILON);
        assert!((lro_wrt_venus - lro).vmag() < std::f64::EPSILON);
    }

    #[test]
    fn test_cosm_frame_change_ssb2luna() {
        let cosm = Cosm::de438();
        let luna = cosm.frame("Luna");
        let ssb = cosm.frame("SSB");

        let jde = Epoch::from_jde_et(2_458_823.5);
        // From JPL HORIZONS
        let lro = Orbit::cartesian(
            4.227_396_973_787_854E7,
            1.305_852_533_250_192E8,
            5.657_002_470_685_254E7,
            -2.964_638_617_895_494E1,
            7.078_704_012_700_072,
            3.779_568_779_111_446,
            jde,
            ssb,
        );

        let lro_jpl = Orbit::cartesian(
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
        assert!(dbg!(lro_moon_earth_delta.rmag()) < 0.25);
        assert!(dbg!(lro_moon_earth_delta.vmag()) < 1e-5);
        // And the converse
        let lro_wrt_ssb = cosm.frame_chg(&lro_wrt_moon, ssb);
        assert!((lro_wrt_ssb - lro).rmag() < std::f64::EPSILON);
        assert!((lro_wrt_ssb - lro).vmag() < std::f64::EPSILON);
    }

    #[test]
    fn test_cosm_lt_corr() {
        let cosm = Cosm::de438();

        let jde = Epoch::from_jde_et(2_452_312.500_742_881);

        let mars2k = cosm.frame("Mars Barycenter J2000");

        let out_state = cosm.celestial_state(
            Bodies::EarthBarycenter.ephem_path(),
            jde,
            mars2k,
            LightTimeCalc::LightTime,
        );

        // Note that the following data comes from SPICE (via spiceypy).
        // There is currently a difference in computation for de438s: https://github.com/brandon-rhodes/python-jplephem/issues/33 .
        // However, in writing this test, I also checked the computed light time, which matches SPICE to 2.999058779096231e-10 seconds.
        assert!(dbg!(out_state.x - -2.577_185_470_734_315_8e8).abs() < 1e-3);
        assert!(dbg!(out_state.y - -5.814_057_247_686_307e7).abs() < 1e-2);
        assert!(dbg!(out_state.z - -2.493_960_187_215_911_6e7).abs() < 1e-3);
        assert!(dbg!(out_state.vx - -3.460_563_654_257_750_7).abs() < 1e-7);
        assert!(dbg!(out_state.vy - -3.698_207_386_702_523_5e1).abs() < 1e-7);
        assert!(dbg!(out_state.vz - -1.690_807_917_994_789_7e1).abs() < 1e-7);
    }

    #[test]
    fn test_cosm_aberration_corr() {
        let cosm = Cosm::de438();

        let jde = Epoch::from_jde_et(2_452_312.500_742_881);

        let mars2k = cosm.frame("Mars Barycenter J2000");

        let out_state = cosm.celestial_state(
            Bodies::EarthBarycenter.ephem_path(),
            jde,
            mars2k,
            LightTimeCalc::Abberation,
        );

        assert!(dbg!(out_state.x - -2.577_231_712_700_484_4e8).abs() < 1e-3);
        assert!(dbg!(out_state.y - -5.812_356_237_533_56e7).abs() < 1e-2);
        assert!(dbg!(out_state.z - -2.493_146_410_521_204_8e7).abs() < 1e-3);
        // Reenable this test after #96 is implemented.
        dbg!(out_state.vx - -3.463_585_965_206_417);
        dbg!(out_state.vy - -3.698_169_177_803_263e1);
        dbg!(out_state.vz - -1.690_783_648_756_073e1);
    }

    #[test]
    fn test_cosm_rotation_validation() {
        let jde = Epoch::from_gregorian_utc_at_midnight(2002, 2, 7);
        let cosm = Cosm::de438();

        println!("Available ephems: {:?}", cosm.xb.ephemeris_get_names());
        println!("Available frames: {:?}", cosm.frames_get_names());

        let sun2k = cosm.frame("Sun J2000");
        let sun_iau = cosm.frame("IAU Sun");
        let ear_sun_2k =
            cosm.celestial_state(Bodies::Earth.ephem_path(), jde, sun2k, LightTimeCalc::None);
        let ear_sun_iau = cosm.frame_chg(&ear_sun_2k, sun_iau);
        let ear_sun_2k_prime = cosm.frame_chg(&ear_sun_iau, sun2k);

        assert!(
            (ear_sun_2k.rmag() - ear_sun_iau.rmag()).abs() <= 2e-16,
            "a single rotation changes rmag"
        );
        assert!(
            (ear_sun_2k_prime - ear_sun_2k).rmag() <= 1e-7,
            "reverse rotation does not match initial state"
        );

        // Test an EME2k to Earth IAU rotation

        let eme2k = cosm.frame("EME2000");
        let earth_iau = cosm.frame("IAU Earth");
        println!("{:?}\n{:?}", eme2k, earth_iau);

        let dt = Epoch::from_gregorian_utc(2023, 11, 16, 06, 11, 19, 146200000);

        // Case 1: Initialize from EME2000
        let state_eme2k = Orbit::cartesian(
            2196.260617879428,
            5161.156645108604,
            3026.122639999121,
            -0.376262062797,
            0.159851964260,
            0.000983174100,
            dt,
            eme2k,
        );

        let state_iau_earth = Orbit::cartesian(
            925.651012109672,
            -5529.343261192083,
            3031.179663596684,
            0.000032603895,
            -0.000229009865,
            0.000108943879,
            dt,
            earth_iau,
        );

        let state_iau_earth_computed = cosm.frame_chg(&state_eme2k, earth_iau);
        let delta_state = cosm.frame_chg(&state_iau_earth_computed, eme2k) - state_eme2k;

        println!("{}", delta_state);

        assert!(
            delta_state.rmag().abs() < 1e-11,
            "Inverse rotation is broken"
        );
        assert!(
            delta_state.vmag().abs() < 1e-11,
            "Inverse rotation is broken"
        );

        // Check the state matches
        let dcm = cosm
            .try_position_dcm_from_to(&cosm.frame("EME2000"), &cosm.frame("IAU_Moon"), dt)
            .unwrap();
        // The DCM matches SPICE very well

        let spice_dcm = Matrix3::new(
            -0.0517208410530453,
            0.9259079323890239,
            0.3741917360656812,
            -0.9986038252635100,
            -0.0439204209304354,
            -0.0293495620815480,
            -0.0107403337867543,
            -0.3751872830525792,
            0.9268868042354325,
        );
        assert!(
            (dcm - spice_dcm).norm() < 1e-9,
            "DCM EME2000 -> IAU Moon error too large"
        );

        let (pos_rss, vel_rss) = state_iau_earth.rss(&state_iau_earth_computed);
        dbg!(pos_rss, vel_rss);

        assert!(
            state_iau_earth.eq_within(&state_iau_earth_computed, 1e-4, 1e-4),
            "EME2000 -> IAU Earth failed\nComputed: {}\nExpected: {}",
            state_iau_earth_computed,
            state_iau_earth,
        );

        // Case 2
        // Earth Body Fixed state:
        // State (km, km/sec)
        // 'Earth' -> 'test' in 'Earth Body Fixed' at '31-JAN-2000 12:00:00.0000 TAI'
        // Pos:  3.092802381110541e+02 -3.431791232988777e+03  6.891017545171710e+03
        // Vel:  6.917077556761001e+00  6.234631407415389e-01  4.062487128428244e-05
        let state_eme2k = Orbit::cartesian(
            -2436.45,
            -2436.45,
            6891.037,
            5.088_611,
            -5.088_611,
            0.0,
            Epoch::from_gregorian_tai_at_noon(2000, 1, 31),
            eme2k,
        );

        let state_ecef = cosm.frame_chg(&state_eme2k, earth_iau);
        println!("{}\n{}", state_eme2k, state_ecef);
        assert!((state_ecef.x - 309.280_238_111_054_1).abs() < 1e-5);
        assert!((state_ecef.y - -3_431.791_232_988_777).abs() < 1e-5);
        assert!((state_ecef.z - 6_891.017_545_171_71).abs() < 1e-5);

        // Case 3
        // Earth Body Fixed state:
        // State (km, km/sec)
        // 'Earth' -> 'test' in 'Earth Body Fixed' at '01-MAR-2000 12:00:00.0000 TAI'
        // Pos: -1.424497118292030e+03 -3.137502417055381e+03  6.890998090503171e+03
        // Vel:  6.323912379829687e+00 -2.871020900962905e+00  8.125749038014632e-05
        let state_eme2k = Orbit::cartesian(
            -2436.45,
            -2436.45,
            6891.037,
            5.088_611,
            -5.088_611,
            0.0,
            Epoch::from_gregorian_tai_at_noon(2000, 3, 1),
            eme2k,
        );

        let state_ecef = cosm.frame_chg(&state_eme2k, earth_iau);
        println!("{}\n{}", state_eme2k, state_ecef);
        assert!((state_ecef.x - -1_424.497_118_292_03).abs() < 1e-5);
        assert!((state_ecef.y - -3_137.502_417_055_381).abs() < 1e-5);
        assert!((state_ecef.z - 6_890.998_090_503_171).abs() < 1e-5);

        // Ground station example
        let dt = Epoch::from_gregorian_tai_hms(2020, 1, 1, 0, 0, 20);
        let gs = Orbit::cartesian(
            -4461.153491497329,
            2682.445251105359,
            -3674.3793821716713,
            -0.19560699645796042,
            -0.32531244947129817,
            0.0,
            dt,
            earth_iau,
        );
        let gs_eme = cosm.frame_chg(&gs, eme2k);
        println!("{}\n{}", gs, gs_eme);
    }

    #[test]
    fn test_cosm_fix_frame_name() {
        assert_eq!(
            Cosm::fix_frame_name("Mars barycenter j2000"),
            "Mars Barycenter J2000"
        );
    }

    #[test]
    fn test_cosm_rotation_spiceypy_pos_dcm() {
        // These validation tests are from tests/spiceypy/rotations.py
        use crate::linalg::{Matrix3, Vector3};
        use std::f64::EPSILON;
        let cosm = Cosm::de438();

        let et0 = Epoch::from_gregorian_utc_at_noon(2022, 11, 30);
        // Moon Earth J2000
        let out_state = cosm.celestial_state(
            Bodies::Luna.ephem_path(),
            et0,
            cosm.frame("EME2000"),
            LightTimeCalc::None,
        );
        println!("{}", out_state);
        let sp_ex = Vector3::new(341456.50984349, -123580.72487638, -87633.23833676);
        println!(
            "{}\n{:.3} m",
            out_state.radius() - sp_ex,
            (out_state.radius() - sp_ex).norm() * 1e3,
        );
        // Check that we are below 1 meter of error
        assert!(
            (out_state.radius() - sp_ex).norm() * 1e3 < 1.0,
            "Larger than 1 meter error"
        );

        // Moon IAU Earth
        let out_state = cosm.celestial_state(
            Bodies::Luna.ephem_path(),
            et0,
            cosm.frame("IAU Earth"),
            LightTimeCalc::None,
        );
        println!("{}", out_state);
        let sp_ex = Vector3::new(-7239.63398824, 363242.64460769, -86871.72713248);
        println!(
            "{}\n{:.3} m",
            out_state.radius() - sp_ex,
            (out_state.radius() - sp_ex).norm() * 1e3,
        );

        assert!(
            (out_state.radius() - sp_ex).norm() * 1e3 < 1.0,
            "Larger than 1 meter error"
        );

        // IAU Earth <-> EME2000
        println!("\n=== IAU Earth <-> EME2000 ===\n");
        let dcm = cosm
            .try_position_dcm_from_to(&cosm.frame("iau_earth"), &cosm.frame("EME2000"), et0)
            .unwrap();

        // Position rotation first
        let sp_ex = Matrix3::new(
            // First row
            -3.58819172e-01,
            9.33404435e-01,
            2.22748179e-03,
            // Second row
            -9.33406755e-01,
            -3.58820051e-01,
            -5.70997090e-06,
            // Third row
            7.93935415e-04,
            -2.08119539e-03,
            9.99997519e-01,
        );

        println!("{}", (dcm - sp_ex).norm());
        assert!((dcm - sp_ex).norm() < 1e-6, "3x3 DCM error");

        let dcm_return = cosm
            .try_position_dcm_from_to(&cosm.frame("EME2000"), &cosm.frame("iau_earth"), et0)
            .unwrap();
        println!("{}", (dcm.transpose() - dcm_return).norm()); // TODO: Add as test!

        // IAU Earth <-> IAU Mars
        println!("\n=== IAU Earth <-> IAU Mars ===\n");
        let dcm = cosm
            .try_position_dcm_from_to(&cosm.frame("iau_earth"), &cosm.frame("iau mars"), et0)
            .unwrap();

        // Position rotation first
        let sp_ex = Matrix3::new(
            2.58117084e-01,
            7.55715666e-01,
            -6.01888198e-01,
            -9.40722808e-01,
            3.38489522e-01,
            2.15741174e-02,
            2.20036747e-01,
            5.60641307e-01,
            7.98288891e-01,
        );

        println!("{}", (dcm - sp_ex).norm());
        assert!((dcm - sp_ex).norm() < 1e-1, "3x3 DCM error"); // Error larger here because IAU Mars defined differently in SPICE than IAU

        let dcm_return = cosm
            .try_position_dcm_from_to(&cosm.frame("iau mars"), &cosm.frame("iau_earth"), et0)
            .unwrap();
        assert!(
            (dcm.transpose() - dcm_return).norm() < EPSILON,
            "Return DCM is not the transpose of the forward"
        );
    }

    #[test]
    fn test_cosm_rotation_spiceypy_dcm() {
        // These validation tests are from tests/spiceypy/rotations.py
        use crate::linalg::Matrix6;
        let cosm = Cosm::de438();

        let et0 = Epoch::from_gregorian_utc_at_noon(2022, 11, 30);

        // IAU Earth <-> EME2000
        println!("\n=== IAU Earth <-> EME2000 ===\n");
        let dcm = cosm
            .try_dcm_from_to(&cosm.frame("iau_earth"), &cosm.frame("EME2000"), et0)
            .unwrap();

        let sp_iau_earth_2_eme2k = Matrix6::new(
            // First row
            -3.58819172e-01,
            9.33404435e-01,
            2.22748179e-03,
            0.0,
            0.0,
            0.0,
            // Second row
            -9.33406755e-01,
            -3.58820051e-01,
            -5.70997090e-06,
            0.0,
            0.0,
            0.0,
            // Third row
            7.93935415e-04,
            -2.08119539e-03,
            9.99997519e-01,
            0.0,
            0.0,
            0.0,
            // Fourth row
            6.80649250e-05,
            2.61655068e-05,
            3.08051436e-12,
            -3.58819172e-01,
            9.33404435e-01,
            2.22748179e-03,
            // Fifth row
            -2.61655708e-05,
            6.80650942e-05,
            -1.57934025e-14,
            -9.33406755e-01,
            -3.58820051e-01,
            -5.70997090e-06,
            // Sixth row
            -1.51762071e-07,
            -5.78975647e-08,
            -6.86189683e-15,
            7.93935415e-04,
            -2.08119539e-03,
            9.99997519e-01,
        );

        assert!(
            (dcm - sp_iau_earth_2_eme2k).norm() < 1e-8,
            "IAU Earth -> EME2000"
        );

        let sp_eme2k_2_iau_earth = Matrix6::new(
            // First row
            -3.58819172e-01,
            -9.33406755e-01,
            7.93935415e-04,
            0.0,
            0.0,
            0.0,
            // Second row
            9.33404435e-01,
            -3.58820051e-01,
            -2.08119539e-03,
            0.0,
            0.0,
            0.0,
            // Third row
            2.22748179e-03,
            -5.70997090e-06,
            9.99997519e-01,
            0.0,
            0.0,
            0.0,
            // Fourth row
            6.80649250e-05,
            -2.61655708e-05,
            -1.51762071e-07,
            -3.58819172e-01,
            -9.33406755e-01,
            7.93935415e-04,
            // Fifth row
            2.61655068e-05,
            6.80650942e-05,
            -5.78975647e-08,
            9.33404435e-01,
            -3.58820051e-01,
            -2.08119539e-03,
            // Sixth row
            3.08051436e-12,
            -1.57934025e-14,
            -6.86189683e-15,
            2.22748179e-03,
            -5.70997090e-06,
            9.99997519e-01,
        );

        let dcm_return = cosm
            .try_dcm_from_to(&cosm.frame("EME2000"), &cosm.frame("iau_earth"), et0)
            .unwrap();

        assert!(
            (dcm_return - sp_eme2k_2_iau_earth).norm() < 1e-8,
            "EME2000 -> IAU Earth"
        );

        // IAU Earth <-> IAU Mars
        println!("\n=== IAU Earth <-> IAU Mars ===\n");
        let dcm = cosm
            .try_dcm_from_to(&cosm.frame("iau mars"), &cosm.frame("iau earth"), et0)
            .unwrap();

        // Position rotation first
        let sp_iau_mars_2_iau_earth = Matrix6::new(
            // First row
            2.58117084e-01,
            -9.40722808e-01,
            2.20036747e-01,
            0.0,
            0.0,
            0.0,
            // Second row
            7.55715666e-01,
            3.38489522e-01,
            5.60641307e-01,
            0.0,
            0.0,
            0.0,
            // Third row
            -6.01888198e-01,
            2.15741174e-02,
            7.98288891e-01,
            0.0,
            0.0,
            0.0,
            // Fourth row
            -1.15728287e-05,
            6.38714373e-06,
            4.08826103e-05,
            2.58117084e-01,
            -9.40722808e-01,
            2.20036747e-01,
            // Fifth row
            5.17068221e-06,
            1.50318153e-05,
            -1.60453349e-05,
            7.55715666e-01,
            3.38489522e-01,
            5.60641307e-01,
            // Sixth row
            1.52922212e-06,
            4.26631500e-05,
            1.17187529e-12,
            -6.01888198e-01,
            2.15741174e-02,
            7.98288891e-01,
        );

        assert!(
            (dcm - sp_iau_mars_2_iau_earth).norm() < 1e-1,
            "IAU Mars -> IAU Earth failed"
        );

        // And backward
        // Position rotation first
        let sp_iau_earth_2_iau_mars = Matrix6::new(
            // First row
            2.58117084e-01,
            7.55715666e-01,
            -6.01888198e-01,
            0.0,
            0.0,
            0.0,
            // Second row
            -9.40722808e-01,
            3.38489522e-01,
            2.15741174e-02,
            0.0,
            0.0,
            0.0,
            // Third row
            2.20036747e-01,
            5.60641307e-01,
            7.98288891e-01,
            0.0,
            0.0,
            0.0,
            // Fourth row
            -1.15728287e-05,
            5.17068221e-06,
            1.52922212e-06,
            2.58117084e-01,
            7.55715666e-01,
            -6.01888198e-01,
            // Fifth row
            6.38714373e-06,
            1.50318153e-05,
            4.26631500e-05,
            -9.40722808e-01,
            3.38489522e-01,
            2.15741174e-02,
            // Sixth row
            4.08826103e-05,
            -1.60453349e-05,
            1.17187529e-12,
            2.20036747e-01,
            5.60641307e-01,
            7.98288891e-01,
        );

        let dcm_return = cosm
            .try_dcm_from_to(&cosm.frame("iau earth"), &cosm.frame("iau mars"), et0)
            .unwrap();

        // Allow for "large" error because of different IAU Mars definition
        assert!(
            (dcm_return - sp_iau_earth_2_iau_mars).norm() < 1e-1,
            "IAU Earth -> IAU Mars failed"
        );
    }

    #[test]
    fn cosm_spice_concrete_validation() {
        // Example from sc_query.py from an actual mission design
        let cosm = Cosm::de438();
        let mj2k = cosm.frame("luna");
        let eme2k = cosm.frame("eme2000");

        let et: Epoch = Epoch::from_gregorian_utc(2022, 11, 29, 6, 47, 5, 0);
        println!("{:.6}", et.as_tdb_seconds());

        let moon_state = Orbit::cartesian(
            -274450.055,
            197469.145,
            120487.744,
            7.18981408,
            -2.46517900,
            -1.98864578,
            et,
            mj2k,
        );

        let moon_into_eme2k = cosm.frame_chg(&moon_state, eme2k);
        let sp_state_in_eme2k = Orbit::cartesian(
            4.50272491e+03,
            -8.84403138e+03,
            -5.43867081e+03,
            7.91088598e+00,
            -1.75366723e+00,
            -1.67218499e+00,
            et,
            mj2k,
        );

        println!(
            "err: {} \t {}\nWas: {}\nIs:  {}",
            (moon_into_eme2k - sp_state_in_eme2k).rmag(),
            (moon_into_eme2k - sp_state_in_eme2k).vmag(),
            moon_state,
            moon_into_eme2k
        );
        // SPICE only returns data down to the meter level, so we can't check less than that
        assert!(
            (moon_into_eme2k - sp_state_in_eme2k).rmag() < 1e-3,
            "error greater than expected"
        );
        // SPICE only returns data down to the millimeter per second level, so we can't check less than that
        assert!(
            (moon_into_eme2k - sp_state_in_eme2k).vmag() < 1e-6,
            "error greater than expected"
        );

        // End of transfer
        let et: Epoch =
            Epoch::from_gregorian_utc(2022, 12, 4, 11, 59, 51, 0) + 884000 * Unit::Microsecond;
        println!("{:.6}", et.as_tdb_seconds());

        let moon_state = Orbit::cartesian(
            482.668990,
            -1202.72224,
            1766.54808,
            0.575668933,
            -0.975029706,
            -1.95884493,
            et,
            mj2k,
        );

        let moon_into_eme2k = cosm.frame_chg(&moon_state, eme2k);

        let sp_state_in_eme2k = Orbit::cartesian(
            3.37603919e+05,
            1.79697913e+05,
            7.14382236e+04,
            1.11062935e-01,
            -1.89556775e-01,
            -1.52151242e+00,
            et,
            mj2k,
        );

        println!(
            "err: {} \t {}\nWas: {}\nIs:  {}",
            (moon_into_eme2k - sp_state_in_eme2k).rmag(),
            (moon_into_eme2k - sp_state_in_eme2k).vmag(),
            moon_state,
            moon_into_eme2k
        );

        // SPICE only returns data down to the meter level, so we can't check less than that
        assert!(
            (moon_into_eme2k - sp_state_in_eme2k).rmag() < 1e-3,
            "error greater than expected"
        );
        // SPICE only returns data down to the millimeter per second level, so we can't check less than that
        assert!(
            (moon_into_eme2k - sp_state_in_eme2k).vmag() < 1e-6,
            "error greater than expected"
        );
    }

    #[test]
    fn debug_cosm() {
        dbg!(Cosm::de438_gmat());
    }

    #[test]
    fn why_broken() {
        let e = Epoch::from_gregorian_tai_hms(2002, 02, 14, 0, 0, 0);
        println!("{}", e.as_tdb_seconds());
        println!("{}", e.as_jde_tdb_duration().in_unit(Unit::Second));
    }
}
