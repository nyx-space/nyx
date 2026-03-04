/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

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
use crate::md::trajectory::INTERPOLATION_SAMPLES;
use crate::md::StateParameter;
use crate::od::DynamicsError;
use crate::{cosmic::State, md::prelude::Interpolatable};
use anise::analysis::prelude::OrbitalElement;
use anise::errors::PhysicsError;
use anise::math::interpolation::{hermite_eval, InterpolationError};
use anise::{
    astro::Location,
    prelude::{Frame, Orbit},
};
use core::error::Error;
use core::fmt;
use core::ops::Add;
use hifitime::Epoch;
use nalgebra::{Const, DimName, OMatrix, OVector, Vector3};
pub mod ground_dynamics;
pub mod sensitivity;
pub mod trk_device;

/// Represents a ground position/nav/timing receiver, e.g. a customer
/// Note that we rebuild the Location structure from ANISE but _without_ a terrain mask because the mask is not copyable
/// and this PNTRx must be copyable to implement State.
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct GroundAsset {
    pub latitude_deg: f64,
    pub longitude_deg: f64,
    pub height_km: f64,
    // Velocity in the SEZ frame, South component, in METERS per second
    pub latitude_rate_deg_s: f64,
    // Velocity in the SEZ frame, East component, in METERS per second
    pub longitude_rate_deg_s: f64,
    // Velocity in the SEZ frame, Up/+Z component, in METERS per second
    pub height_rate_km_s: f64,
    // Epoch
    pub epoch: Epoch,
    /// Frame on which this location rests
    pub frame: Frame,
    pub stm: Option<OMatrix<f64, Const<6>, Const<6>>>,
}

impl GroundAsset {
    pub fn from_fixed(
        latitude_deg: f64,
        longitude_deg: f64,
        height_km: f64,
        epoch: Epoch,
        frame: Frame,
    ) -> Self {
        Self {
            latitude_deg,
            longitude_deg,
            height_km,
            epoch,
            frame,
            ..Default::default()
        }
    }

    pub fn with_velocity_sez_m_s(
        mut self,
        vel_s_m_s: f64,
        vel_e_m_s: f64,
        vel_z_m_s: f64,
    ) -> Result<Self, Box<dyn Error>> {
        let (lat_rate_deg_s, long_rate_deg_s, alt_rate_km_s) = latlongalt_rate(
            self.orbit(),
            Vector3::new(vel_s_m_s, vel_e_m_s, vel_z_m_s) * 1e-3,
        )?;

        self.latitude_rate_deg_s = lat_rate_deg_s;
        self.longitude_rate_deg_s = long_rate_deg_s;
        self.height_rate_km_s = alt_rate_km_s;

        Ok(self)
    }

    pub fn to_location(&self) -> Location {
        Location {
            latitude_deg: self.latitude_deg,
            longitude_deg: self.longitude_deg,
            height_km: self.height_km,
            frame: self.frame.into(),
            ..Default::default()
        }
    }

    /// Compute the velocity in m/s in the SEZ frame from the state data stored in latitude deg/s, longitude deg/s, and height in km/s for integration of EOMs
    pub fn velocity_sez_m_s(&self) -> Result<OVector<f64, Const<3>>, PhysicsError> {
        // First, convert from SEZ to body frame.
        let rx = Orbit::try_latlongalt(
            self.latitude_deg,
            self.longitude_deg,
            self.height_km,
            self.epoch,
            self.frame,
        )?;

        velocity_sez_from_latlongalt_rate(
            rx,
            self.latitude_rate_deg_s,
            self.longitude_rate_deg_s,
            self.height_rate_km_s,
        )
        .map(|v| v * 1e3)
    }
}

impl Default for GroundAsset {
    fn default() -> Self {
        Self {
            frame: Frame::from_ephem_j2000(399),
            latitude_deg: 0.,
            longitude_deg: 0.,
            height_km: 0.,
            latitude_rate_deg_s: 0.,
            longitude_rate_deg_s: 0.,
            height_rate_km_s: 0.0,
            epoch: Epoch::from_tdb_seconds(0.0),
            stm: None,
        }
    }
}

impl fmt::Display for GroundAsset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "lat.: {:.3} deg\tlong.: {:.3} deg\talt.: {:.3} km (speed: {:.3} m/s)",
            self.latitude_deg,
            self.longitude_deg,
            self.height_km,
            self.velocity_sez_m_s().unwrap().norm()
        )
    }
}

impl fmt::LowerExp for GroundAsset {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "lat.: {:.e} deg\tlong.: {:.e} deg\talt.: {:.e} km (speed: {:.e} m/s)",
            self.latitude_deg,
            self.longitude_deg,
            self.height_km,
            self.velocity_sez_m_s().unwrap().norm()
        )
    }
}

impl State for GroundAsset {
    // State is purely lat/long/alt
    type Size = Const<6>;
    type VecLength = Const<{ 6 + 36 }>;

    fn to_vector(&self) -> nalgebra::OVector<f64, Self::VecLength> {
        let mut vector = OVector::<f64, Const<42>>::zeros();
        vector[0] = self.latitude_deg;
        vector[1] = self.longitude_deg;
        vector[2] = self.height_km;

        vector[3] = self.latitude_rate_deg_s;
        vector[4] = self.longitude_rate_deg_s;
        vector[5] = self.height_rate_km_s;

        if let Some(stm) = self.stm {
            let stm_slice = stm.as_slice();
            for i in 0..36 {
                vector[6 + i] = stm_slice[i];
            }
        }

        vector
    }

    fn orbit(&self) -> Orbit {
        Orbit::try_latlongalt(
            self.latitude_deg,
            self.longitude_deg,
            self.height_km,
            self.epoch,
            self.frame,
        )
        .unwrap()
    }

    fn epoch(&self) -> Epoch {
        self.epoch
    }

    fn set_epoch(&mut self, epoch: Epoch) {
        self.epoch = epoch;
    }

    /// Vector is expected to be organized as such:
    /// [Latitude, Longitude, Height, ]
    fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Const<42>>) {
        let asset_state =
            OVector::<f64, Self::Size>::from_column_slice(&vector.as_slice()[..Self::Size::dim()]);

        if self.stm.is_some() {
            let sc_full_stm = OMatrix::<f64, Self::Size, Self::Size>::from_column_slice(
                &vector.as_slice()[Self::Size::dim()..],
            );

            self.stm = Some(sc_full_stm);
        }

        self.latitude_deg = asset_state[0];
        self.longitude_deg = asset_state[1];
        self.height_km = asset_state[2];

        self.latitude_rate_deg_s = asset_state[3];
        self.longitude_rate_deg_s = asset_state[4];
        self.height_rate_km_s = asset_state[5];

        self.epoch = epoch;
    }

    fn reset_stm(&mut self) {
        self.stm = Some(OMatrix::<f64, Const<6>, Const<6>>::identity());
    }

    fn unset_stm(&mut self) {
        self.stm = None;
    }

    fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, DynamicsError> {
        match self.stm {
            Some(stm) => Ok(stm),
            None => Err(DynamicsError::StateTransitionMatrixUnset),
        }
    }

    fn with_stm(mut self) -> Self {
        self.stm = Some(OMatrix::<f64, Const<6>, Const<6>>::identity());
        self
    }

    fn add(self, other: OVector<f64, Self::Size>) -> Self {
        self + other
    }
}

impl Interpolatable for GroundAsset {
    fn interpolate(mut self, epoch: Epoch, states: &[Self]) -> Result<Self, InterpolationError> {
        // Interpolate the Orbit first
        // Statically allocated arrays of the maximum number of samples
        let mut epochs_tdb = [0.0; INTERPOLATION_SAMPLES];
        let mut xs = [0.0; INTERPOLATION_SAMPLES];
        let mut ys = [0.0; INTERPOLATION_SAMPLES];
        let mut zs = [0.0; INTERPOLATION_SAMPLES];
        let mut vxs = [0.0; INTERPOLATION_SAMPLES];
        let mut vys = [0.0; INTERPOLATION_SAMPLES];
        let mut vzs = [0.0; INTERPOLATION_SAMPLES];

        for (cno, state) in states.iter().enumerate() {
            xs[cno] = state.latitude_deg;
            ys[cno] = state.longitude_deg;
            zs[cno] = state.height_km;
            vxs[cno] = state.latitude_rate_deg_s;
            vys[cno] = state.longitude_rate_deg_s;
            vzs[cno] = state.height_rate_km_s;
            epochs_tdb[cno] = state.epoch.to_et_seconds();
        }

        // Ensure that if we don't have enough states, we only interpolate using what we have instead of INTERPOLATION_SAMPLES
        let n = states.len();

        let (latitude_deg, latitude_vel_deg_s) =
            hermite_eval(&epochs_tdb[..n], &xs[..n], &vxs[..n], epoch.to_et_seconds())?;

        let (longitude_deg, longitude_vel_deg_s) =
            hermite_eval(&epochs_tdb[..n], &ys[..n], &vys[..n], epoch.to_et_seconds())?;

        let (height_km, height_vel_km_s) =
            hermite_eval(&epochs_tdb[..n], &zs[..n], &vzs[..n], epoch.to_et_seconds())?;

        self.latitude_deg = latitude_deg;
        self.longitude_deg = longitude_deg;
        self.height_km = height_km;
        self.latitude_rate_deg_s = latitude_vel_deg_s;
        self.longitude_rate_deg_s = longitude_vel_deg_s;
        self.height_rate_km_s = height_vel_km_s;

        self.epoch = epoch;

        Ok(self)
    }

    fn frame(&self) -> Frame {
        self.frame
    }

    fn set_frame(&mut self, frame: Frame) {
        self.frame = frame;
    }

    fn export_params() -> Vec<StateParameter> {
        vec![
            StateParameter::Element(OrbitalElement::X),
            StateParameter::Element(OrbitalElement::Y),
            StateParameter::Element(OrbitalElement::Z),
            StateParameter::Element(OrbitalElement::VX),
            StateParameter::Element(OrbitalElement::VY),
            StateParameter::Element(OrbitalElement::VZ),
        ]
    }
}

impl Add<OVector<f64, Const<6>>> for GroundAsset {
    type Output = Self;

    /// Adds the provided state deviation to this orbit
    fn add(mut self, asset_state: OVector<f64, Const<6>>) -> Self {
        self.latitude_deg += asset_state[0];
        self.longitude_deg += asset_state[1];
        self.height_km += asset_state[2];

        self.latitude_rate_deg_s += asset_state[3];
        self.longitude_rate_deg_s += asset_state[4];
        self.height_rate_km_s += asset_state[5];

        self
    }
}

pub fn latlongalt_rate(
    orbit: Orbit,
    velocity_sez_km_s: Vector3<f64>,
) -> Result<(f64, f64, f64), PhysicsError> {
    // Get current lat, long, alt
    let (lat_deg, _long_deg, alt_km) = orbit.latlongalt()?;
    let lat_rad = lat_deg.to_radians();

    // // Get the DCM from SEZ to body-fixed frame
    // let sez2body_dcm = orbit.dcm_from_topocentric_to_body_fixed()?;

    // Extract SEZ velocity components
    // SEZ frame: x=South, y=East, z=Zenith (up)
    let v_south = velocity_sez_km_s.x;
    let v_east = velocity_sez_km_s.y;
    let v_zenith = velocity_sez_km_s.z;

    // Get ellipsoid parameters
    let a_km = orbit.frame.mean_equatorial_radius_km()?;
    let b_km = orbit.frame.shape.unwrap().polar_radius_km;
    let e2 = (a_km.powi(2) - b_km.powi(2)) / a_km.powi(2);

    // Compute radius of curvature in the meridian (N)
    let sin_lat = lat_rad.sin();
    let n = a_km / (1.0 - e2 * sin_lat.powi(2)).sqrt();

    // Compute radius of curvature in the prime vertical (M)
    let m = a_km * (1.0 - e2) / (1.0 - e2 * sin_lat.powi(2)).powf(1.5);

    // Convert SEZ velocities to geodetic coordinate rates
    // Altitude rate (positive = up)
    let alt_rate_km_s = v_zenith;

    // Latitude rate (positive = north, but SEZ x-axis is south)
    let lat_rate_rad_s = -v_south / (m + alt_km);
    let lat_rate_deg_s = lat_rate_rad_s.to_degrees();

    // Longitude rate (positive = east)
    let cos_lat = lat_rad.cos();
    let long_rate_rad_s = if cos_lat.abs() > 1e-10 {
        v_east / ((n + alt_km) * cos_lat)
    } else {
        0.0 // At poles, longitude rate is undefined
    };
    let long_rate_deg_s = long_rate_rad_s.to_degrees();

    Ok((lat_rate_deg_s, long_rate_deg_s, alt_rate_km_s))
}

/// Convert geodetic coordinate rates (lat rate, long rate, alt rate) to velocity in SEZ frame
pub fn velocity_sez_from_latlongalt_rate(
    orbit: Orbit,
    lat_rate_deg_s: f64,
    long_rate_deg_s: f64,
    alt_rate_km_s: f64,
) -> Result<Vector3<f64>, PhysicsError> {
    // Get current lat and alt
    let (lat_deg, _long_deg, alt_km) = orbit.latlongalt()?;
    let lat_rad = lat_deg.to_radians();

    // Get ellipsoid parameters
    let a_km = orbit.frame.mean_equatorial_radius_km()?;
    let b_km = orbit.frame.shape.unwrap().polar_radius_km;
    let e2 = (a_km.powi(2) - b_km.powi(2)) / a_km.powi(2);

    // Compute radius of curvature in the meridian (M)
    let sin_lat = lat_rad.sin();
    let m = a_km * (1.0 - e2) / (1.0 - e2 * sin_lat.powi(2)).powf(1.5);

    // Compute radius of curvature in the prime vertical (N)
    let n = a_km / (1.0 - e2 * sin_lat.powi(2)).sqrt();

    // Convert rates to radians per second
    let lat_rate_rad_s = lat_rate_deg_s.to_radians();
    let long_rate_rad_s = long_rate_deg_s.to_radians();

    // Convert geodetic coordinate rates to SEZ velocities
    // SEZ frame: x=South, y=East, z=Zenith (up)
    let v_south = -lat_rate_rad_s * (m + alt_km); // Negative because x is south
    let v_east = long_rate_rad_s * (n + alt_km) * lat_rad.cos();
    let v_zenith = alt_rate_km_s;

    Ok(Vector3::new(v_south, v_east, v_zenith))
}
