use super::hifitime::SECONDS_PER_DAY;
use super::Dynamics;
use celestia::{CelestialBody, EARTH};
use std::collections::HashMap;

/// `TwoBody` exposes the equations of motion for a simple two body propagation.
#[derive(Copy, Clone)]
pub struct Harmonics {
    julian_days: f64,
    mu: f64,
    radius: f64,
    order: u16,
    degree: u16,
    // data is (degree, order) -> (C_nm, S_nm)
    data: HashMap<(u16, u16), (f64, f64)>,
}

impl Harmonics {
    /// Initialize `Harmonics` as J<sub>2</sub> only using the JGM3 model (available in GMAT)
    ///
    /// Use the embedded Earth parameter. If others are needed, low from `from_gunzip`.
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn J2_JGM3(start_julian_days: f64) -> Harmonics {
        let mut data = HashMap::new();
        data.insert((2, 0), (-4.84165374886470e-04, 0.0));
        Harmonics {
            julian_days: start_julian_days,
            mu: EARTH::gm(),
            radius: EARTH::eq_radius(),
            order: 2,
            degree: 0,
            data: data,
        }
    }

    /// Initialize `Harmonics` as J<sub>2</sub> only using the EGM2008 model (from the GRACE mission, best model as of 2018)
    ///
    /// *WARNING:* This is an EARTH gravity model, and _should not_ be used around any other body.
    pub fn J2_EGM2008(start_julian_days: f64) -> Harmonics {
        let mut data = HashMap::new();
        data.insert((2, 0), (-0.484165143790815e-03, 0.0));
        Harmonics {
            julian_days: start_julian_days,
            mu: EARTH::gm(),
            radius: EARTH::eq_radius(),
            order: 2,
            degree: 0,
            data: data,
        }
    }

    /// Initialize `Harmonics` from the file path (must be a gunzipped file)
    ///
    /// Gravity models provided by `nyx`: TODO: add github links and examples
    /// + EMG2008 to 2190 for Earth (tide free)
    /// + Moon to 1500 (from SHADR file)
    /// + Mars to 120 (from SHADR file)
    /// + Venus to 150 (from SHADR file)
    pub fn from_gunzip(
        start_julian_days: f64,
        degree: u16,
        order: u16,
        filepath: String,
    ) -> Harmonics {
        unimplemented!();
    }
}

impl Dynamics for Harmonics {
    type StateSize = U6;

    fn time(&self) -> f64 {
        0.0
    }

    fn state(&self) -> VectorN<f64, Self::StateSize> {
        // No state is associated with Harmonics, always return zero
        Vector6::zero()
    }

    /// WARNING: The *new_t* parameter is considered a TIME INCREASE from the initial julian days.
    /// This is likely a bad assumption (breaks back prop unless new_t is negative) -- I will need to clarify that in the docs.
    fn set_state(&mut self, new_t: f64, _new_state: &VectorN<f64, Self::StateSize>) {
        self.time += new_t / SECONDS_PER_DAY;
    }

    /// WARNING: This provides a DELTA of the state, which must be added to the result of the TwoBody propagator being used.
    /// However, the provided `state` must be the position and velocity.
    fn eom(&self, _t: f64, state: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let radius = state.fixed_slice::<U3, U1>(0, 0);
        // Using the GMAT notation, with extra character for ease of highlight
        let r_ = radius.norm();
        let s_ = radius[(0, 0)] / radius.norm();
        let t_ = radius[(1, 0)] / radius.norm();
        let u_ = radius[(2, 0)] / radius.norm();
        // TODO: This
        let velocity = state.fixed_slice::<U3, U1>(3, 0);
        let body_acceleration = (-self.mu / radius.norm().powi(3)) * radius;
        Vector6::from_iterator(velocity.iter().chain(body_acceleration.iter()).cloned())
    }
}
