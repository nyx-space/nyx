use crate::celestia::{Frame, Orbit};
use crate::dimensions::Vector3;
use crate::time::{Epoch, TimeUnit};
use std::f64::consts::FRAC_PI_2 as half_pi;

#[derive(Debug)]
pub enum ThrustingError {
    NoObjectiveDefined,
}

/// The `ThrustControl` trait handles control laws, optimizations, and other such methods for
/// controlling the overall thrust direction when tied to a `Spacecraft`. For delta V control,
/// tie the DeltaVctrl to a MissionArc.
pub trait ThrustControl: Send + Sync {
    /// Returns a unit vector corresponding to the thrust direction in the inertial frame.
    fn direction(&self, state: &Orbit) -> Vector3<f64>;

    /// Returns a number between [0;1] corresponding to the engine throttle level.
    /// For example, 0 means coasting, i.e. no thrusting, and 1 means maximum thrusting.
    fn throttle(&self, state: &Orbit) -> f64;

    /// Prepares the controller for the next maneuver (called from set_state of the dynamics).
    fn next(&mut self, state: &Orbit);

    /// Returns whether this thrust control has been achieve, if it has an objective
    fn achieved(&self, _state: &Orbit) -> Result<bool, ThrustingError> {
        Err(ThrustingError::NoObjectiveDefined)
    }
}

/// Goals used for sub-optimal controls
#[derive(Copy, Clone, Debug)]
pub enum Achieve {
    Sma { target: f64, tol: f64 },
    Ecc { target: f64, tol: f64 },
    Inc { target: f64, tol: f64 },
    Raan { target: f64, tol: f64 },
    Aop { target: f64, tol: f64 },
}

impl Achieve {
    pub fn achieved(&self, state: &Orbit) -> bool {
        match *self {
            Achieve::Sma { target, tol } => (state.sma() - target).abs() < tol,
            Achieve::Ecc { target, tol } => (state.ecc() - target).abs() < tol,
            Achieve::Inc { target, tol } => (state.inc() - target).abs() < tol,
            Achieve::Raan { target, tol } => (state.raan() - target).abs() < tol,
            Achieve::Aop { target, tol } => (state.aop() - target).abs() < tol,
        }
    }
}

/// Mnvr defined a single maneuver. Direction MUST be in the VNC frame (Velocity / Normal / Cross).
/// It may be used with a maneuver scheduler.
#[derive(Copy, Clone, Debug)]
pub struct Mnvr {
    /// Start epoch of the maneuver
    pub start: Epoch,
    /// End epoch of the maneuver
    pub end: Epoch,
    /// Thrust level, if 1.0 use all thruster available at full power
    pub thrust_lvl: f64,
    /// Direction of the thrust in the VNC frame
    pub vector: Vector3<f64>,
}

impl Mnvr {
    /// Creates an instantaneous maneuver whose vector is the deltaV.
    pub fn instantaneous(dt: Epoch, vector: Vector3<f64>) -> Self {
        Self {
            start: dt,
            end: dt + TimeUnit::Microsecond,
            thrust_lvl: 1.0,
            vector,
        }
    }
}

/// Ruggiero defines the closed loop control law from IEPC 2011-102
#[derive(Copy, Clone, Debug)]
pub struct Ruggiero {
    /// Stores the objectives
    objectives: [Option<Achieve>; 5],
    init_state: Orbit,
    achieved: bool,
}

/// The Ruggiero is a locally optimal control of a state for specific osculating elements.
/// WARNING: Objectives must be in degrees!
impl Ruggiero {
    pub fn new(objectives: Vec<Achieve>, initial: Orbit) -> Self {
        let mut objs: [Option<Achieve>; 5] = [None, None, None, None, None];
        for i in 0..objectives.len() {
            objs[i] = Some(objectives[i]);
        }
        Self {
            objectives: objs,
            init_state: initial,
            achieved: false,
        }
    }

    fn weighting(init: f64, target: f64, osc: f64, tol: f64) -> f64 {
        if (osc - target).abs() < tol {
            0.0
        } else {
            // Let's add the tolerance to the initial value if we want to keep a parameter fixed (i.e. target and initial are equal)
            (target - osc)
                / (target
                    - if (init - target).abs() < tol {
                        init + tol
                    } else {
                        init
                    })
                .abs()
        }
    }
}

impl ThrustControl for Ruggiero {
    /// Returns whether the control law has achieved all goals
    fn achieved(&self, state: &Orbit) -> Result<bool, ThrustingError> {
        for maybe_obj in &self.objectives {
            if let Some(obj) = maybe_obj {
                if !obj.achieved(state) {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }

    fn direction(&self, osc: &Orbit) -> Vector3<f64> {
        if self.achieved {
            Vector3::zeros()
        } else {
            let mut ctrl = Vector3::zeros();
            for maybe_obj in &self.objectives {
                if let Some(obj) = maybe_obj {
                    match *obj {
                        Achieve::Sma { target, tol } => {
                            let weight =
                                Self::weighting(self.init_state.sma(), target, osc.sma(), tol);
                            if weight.abs() > 0.0 {
                                let num = osc.ecc() * osc.ta().to_radians().sin();
                                let denom = 1.0 + osc.ecc() * osc.ta().to_radians().cos();
                                let alpha = num.atan2(denom);
                                // TODO: Add efficiency
                                ctrl += unit_vector_from_angles(alpha, 0.0) * weight;
                            }
                        }
                        Achieve::Ecc { target, tol } => {
                            let weight =
                                Self::weighting(self.init_state.ecc(), target, osc.ecc(), tol);
                            if weight.abs() > 0.0 {
                                let num = osc.ta().to_radians().sin();
                                let denom =
                                    osc.ta().to_radians().cos() + osc.ea().to_radians().cos();
                                let alpha = num.atan2(denom);
                                // TODO: Add efficiency
                                ctrl += unit_vector_from_angles(alpha, 0.0) * weight;
                            }
                        }
                        Achieve::Inc { target, tol } => {
                            let weight =
                                Self::weighting(self.init_state.inc(), target, osc.inc(), tol);
                            if weight.abs() > 0.0 {
                                let beta =
                                    half_pi.copysign(((osc.ta() + osc.aop()).to_radians()).cos());
                                // TODO: Add efficiency
                                ctrl += unit_vector_from_angles(0.0, beta) * weight;
                            }
                        }
                        Achieve::Raan { target, tol } => {
                            // BUG: https://gitlab.com/chrisrabotin/nyx/issues/83
                            let weight =
                                Self::weighting(self.init_state.raan(), target, osc.raan(), tol);
                            if weight.abs() > 0.0 {
                                let beta =
                                    half_pi.copysign(((osc.ta() + osc.aop()).to_radians()).sin());
                                // TODO: Add efficiency
                                ctrl += unit_vector_from_angles(0.0, beta) * weight;
                            }
                        }
                        Achieve::Aop { target, tol } => {
                            let weight =
                                Self::weighting(self.init_state.aop(), target, osc.aop(), tol);
                            let oe2 = 1.0 - osc.ecc().powi(2);
                            let e3 = osc.ecc().powi(3);
                            // Compute the optimal true anomaly for in-plane thrusting
                            let sqrt_val = (0.25 * (oe2 / e3).powi(2) + 1.0 / 27.0).sqrt();
                            let opti_ta_alpha = ((oe2 / (2.0 * e3) + sqrt_val).powf(1.0 / 3.0)
                                - (-oe2 / (2.0 * e3) + sqrt_val).powf(1.0 / 3.0)
                                - 1.0 / osc.ecc())
                            .acos();
                            // Compute the optimal true anomaly for out of plane thrusting
                            let opti_ta_beta = (-osc.ecc() * osc.aop().to_radians().cos()).acos()
                                - osc.aop().to_radians();
                            // And choose whether to do an in-plane or out of plane thrust
                            if (osc.ta().to_radians() - opti_ta_alpha).abs()
                                < (osc.ta().to_radians() - opti_ta_beta).abs()
                            {
                                // In plane
                                let p = osc.semi_parameter();
                                let (sin_ta, cos_ta) = osc.ta().to_radians().sin_cos();
                                let alpha = (-p * cos_ta).atan2((p + osc.rmag()) * sin_ta);
                                ctrl += unit_vector_from_angles(alpha, 0.0) * weight;
                            } else {
                                // Out of plane
                                let beta = half_pi.copysign(
                                    -(osc.ta().to_radians() + osc.aop().to_radians()).sin(),
                                ) * osc.inc().to_radians().cos();
                                ctrl += unit_vector_from_angles(0.0, beta) * weight;
                            };
                        }
                    }
                }
            }
            // Return a normalized vector
            ctrl = if ctrl.norm() > 0.0 {
                ctrl / ctrl.norm()
            } else {
                ctrl
            };
            // Convert to inertial
            osc.dcm_to_inertial(Frame::RCN) * ctrl
        }
    }

    // Either thrust full power or not at all
    fn throttle(&self, osc: &Orbit) -> f64 {
        if self.achieved {
            0.0
        } else {
            for maybe_obj in &self.objectives {
                if let Some(obj) = maybe_obj {
                    match *obj {
                        Achieve::Sma { target, tol } => {
                            let weight =
                                Self::weighting(self.init_state.sma(), target, osc.sma(), tol);
                            if weight.abs() > 0.0 {
                                return 1.0;
                            }
                        }
                        Achieve::Ecc { target, tol } => {
                            let weight =
                                Self::weighting(self.init_state.ecc(), target, osc.ecc(), tol);
                            if weight.abs() > 0.0 {
                                return 1.0;
                            }
                        }
                        Achieve::Inc { target, tol } => {
                            let weight =
                                Self::weighting(self.init_state.inc(), target, osc.inc(), tol);
                            if weight.abs() > 0.0 {
                                return 1.0;
                            }
                        }
                        Achieve::Raan { target, tol } => {
                            let weight =
                                Self::weighting(self.init_state.raan(), target, osc.raan(), tol);
                            if weight.abs() > 0.0 {
                                return 1.0;
                            }
                        }
                        Achieve::Aop { target, tol } => {
                            let weight =
                                Self::weighting(self.init_state.aop(), target, osc.aop(), tol);
                            if weight.abs() > 0.0 {
                                return 1.0;
                            }
                        }
                    }
                }
            }
            0.0
        }
    }

    /// Update the state for the next iteration
    fn next(&mut self, osc: &Orbit) {
        if self.throttle(osc) > 0.0 {
            if self.achieved {
                info!("enabling control: {:o}", osc);
            }
            self.achieved = false;
        } else {
            if !self.achieved {
                info!("disabling control: {:o}", osc);
            }
            self.achieved = true;
        }
    }
}

/// A controller for a set of pre-determined maneuvers.
#[derive(Clone, Debug)]
pub struct FiniteBurns {
    /// Maneuvers should be provided in chronological order, first maneuver first in the list
    pub mnvrs: Vec<Mnvr>,
    pub mnvr_no: usize,
    /// The frame in which the maneuvers are defined.
    pub frame: Frame,
}

impl FiniteBurns {
    /// Builds a schedule from the vector of maneuvers, must be provided in chronological order.
    pub fn from_mnvrs(mnvrs: Vec<Mnvr>, frame: Frame) -> Self {
        assert!(
            matches!(frame, Frame::Inertial | Frame::VNC),
            "Maneuvers must be either in the inertial frame or in a body frame"
        );
        Self {
            mnvrs,
            mnvr_no: 0,
            frame,
        }
    }
}

impl ThrustControl for FiniteBurns {
    fn direction(&self, osc: &Orbit) -> Vector3<f64> {
        // NOTE: We do not increment the mnvr number here. The power function is called first,
        // so we let that function handle starting and stopping of the maneuver.
        if self.mnvr_no >= self.mnvrs.len() {
            Vector3::zeros()
        } else {
            let next_mnvr = self.mnvrs[self.mnvr_no];
            if next_mnvr.start <= osc.dt {
                if matches!(self.frame, Frame::Inertial) {
                    next_mnvr.vector
                } else {
                    osc.dcm_to_inertial(self.frame) * next_mnvr.vector
                }
            } else {
                Vector3::zeros()
            }
        }
    }

    fn throttle(&self, osc: &Orbit) -> f64 {
        if self.mnvr_no >= self.mnvrs.len() {
            0.0
        } else {
            let next_mnvr = self.mnvrs[self.mnvr_no];
            if next_mnvr.start <= osc.dt {
                next_mnvr.thrust_lvl
            } else {
                0.0
            }
        }
    }

    fn next(&mut self, osc: &Orbit) {
        if self.mnvr_no < self.mnvrs.len() {
            let cur_mnvr = self.mnvrs[self.mnvr_no];
            if osc.dt >= cur_mnvr.end {
                self.mnvr_no += 1;
            }
        }
    }
}

fn unit_vector_from_angles(alpha: f64, beta: f64) -> Vector3<f64> {
    Vector3::new(
        alpha.sin() * beta.cos(),
        alpha.cos() * beta.cos(),
        beta.sin(),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::celestia::Cosm;
    #[test]
    fn ruggiero_weight() {
        let mut cosm = Cosm::from_xb("./de438s").unwrap();
        cosm.mut_gm_for_frame("EME2000", 398_600.433);
        let eme2k = cosm.frame("EME2000");
        let start_time = Epoch::from_gregorian_tai_at_midnight(2020, 1, 1);
        let orbit = Orbit::keplerian(7378.1363, 0.01, 0.05, 0.0, 0.0, 1.0, start_time, eme2k);

        // Define the objectives
        let objectives = vec![
            Achieve::Sma {
                target: 42164.0,
                tol: 1.0,
            },
            Achieve::Ecc {
                target: 0.01,
                tol: 5e-5,
            },
        ];

        let ruggiero = Ruggiero::new(objectives, orbit);
        // 7301.597157 201.699933 0.176016 -0.202974 7.421233 0.006476 298.999726
        let osc = Orbit::cartesian(
            7_303.253_461_441_64f64,
            127.478_714_816_381_75,
            0.111_246_193_227_445_4,
            -0.128_284_025_765_195_6,
            7.422_889_151_816_439,
            0.006_477_694_429_837_2,
            start_time,
            eme2k,
        );
        let expected = Vector3::new(
            -0.017_279_636_133_108_3,
            0.999_850_315_226_803,
            0.000_872_534_222_883_2,
        );

        let got = ruggiero.direction(&osc);

        println!("{}", expected - got);
        assert!(
            (expected - got).norm() < 1e-12,
            "incorrect direction computed"
        );
    }
}
