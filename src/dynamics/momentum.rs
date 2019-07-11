use super::na::{Matrix3, Vector3, VectorN, U3};
use super::Dynamics;
use std::f64;
use utils::is_diagonal;
// use alga::general::operator::Inverse;

/// `AngularMom` exposes the equations of motion for the angular momentum of a **rigid body**.
#[derive(Copy, Clone, Debug)]
pub struct AngularMom {
    time: f64, // Needed to stop the integration by calling dyn.time()
    tensor: Matrix3<f64>,
    velocity: Vector3<f64>,
}

impl AngularMom {
    /// Initializes a new AngularMom struct with an time-invariant inertia tensor of a rigid body and its original angular velocity.
    ///
    /// Throughout this documentation [I] refers to the inertia tensor and ω to the angular velocity.
    /// NOTE: The provided inertia tensor **must** be expressed in a frame such that it is a diagonal matrix.
    /// There is always such a frame. If this is _not_ the primary frame desired, use the parallel axis theorem.
    /// This theorem is developped in Schaub & Junkins, "Analytical Mechanics of Space Systems", 3th ed., page 163,
    /// section 4.2.2 "Inertia Matrix Properties". This function will **panic!** if the inertia tensor is not diagonal.
    pub fn from_tensor_matrix(tensor: &Matrix3<f64>, velocity: &Vector3<f64>) -> AngularMom {
        if !is_diagonal(tensor) {
            panic!("The provided inertia tensor is not diagonal.");
        }
        AngularMom {
            time: 0.0,
            tensor: *tensor,
            velocity: *velocity,
        }
    }

    /// Returns the angular momentum of the system, i.e. [I]ω
    pub fn momentum(&self) -> Vector3<f64> {
        self.tensor * self.velocity
    }
}

impl Dynamics for AngularMom {
    type StateSize = U3;

    fn time(&self) -> f64 {
        self.time
    }

    /// Returns the **angular velocity** ω of the system, not its momentum.
    fn state(&self) -> VectorN<f64, Self::StateSize> {
        self.velocity
    }

    /// Set the **angular velocity** ω of the system and the time.
    fn set_state(&mut self, new_t: f64, new_angular_velocity: &VectorN<f64, Self::StateSize>) {
        self.time = new_t;
        self.velocity = *new_angular_velocity;
    }

    /// Computes the instantaneous equations of motion of the angular velocity of a tensor (i.e. the angular acceleration).
    /// [I]̲̇ω = -[̃ω][I]̲ω + ̲L
    ///
    /// Source: Schaub & Junkins, 3th ed., eq. 4.32.
    fn eom(&self, _t: f64, omega: &VectorN<f64, Self::StateSize>) -> VectorN<f64, Self::StateSize> {
        let omega_dot_x = (self.tensor[(1, 1)] - self.tensor[(2, 2)]) / self.tensor[(0, 0)] * omega[(1, 0)] * omega[(2, 0)];
        let omega_dot_y = (self.tensor[(2, 2)] - self.tensor[(0, 0)]) / self.tensor[(1, 1)] * omega[(2, 0)] * omega[(0, 0)];
        let omega_dot_z = (self.tensor[(0, 0)] - self.tensor[(1, 1)]) / self.tensor[(2, 2)] * omega[(0, 0)] * omega[(1, 0)];
        Vector3::new(omega_dot_x, omega_dot_y, omega_dot_z)
    }
}
