extern crate hifitime;
extern crate nalgebra as na;

use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, Dim, DimName, MatrixNM, VectorN};
use celestia::State;
use hifitime::instant::Instant;

pub mod kalman;
pub mod ranging;

pub trait Linearization
where
    Self: Sized,
{
    /// Defines the state size of the estimated state
    type StateSize: Dim + DimName;
    /// Defines the gradient of the equations of motion for these dynamics.
    fn gradient(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> MatrixNM<f64, Self::StateSize, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize>;
}

pub trait Measurement
where
    Self: Sized,
{
    /// Defines the state size of the estimated state
    type StateSize: Dim + DimName;
    /// Defines how much data is measured. For example, if measuring range and range rate, this should be of size 2 (nalgebra::U2).
    type MeasurementSize: Dim + DimName;

    /// Computes a new measurement from the provided information.
    fn new(tx: State, rx: State, elevation_mask: f64, noise: Vector<f64, Self::MeasurementSize>) -> Self;

    /// Returns the measurement/observation as a vector.
    fn observation(&self) -> &Vector<f64, Self::MeasurementSize>;

    /// Returns the measurement sensitivity (often referred to as H tilde).
    fn sensitivity(&self) -> &MatrixNM<f64, Self::MeasurementSize, Self::StateSize>;

    /// Returns whether the transmitter and receiver where in line of sight.
    fn visible(&self) -> bool;

    /// Returns the time at which the measurement was performed.
    fn at(&self) -> Instant;
}
