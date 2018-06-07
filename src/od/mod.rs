extern crate hifitime;
extern crate nalgebra as na;

use self::hifitime::instant::Instant;
use self::na::allocator::Allocator;
use self::na::{DefaultAllocator, Dim, DimName, MatrixMN, VectorN};
use celestia::{CoordinateFrame, State};

pub mod kalman;
pub mod ranging;

pub trait Linearization
where
    Self: Sized,
{
    /// Defines the state size of the estimated state
    type StateSize: Dim + DimName;
    /// Defines the gradient of the equations of motion for these dynamics.
    fn gradient(&self, t: f64, state: &VectorN<f64, Self::StateSize>) -> MatrixMN<f64, Self::StateSize, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize> + Allocator<f64, Self::StateSize, Self::StateSize>;
}

pub trait Measurement
where
    Self: Sized,
    DefaultAllocator: Allocator<f64, Self::MeasurementSize> + Allocator<f64, Self::MeasurementSize, Self::StateSize>,
{
    /// Defines the state size of the estimated state
    type StateSize: Dim + DimName;
    /// Defines how much data is measured. For example, if measuring range and range rate, this should be of size 2 (nalgebra::U2).
    type MeasurementSize: Dim + DimName;

    /// Computes a new measurement from the provided information.
    fn new<F: CoordinateFrame>(
        dt: Instant,
        tx: State<F>,
        rx: State<F>,
        obs: VectorN<f64, Self::MeasurementSize>,
        visible: bool,
    ) -> Self;

    /// Returns the measurement/observation as a vector.
    fn observation(&self) -> &VectorN<f64, Self::MeasurementSize>
    where
        DefaultAllocator: Allocator<f64, Self::MeasurementSize>;

    /// Returns the measurement sensitivity (often referred to as H tilde).
    fn sensitivity(&self) -> &MatrixMN<f64, Self::MeasurementSize, Self::StateSize>
    where
        DefaultAllocator: Allocator<f64, Self::StateSize, Self::MeasurementSize>;

    /// Returns whether the transmitter and receiver where in line of sight.
    fn visible(&self) -> bool;

    /// Returns the time at which the measurement was performed.
    fn at(&self) -> Instant;
}
