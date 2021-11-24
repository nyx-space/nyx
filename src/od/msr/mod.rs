extern crate rand_pcg;

// use super::{Measurement, MeasurementDevice};
use crate::linalg::allocator::Allocator;
use crate::linalg::{Const, DMatrix, DVector, DefaultAllocator, DimName, OMatrix, SVector};
use crate::{time::Epoch, NyxError, State, TimeTagged};

pub mod gs_setup;

/// Measurement Arc
pub mod arc;

pub mod sim;

/*

// How to use

let mut msp = MsrSim {
    seed: 0 (u64),
    light_time: LightTimeCalc::None,
    arcs: Vec<Arc>,
};


*/

/// A trait defining a measurement where Observable is the state that can be observed (e.g. Orbit or Spacecraft)
/// and MeasurementSize is the size of the measurement (e.g. Const<1> for TwoWayDoppler or Const<2> for TwoWayRangeRangeRate).
pub trait ObservationOf<Observable: State>: TimeTagged
where
    Observable: State,
    DefaultAllocator: Allocator<f64, Observable::Size>
        + Allocator<f64, Observable::Size, Observable::Size>
        + Allocator<f64, Const<1>, Observable::Size>
        + Allocator<f64, Observable::VecLength>,
{
    /// Returns the measurement/observation as a vector, if partner is visible
    fn observation(&self) -> Option<f64>;

    /// Returns the measurement sensitivity (often referred to as H tilde)
    fn sensitivity(&self) -> Option<OMatrix<f64, Const<1>, Observable::Size>>;
}

/// A trait to generalize measurement devices such as a ground station
pub trait Observer<Observable>
where
    Observable: State,
    DefaultAllocator: Allocator<f64, Observable::Size>
        + Allocator<f64, Observable::Size, Observable::Size>
        + Allocator<f64, Const<1>, Observable::Size>
        + Allocator<f64, Observable::VecLength>,
{
    /// Returns the measurement if the device can measure the observable (visibility), or an error if the computation failed
    fn measure(
        &self,
        input: &Observable,
        kind: MeasurementKind,
    ) -> Result<Option<Measurement<Observable>>, NyxError>;
}

pub struct Measurement<Observable: State>
where
    Observable: State,
    DefaultAllocator: Allocator<f64, Observable::Size>
        + Allocator<f64, Observable::Size, Observable::Size>
        + Allocator<f64, Const<1>, Observable::Size>
        + Allocator<f64, Observable::VecLength>,
{
    pub epoch: Epoch,
    pub kind: MeasurementKind,
    pub value: f64,
    pub sensitivity: OMatrix<f64, Const<1>, Observable::Size>,
}

/// The problem with this design is that only some measurements are supported, and I'm the king of those... Or, I could add a bunch of stuff (everything from TDM plus IMU readings) but not implement them here.
#[derive(Copy, Clone, Debug)]
pub enum MeasurementKind {
    OneWayRange,
    OneWayDoppler,
    CustomMeasurementKind(u8),
}
