use super::{ObservationOf, Observer};
use crate::linalg::allocator::Allocator;
use crate::linalg::{Const, DVector, DefaultAllocator, Dim, DimName, OMatrix, OVector, SVector};
use crate::State;
use std::sync::Arc;

/// This structure allows the simulation of measurements from a list of observers, their respective cadences, hand-over strategy between stations, noise models, and the seed.
/// A new measurement arc is generated whenever there is a handoff between one observer to the next. Observers do _not_ need to support the same observation types.
pub struct MsrSim<S: State>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::VecLength>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<f64, Const<1>, S::Size>,
{
    participants: Vec<Arc<dyn Observer<S>>>,
}
