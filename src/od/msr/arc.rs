extern crate nalgebra;
use self::nalgebra::vector;
use super::{ObservationOf, Observer};
use crate::linalg::allocator::Allocator;
use crate::linalg::{Const, DefaultAllocator, DimName, OMatrix, SVector};
use crate::State;
use std::marker::PhantomData;

pub struct MeasurementArc<S: State, Src: Observer<S>>
where
    DefaultAllocator: Allocator<f64, S::Size>
        + Allocator<f64, S::VecLength>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<f64, Const<1>, S::Size>,
{
    source: Src,
    _state: PhantomData<S>,
}

#[test]
fn name() {
    let zero: SVector<f64, 0> = vector![];
    println!("{:?}", 2.0 * zero);
}
