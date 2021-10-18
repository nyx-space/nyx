use crate::cosmic::{Orbit, Spacecraft};
use crate::linalg::allocator::Allocator;
use crate::linalg::{
    Const, DefaultAllocator, DimAdd, DimMul, DimName, DimProd, DimSum, OMatrix, OVector, SVector,
};
use crate::time::{Duration, Epoch};
use crate::{NyxError, State};
use std::fmt;

#[derive(Copy, Clone, Debug, PartialEq)]
pub struct NoParameters {}

impl fmt::Display for NoParameters {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[no estimated parameters]",)
    }
}

impl fmt::LowerExp for NoParameters {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(f, "[no estimated parameters]",)
    }
}

impl State for NoParameters {
    type Size = Const<0>;
    type VecLength = Const<0>;

    fn reset_stm(&mut self) {}

    fn zeros() -> Self {
        Self {}
    }

    fn as_vector(&self) -> Result<OVector<f64, Const<0>>, NyxError> {
        Ok(OVector::<f64, Const<0>>::zeros())
    }

    fn set(&mut self, _: Epoch, _: &OVector<f64, Const<0>>) -> Result<(), NyxError> {
        Ok(())
    }

    fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, NyxError> {
        Ok(OMatrix::<f64, Const<0>, Const<0>>::zeros())
    }

    fn epoch(&self) -> Epoch {
        unimplemented!();
    }

    fn set_epoch(&mut self, _: Epoch) {}

    fn add(self, _: OVector<f64, Self::Size>) -> Self {
        self
    }
}

pub struct ClockBiasDrift {
    bias: f64,
    drift: f64,
    epoch: Epoch,
}
