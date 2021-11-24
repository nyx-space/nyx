use super::parameters::NoParameters;
use crate::cosmic::{Orbit, Spacecraft};
use crate::linalg::allocator::Allocator;
use crate::linalg::{
    Const, DefaultAllocator, DimAdd, DimMul, DimName, DimProd, DimSum, OMatrix, OVector, SVector,
};
use crate::time::{Duration, Epoch};
use crate::{NyxError, State};
use std::fmt;
use std::marker::PhantomData;
use std::sync::Arc;

/// The EmbedState trait allows specifying that a super state contains several sub-states, e.g. a Spacecraft contains an Orbit and a Quaternion.
pub trait EmbedState<O: State>
where
    Self: State,
    DefaultAllocator: Allocator<f64, O::Size>
        + Allocator<f64, O::VecLength>
        + Allocator<f64, O::Size, O::Size>
        + Allocator<f64, Self::Size>
        + Allocator<f64, Self::VecLength>
        + Allocator<f64, Self::Size, Self::Size>,
{
    fn embed(&mut self, other: O);
    fn extract(&self) -> O;
}

impl<S: State> EmbedState<S> for S
where
    DefaultAllocator: Allocator<f64, Self::Size>
        + Allocator<f64, Self::VecLength>
        + Allocator<f64, Self::Size, Self::Size>,
{
    fn extract(&self) -> S {
        *self
    }

    fn embed(&mut self, s: S) {
        *self = s;
    }
}

// Yet another reason why I would want a Parameter trait! The following "overwrites" the generic impl above, but it should not really.`
// impl<S: State> EmbedState<NoParameters> for S
// where
//     DefaultAllocator: Allocator<f64, Self::Size>
//         + Allocator<f64, Self::VecLength>
//         + Allocator<f64, Self::Size, Self::Size>,
// {
//     fn extract(&self) -> NoParameters {
//         NoParameters::zeros()
//     }

//     fn embed(&mut self, _: NoParameters) {}
// }

pub trait NavState<X: State, P: State>
where
    Self: EmbedState<X> + EmbedState<P>,
    DefaultAllocator: Allocator<f64, X::Size>
        + Allocator<f64, X::VecLength>
        + Allocator<f64, X::Size, X::Size>
        + Allocator<f64, P::Size>
        + Allocator<f64, P::VecLength>
        + Allocator<f64, P::Size, P::Size>
        + Allocator<f64, Self::Size>
        + Allocator<f64, Self::VecLength>
        + Allocator<f64, Self::Size, Self::Size>,
{
    /// Extract the estimated physical state
    fn state(&self) -> X {
        self.extract()
    }

    /// Extract the estimated parameters
    fn parameters(&self) -> P {
        self.extract()
    }
}

/**** End trait definitions ****/

impl EmbedState<Orbit> for Spacecraft {
    fn extract(&self) -> Orbit {
        self.orbit
    }

    fn embed(&mut self, orbit: Orbit) {
        self.orbit = orbit;
    }
}

impl EmbedState<NoParameters> for Orbit {
    fn extract(&self) -> NoParameters {
        NoParameters::zeros()
    }

    fn embed(&mut self, _: NoParameters) {}
}

pub type OrbitNavState = Orbit;

impl NavState<Orbit, NoParameters> for OrbitNavState {}

// #[derive(Copy, Clone, Debug, PartialEq)]
// pub struct NavState<X: State, P: State>
// where
//     DefaultAllocator: Allocator<f64, X::Size>
//         + Allocator<f64, X::VecLength>
//         + Allocator<f64, X::Size, X::Size>
//         + Allocator<f64, P::Size>
//         + Allocator<f64, P::VecLength>
//         + Allocator<f64, P::Size, P::Size>,
// {
//     /// Estimated physical state
//     state: X,
//     /// Estimated paramaters
//     parameters: P,
// }

// impl<X: State, P: State> fmt::Display for NavState<X, P>
// where
//     DefaultAllocator: Allocator<f64, X::Size>
//         + Allocator<f64, X::VecLength>
//         + Allocator<f64, X::Size, X::Size>
//         + Allocator<f64, P::Size>
//         + Allocator<f64, P::VecLength>
//         + Allocator<f64, P::Size, P::Size>,
// {
//     // Prints as Cartesian in floating point with units
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         let decimals = f.precision().unwrap_or(6);
//         write!(
//             f,
//             "state = {}\tparam = {}",
//             format!("{:.*}", decimals, self.state),
//             format!("{:.*}", decimals, self.parameters),
//         )
//     }
// }

// impl<X: State, P: State> fmt::LowerExp for NavState<X, P>
// where
//     DefaultAllocator: Allocator<f64, X::Size>
//         + Allocator<f64, X::VecLength>
//         + Allocator<f64, X::Size, X::Size>
//         + Allocator<f64, P::Size>
//         + Allocator<f64, P::VecLength>
//         + Allocator<f64, P::Size, P::Size>,
// {
//     // Prints as Cartesian in scientific notation with units
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         let decimals = f.precision().unwrap_or(6);
//         write!(
//             f,
//             "state = {}\tparam = {}",
//             format!("{:.*e}", decimals, self.state),
//             format!("{:.*e}", decimals, self.parameters),
//         )
//     }
// }

// impl<X: State, P: State> State for NavState<X, P>
// where
//     X::Size: DimAdd<P::Size>,
//     <X::Size as DimAdd<P::Size>>::Output: DimName,
//     DefaultAllocator: Allocator<f64, X::Size>
//         + Allocator<f64, X::VecLength>
//         + Allocator<f64, X::Size, X::Size>
//         + Allocator<f64, P::Size>
//         + Allocator<f64, P::VecLength>
//         + Allocator<f64, P::Size>
//         + Allocator<f64, P::Size, P::Size>
//         + Allocator<f64, <X::Size as DimAdd<P::Size>>::Output, <X::Size as DimAdd<P::Size>>::Output>
//         + Allocator<f64, <X::Size as DimAdd<P::Size>>::Output>,
// {
//     type Size = DimSum<X::Size, P::Size>;
//     type VecLength = DimSum<X::Size, P::Size>;

//     fn zeros() -> Self {
//         todo!()
//     }

//     fn as_vector(&self) -> Result<OVector<f64, Self::VecLength>, NyxError> {
//         todo!()
//     }

//     fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, NyxError> {
//         todo!()
//     }

//     fn set(
//         &mut self,
//         _: hifitime::Epoch,
//         _: &na::Matrix<
//             f64,
//             Self::VecLength,
//             Const<1_usize>,
//             <DefaultAllocator as Allocator<f64, Self::VecLength>>::Buffer,
//         >,
//     ) -> Result<(), NyxError> {
//         todo!()
//     }

//     fn add(self, other: OVector<f64, Self::Size>) -> Self {
//         todo!()
//     }
//     fn reset_stm(&mut self) {
//         todo!()
//     }
//     fn epoch(&self) -> hifitime::Epoch {
//         todo!()
//     }
//     fn set_epoch(&mut self, _: hifitime::Epoch) {
//         todo!()
//     }
// }

#[test]
fn name() {
    let sc = Spacecraft::zeros();
    let orbit: Orbit = sc.extract();
    println!("{}", orbit);
}

// /// A trait for generate propagation and estimation state.
// /// The first parameter is the size of the state, the second is the size of the propagated state including STM and extra items.
// pub trait State: Copy + PartialEq + fmt::Display + fmt::LowerExp + Send + Sync
// where
//     Self: Sized,
//     Self::Size: DimMul<Self::Size>,
//     DefaultAllocator: Allocator<f64, Self::Size>
//         + Allocator<f64, Self::Size, Self::Size>
//         + Allocator<f64, Self::VecLength>
//         + Allocator<f64, <Self::Size as DimMul<Self::Size>>::Output>,
// {
//     /// Size of the state and its STM
//     type Size: DimName;
//     type VecLength: DimName;
//     /// Initialize an empty state
//     fn zeros() -> Self;

//     fn epoch(&self) -> Epoch;

//     /// Return this state as a vector for the propagation/estimation
//     fn as_vector(&self) -> Result<OVector<f64, DimProd<Self::Size, Self::Size>>, NyxError>;

//     /// Return this state as a vector for the propagation/estimation
//     fn stm(&self) -> Result<OMatrix<f64, Self::Size, Self::Size>, NyxError>;

//     /// Return this state as a vector for the propagation/estimation
//     fn reset_stm(&mut self);

//     /// Set this state
//     fn set(&mut self, epoch: Epoch, vector: &OVector<f64, Self::VecLength>)
//         -> Result<(), NyxError>;

//     /// Reconstruct a new State from the provided delta time in seconds compared to the current state
//     /// and with the provided vector.
//     fn set_with_delta_seconds(
//         self,
//         delta_t_s: f64,
//         vector: &OVector<f64, Self::VecLength>,
//     ) -> Self {
//         let mut me = self;
//         me.set(me.epoch() + delta_t_s, vector).unwrap();
//         me
//     }

//     fn add(self, other: OVector<f64, Self::Size>) -> Self;
// }

// pub trait Embed<S: State> {}
