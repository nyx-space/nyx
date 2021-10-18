use super::parameters::NoParameters;
use super::state::{EmbedState, NavState};
use super::SeedableRng;
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

/// A NavModel returns a full sensor reading of all sensors at the requested time.
/// For example, for a ground station that support simultaneous two way doppler (measurement/obs of size 1) and angles observations (obs of size 2)
/// in order to estimate a navigation state composed of an Orbit and a ClockBias then an implementation of NavModel would be a structure
/// that does: `impl NavModel<Const<3>, NavState<Orbit, ClockBias>, Orbit, ClockBias> {}`
pub trait NavModel<N, S, X, P>
where
    N: DimName,
    S: NavState<X, P> + EmbedState<P>,
    P: State,
    X: State,
    DefaultAllocator: Allocator<f64, X::Size>
        + Allocator<f64, X::VecLength>
        + Allocator<f64, X::Size, X::Size>
        + Allocator<f64, P::Size>
        + Allocator<f64, P::VecLength>
        + Allocator<f64, P::Size, P::Size>
        + Allocator<f64, S::Size>
        + Allocator<f64, S::VecLength>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<f64, N>
        + Allocator<f64, N, S::Size>,
{
    /// Returns the latest observation from the sensor and the epoch of that observation
    fn latest_observation(&self) -> (Epoch, OVector<f64, N>);
    /// Returns the relevant observation for the requested Epoch.
    /// This method allows supporting lag in observations. There are two main strategies to handle this based on the Epoch (embedded in the NavState).
    /// **First**, one may simply wish to null the observation (and the relevant computed observation!) if the time difference
    /// between the requested epoch and the latest measurement is too big (thereby nulling any contribution of that sensor to the nav state).
    /// **Second** strategy consists in propagating the measurement forward using a measurement model. If so, this needs to be handled by the implementor of NavSensor.
    fn relevant_observation(&self, current_state: S) -> Result<OVector<f64, N>, NyxError>;
    /// Returns the computed observation and its associated sensitivity matrix for the full state S, for the current state (which also contains an Epoch).
    fn computed_observation(
        &self,
        current_state: S,
    ) -> Result<(OVector<f64, N>, OMatrix<f64, N, S::Size>), NyxError>;
}

/*
// TwoWayDoppler + OneWayAngles + ClockBias ground station
struct GroundStation {
    twowaydoppler: TwoWayDoppler, // Size 1
    onewayangles: OneWayAngles, // Size 2
    clockbias: ClockBias // Size 2
}

pub struct ODEstState = NavState<Orbit, ClockBias>

impl NavModel<Const<5>, ODEstState, Orbit, ClockBias> {}
*/

/// A SimulatedNavModel allows the simulation of a navigation model
/// For example, for a ground station that support simultaneous two way doppler (measurement/obs of size 1) and angles observations (obs of size 2)
/// in order to estimate a navigation state composed of an Orbit and a ClockBias then an implementation of NavModel would be a structure
/// that does: `impl SimulatedNavModel<Const<3>, NavState<Orbit, ClockBias>, Orbit, ClockBias> {}`
pub trait SimulatedNavModel<N, S, X, P, R>
where
    N: DimName,
    S: NavState<X, P> + EmbedState<P>,
    P: State,
    X: State,
    R: SeedableRng,
    DefaultAllocator: Allocator<f64, X::Size>
        + Allocator<f64, X::VecLength>
        + Allocator<f64, X::Size, X::Size>
        + Allocator<f64, P::Size>
        + Allocator<f64, P::VecLength>
        + Allocator<f64, P::Size, P::Size>
        + Allocator<f64, S::Size>
        + Allocator<f64, S::VecLength>
        + Allocator<f64, S::Size, S::Size>
        + Allocator<f64, N>
        + Allocator<f64, N, S::Size>,
{
    /// Generate a simulated observation from a current navigation state and a seedable random number generator
    fn generate_observation(
        &self,
        current_state: S,
        rng: &mut R,
    ) -> Result<OVector<f64, N>, NyxError>;
}
