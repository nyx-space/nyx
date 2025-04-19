/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#![allow(clippy::type_complexity)] // Allow complex types for generics
#![allow(unused_imports)] // Keep imports for context even if slightly unused in snippet

use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName, OMatrix, OVector, U1}; // Use U1 for MsrSize
use crate::md::trajectory::{Interpolatable, Traj}; // May not need Traj if we propagate point-to-point
pub use crate::od::estimate::*;
pub use crate::od::ground_station::*;
pub use crate::od::snc::*; // SNC not typically used in BLS, but keep context
pub use crate::od::*;
use crate::propagators::Propagator;
pub use crate::time::{Duration, Epoch, Unit};
use anise::prelude::Almanac;
use indexmap::IndexSet;
use log::{debug, info, trace, warn};
use msr::sensitivity::TrackerSensitivity; // Assuming this is the correct path
use snafu::prelude::*;
use std::collections::BTreeMap;
use std::fmt;
use std::marker::PhantomData;
use std::ops::Add;
use std::sync::Arc;
use typed_builder::TypedBuilder;

// Define a simple BLS solution struct
#[derive(Debug, Clone)]
pub struct BLSSolution<StateType: State>
where
    DefaultAllocator: Allocator<<StateType as State>::Size>
        + Allocator<<StateType as State>::Size, <StateType as State>::Size>
        + Allocator<<StateType as State>::VecLength>,
{
    pub estimated_state: StateType,
    pub covariance: OMatrix<f64, <StateType as State>::Size, <StateType as State>::Size>,
    pub num_iterations: usize,
    pub final_rms: f64,
    pub final_corr_pos_km: f64,
    pub converged: bool,
}

impl<StateType> fmt::Display for BLSSolution<StateType>
where
    StateType: State,
    DefaultAllocator: Allocator<<StateType as State>::Size>
        + Allocator<<StateType as State>::Size, <StateType as State>::Size>
        + Allocator<<StateType as State>::VecLength>,
{
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        writeln!(f, "Converged: {}", self.converged)?;
        writeln!(f, "Iterations: {}", self.num_iterations)?;
        writeln!(f, "Final RMS: {}", self.final_rms)?;
        writeln!(f, "Final State: {}", self.estimated_state.orbit())?;
        write!(f, "Final Covariance:\n{:.3e}", self.covariance)
    }
}

impl<StateType: State> From<BLSSolution<StateType>> for KfEstimate<StateType>
where
    DefaultAllocator: Allocator<<StateType as State>::Size>
        + Allocator<<StateType as State>::Size, <StateType as State>::Size>
        + Allocator<<StateType as State>::VecLength>,
    <DefaultAllocator as Allocator<<StateType as State>::Size>>::Buffer<f64>: Copy,
    <DefaultAllocator as Allocator<<StateType as State>::Size, <StateType as State>::Size>>::Buffer<
        f64,
    >: Copy,
{
    fn from(mut bls: BLSSolution<StateType>) -> Self {
        // Set the uncertainty on Cr, Cd, mass to zero
        bls.covariance[(6, 6)] = 0.0;
        bls.covariance[(7, 7)] = 0.0;
        bls.covariance[(8, 8)] = 0.0;

        Self::from_covar(bls.estimated_state, bls.covariance)
    }
}
