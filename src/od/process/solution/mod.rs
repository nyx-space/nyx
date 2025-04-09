/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2018-onwards Christopher Rabotin <christopher.rabotin@gmail.com>

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Affero General Public License as published
    by the Free Software Foundation either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Affero General Public License for more details.

    You should have received a copy of the GNU Affero General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

use crate::linalg::allocator::Allocator;
use crate::linalg::{DefaultAllocator, DimName};
use crate::md::trajectory::{Interpolatable, Traj};
pub use crate::od::estimate::*;
pub use crate::od::*;
use indexmap::IndexSet;
use msr::sensitivity::TrackerSensitivity;
use nalgebra::OMatrix;
use std::collections::BTreeMap;
use std::iter::Zip;
use std::ops::Add;
use std::slice::Iter;

use self::msr::MeasurementType;

mod display;
mod export;
mod filter_data;
mod import;
mod smooth;
mod stats;

/// The `ODSolution` structure is designed to manage and analyze the results of an OD process, including
/// smoothing. It provides various functionalities such as splitting solutions by tracker or measurement type,
/// joining solutions, and performing statistical analyses.
///
/// **Note:** Many methods in this structure assume that the solution has been split into subsets using the `split()` method.
/// Calling these methods without first splitting will make analysis of operations results less obvious.
///
/// # Fields
/// - `estimates`: A vector of state estimates generated during the OD process.
/// - `residuals`: A vector of residuals corresponding to the state estimates.
/// - `gains`: Filter gains used for measurement updates. These are set to `None` after running the smoother.
/// - `filter_smoother_ratios`: Filter-smoother consistency ratios. These are set to `None` before running the smoother.
/// - `devices`: A map of tracking devices used in the OD process.
/// - `measurement_types`: A set of unique measurement types used in the OD process.
///
/// Implementation detail: these are not stored in vectors to allow for multiple estimates at the same time, e.g. when
/// there are simultaneous measurements of angles and the filter processes each as a scalar.
///
#[derive(Clone, Debug)]
#[allow(clippy::upper_case_acronyms)]
pub struct ODSolution<StateType, EstType, MsrSize, Trk>
where
    StateType: Interpolatable + Add<OVector<f64, <StateType as State>::Size>, Output = StateType>,
    EstType: Estimate<StateType>,
    MsrSize: DimName,
    Trk: TrackerSensitivity<StateType, StateType>,
    <DefaultAllocator as Allocator<<StateType as State>::VecLength>>::Buffer<f64>: Send,
    DefaultAllocator: Allocator<<StateType as State>::Size>
        + Allocator<<StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <StateType as State>::Size>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<StateType as State>::Size, <StateType as State>::Size>
        + Allocator<<StateType as State>::Size, MsrSize>,
{
    /// Vector of estimates available after a pass
    pub estimates: Vec<EstType>,
    /// Vector of residuals available after a pass
    pub residuals: Vec<Option<Residual<MsrSize>>>,
    /// Vector of filter gains used for each measurement update, all None after running the smoother.
    pub gains: Vec<Option<OMatrix<f64, <StateType as State>::Size, MsrSize>>>,
    /// Filter-smoother consistency ratios, all None before running the smoother.
    pub filter_smoother_ratios: Vec<Option<OVector<f64, <StateType as State>::Size>>>,
    /// Tracking devices
    pub devices: BTreeMap<String, Trk>,
    pub measurement_types: IndexSet<MeasurementType>,
}

impl<StateType, EstType, MsrSize, Trk> ODSolution<StateType, EstType, MsrSize, Trk>
where
    StateType: Interpolatable + Add<OVector<f64, <StateType as State>::Size>, Output = StateType>,
    EstType: Estimate<StateType>,
    MsrSize: DimName,
    Trk: TrackerSensitivity<StateType, StateType>,
    <DefaultAllocator as Allocator<<StateType as State>::VecLength>>::Buffer<f64>: Send,
    DefaultAllocator: Allocator<<StateType as State>::Size>
        + Allocator<<StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <StateType as State>::Size>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<StateType as State>::Size, <StateType as State>::Size>
        + Allocator<<StateType as State>::Size, MsrSize>,
{
    pub fn new(
        devices: BTreeMap<String, Trk>,
        measurement_types: IndexSet<MeasurementType>,
    ) -> Self {
        Self {
            estimates: Vec::new(),
            residuals: Vec::new(),
            gains: Vec::new(),
            filter_smoother_ratios: Vec::new(),
            devices,
            measurement_types,
        }
    }

    /// Pushes a new measurement update result, ensuring proper sizes of the arrays.
    pub(crate) fn push_measurement_update(
        &mut self,
        estimate: EstType,
        residual: Residual<MsrSize>,
        gain: Option<OMatrix<f64, <StateType as State>::Size, MsrSize>>,
    ) {
        self.estimates.push(estimate);
        self.residuals.push(Some(residual));
        self.gains.push(gain);
        self.filter_smoother_ratios.push(None);
    }

    /// Pushes a new time update result, ensuring proper sizes of the arrays.
    pub(crate) fn push_time_update(&mut self, estimate: EstType) {
        self.estimates.push(estimate);
        self.residuals.push(None);
        self.gains.push(None);
        self.filter_smoother_ratios.push(None);
    }

    /// Returns a zipper iterator on the estimates and the associated residuals.
    pub fn results(&self) -> Zip<Iter<'_, EstType>, Iter<'_, Option<Residual<MsrSize>>>> {
        self.estimates.iter().zip(self.residuals.iter())
    }

    /// Returns True if this is the result of a filter run
    pub fn is_filter_run(&self) -> bool {
        self.gains.iter().flatten().count() > 0
    }

    /// Returns True if this is the result of a smoother run
    pub fn is_smoother_run(&self) -> bool {
        self.filter_smoother_ratios.iter().flatten().count() > 0
    }

    /// Builds the navigation trajectory for the estimated state only
    pub fn to_traj(&self) -> Result<Traj<StateType>, NyxError>
    where
        DefaultAllocator: Allocator<StateType::VecLength>,
    {
        if self.estimates.is_empty() {
            Err(NyxError::NoStateData {
                msg: "No navigation trajectory to generate: run the OD process first".to_string(),
            })
        } else {
            // Make sure to remove duplicate entries.
            let mut traj = Traj {
                states: self.estimates.iter().map(|est| est.state()).collect(),
                name: None,
            };
            traj.finalize();
            Ok(traj)
        }
    }

    /// Returns the accepted residuals.
    pub fn accepted_residuals(&self) -> Vec<Residual<MsrSize>> {
        self.residuals
            .iter()
            .flatten()
            .filter(|resid| !resid.rejected)
            .cloned()
            .collect::<Vec<Residual<MsrSize>>>()
    }

    /// Returns the rejected residuals.
    pub fn rejected_residuals(&self) -> Vec<Residual<MsrSize>> {
        self.residuals
            .iter()
            .flatten()
            .filter(|resid| resid.rejected)
            .cloned()
            .collect::<Vec<Residual<MsrSize>>>()
    }
}

impl<StateType, EstType, MsrSize, Trk> PartialEq for ODSolution<StateType, EstType, MsrSize, Trk>
where
    StateType: Interpolatable + Add<OVector<f64, <StateType as State>::Size>, Output = StateType>,
    EstType: Estimate<StateType>,
    MsrSize: DimName,
    Trk: TrackerSensitivity<StateType, StateType> + PartialEq,
    <DefaultAllocator as Allocator<<StateType as State>::VecLength>>::Buffer<f64>: Send,
    DefaultAllocator: Allocator<<StateType as State>::Size>
        + Allocator<<StateType as State>::VecLength>
        + Allocator<MsrSize>
        + Allocator<MsrSize, <StateType as State>::Size>
        + Allocator<MsrSize, MsrSize>
        + Allocator<<StateType as State>::Size, <StateType as State>::Size>
        + Allocator<<StateType as State>::Size, MsrSize>,
{
    fn eq(&self, other: &Self) -> bool {
        self.estimates.len() == other.estimates.len()
            && self.residuals == other.residuals
            && self.gains.len() == other.gains.len()
            && self.filter_smoother_ratios.len() == other.filter_smoother_ratios.len()
            && self.devices == other.devices
            && self.measurement_types == other.measurement_types
            // Now check for near equality of gains
            && self.gains.iter().zip(other.gains.iter()).all(|(my_k, other_k)| {
                if let Some(my_k) = my_k {
                    if let Some(other_k) = other_k {
                        (my_k - other_k).norm() < 1e-9
                    } else {
                        false
                    }
                } else {
                    other_k.is_none()
                }
            })
            // Now check for near equality of F-S ratios
            && self.filter_smoother_ratios.iter().zip(other.filter_smoother_ratios.iter()).all(|(my_fs, other_fs)| {
                if let Some(my_fs) = my_fs {
                    if let Some(other_fs) = other_fs {
                        (my_fs - other_fs).norm() < 1e-9
                    } else {
                        false
                    }
                } else {
                    other_fs.is_none()
                }
            })
    }
}
