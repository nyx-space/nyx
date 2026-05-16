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
use crate::linalg::{DefaultAllocator, DimName, OMatrix};
use crate::md::trajectory::Interpolatable;
pub use crate::od::estimate::*;
pub use crate::od::*;
use indexmap::IndexSet;
use msr::sensitivity::TrackerSensitivity;
use std::ops::Add;

use self::msr::MeasurementType;

use super::ODSolution;

#[derive(Clone, Debug, PartialEq)]
#[allow(clippy::upper_case_acronyms)]
pub struct ODRecord<StateType, EstType, MsrSize>
where
    StateType: Interpolatable + Add<OVector<f64, <StateType as State>::Size>, Output = StateType>,
    EstType: Estimate<StateType>,
    MsrSize: DimName,
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
    pub estimate: EstType,
    /// Vector of residuals available after a pass
    pub residual: Option<Residual<MsrSize>>,
    /// Vector of filter gains used for each measurement update, all None after running the smoother.
    pub gain: Option<OMatrix<f64, <StateType as State>::Size, MsrSize>>,
    /// Filter-smoother consistency ratios, all None before running the smoother.
    pub filter_smoother_ratio: Option<OVector<f64, <StateType as State>::Size>>,
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
    /// Returns a set of tuples of tracker and measurement types in this OD solution, e.g. `{(Canberra, Range), (Canberra, Doppler)}`.
    pub fn unique(&self) -> IndexSet<(String, MeasurementType)> {
        let mut mapping = IndexSet::new();
        for resid in self.residuals.iter().flatten() {
            for k in &resid.msr_types {
                mapping.insert((
                    resid.tracker.clone().expect("residual tracker not set!"),
                    *k,
                ));
            }
        }
        mapping
    }

    /// Returns this OD solution without any time update
    pub fn drop_time_updates(mut self) -> Self {
        let mut estimates = Vec::new();
        let mut residuals = Vec::new();
        let mut gains = Vec::new();
        let mut filter_smoother_ratios = Vec::new();

        for (est, (resid_opt, (gain_opt, fsr_opt))) in self.estimates.iter().zip(
            self.residuals
                .iter()
                .zip(self.gains.iter().zip(self.filter_smoother_ratios.iter())),
        ) {
            if resid_opt.is_none() {
                continue;
            }

            estimates.push(est.clone());
            residuals.push(resid_opt.clone());
            gains.push(gain_opt.clone());
            filter_smoother_ratios.push(fsr_opt.clone());
        }

        self.estimates = estimates;
        self.residuals = residuals;
        self.gains = gains;
        self.filter_smoother_ratios = filter_smoother_ratios;

        self
    }

    /// Returns this OD solution with only data from the desired measurement type, dropping all time updates.
    pub fn filter_by_msr_type(mut self, msr_type: MeasurementType) -> Self {
        let mut estimates = Vec::new();
        let mut residuals = Vec::new();
        let mut gains = Vec::new();
        let mut filter_smoother_ratios = Vec::new();

        for (est, (resid_opt, (gain_opt, fsr_opt))) in self.estimates.iter().zip(
            self.residuals
                .iter()
                .zip(self.gains.iter().zip(self.filter_smoother_ratios.iter())),
        ) {
            match resid_opt {
                None => continue, // Drop all time updates
                Some(resid) => {
                    if resid.msr_types.contains(&msr_type) {
                        estimates.push(est.clone());
                        residuals.push(Some(resid.clone()));
                        gains.push(gain_opt.clone());
                        filter_smoother_ratios.push(fsr_opt.clone());
                    }
                }
            }
        }

        self.estimates = estimates;
        self.residuals = residuals;
        self.gains = gains;
        self.filter_smoother_ratios = filter_smoother_ratios;

        self
    }

    /// Returns this OD solution with only data from the desired tracker, dropping all time updates.
    pub fn filter_by_tracker(mut self, tracker: String) -> Self {
        let mut estimates = Vec::new();
        let mut residuals = Vec::new();
        let mut gains = Vec::new();
        let mut filter_smoother_ratios = Vec::new();

        for (est, (resid_opt, (gain_opt, fsr_opt))) in self.estimates.iter().zip(
            self.residuals
                .iter()
                .zip(self.gains.iter().zip(self.filter_smoother_ratios.iter())),
        ) {
            match resid_opt {
                None => continue, // Drop all time updates
                Some(resid) => {
                    if resid.tracker == Some(tracker.clone()) {
                        estimates.push(est.clone());
                        residuals.push(Some(resid.clone()));
                        gains.push(gain_opt.clone());
                        filter_smoother_ratios.push(fsr_opt.clone());
                    }
                }
            }
        }

        self.estimates = estimates;
        self.residuals = residuals;
        self.gains = gains;
        self.filter_smoother_ratios = filter_smoother_ratios;

        self
    }

    /// Returns this OD solution with all data except from the desired tracker, including all time updates
    pub fn exclude_tracker(mut self, excluded_tracker: String) -> Self {
        let mut estimates = Vec::new();
        let mut residuals = Vec::new();
        let mut gains = Vec::new();
        let mut filter_smoother_ratios = Vec::new();

        for (est, (resid_opt, (gain_opt, fsr_opt))) in self.estimates.iter().zip(
            self.residuals
                .iter()
                .zip(self.gains.iter().zip(self.filter_smoother_ratios.iter())),
        ) {
            if let Some(resid) = resid_opt {
                if resid.tracker.is_none() || resid.tracker.as_ref().unwrap() == &excluded_tracker {
                    continue;
                }
            }
            // Otherwise, include in the result.
            estimates.push(est.clone());
            residuals.push(resid_opt.clone());
            gains.push(gain_opt.clone());
            filter_smoother_ratios.push(fsr_opt.clone());
        }

        self.estimates = estimates;
        self.residuals = residuals;
        self.gains = gains;
        self.filter_smoother_ratios = filter_smoother_ratios;

        self
    }

    /// Split this OD solution per tracker and per measurement type, dropping all time updates.
    pub fn split(self) -> Vec<Self> {
        let uniques = self.unique();
        let mut splt = Vec::with_capacity(uniques.len());
        for (tracker, msr_type) in uniques {
            splt.push(
                self.clone()
                    .filter_by_tracker(tracker)
                    .filter_by_msr_type(msr_type),
            );
        }
        splt
    }

    /// Merge this OD solution with another one, returning a new OD solution.
    pub fn merge(mut self, mut other: Self) -> Self {
        self.estimates.append(&mut other.estimates);
        self.residuals.append(&mut other.residuals);
        self.gains.append(&mut other.gains);
        self.filter_smoother_ratios
            .append(&mut other.filter_smoother_ratios);

        // Sort to ensure chronological order using indices based permutations.
        // Generate indices representing original positions
        let mut indices: Vec<usize> = (0..self.estimates.len()).collect();

        // Sort indices based on estimates' epochs
        indices.sort_by(|&a, &b| {
            self.estimates[a]
                .epoch()
                .partial_cmp(&self.estimates[b].epoch())
                .expect("Epoch comparison failed")
        });

        // Apply permutation to both vectors in-place
        for i in 0..indices.len() {
            let current = i;
            while indices[current] != current {
                let target = indices[current];
                indices.swap(current, target);
                self.estimates.swap(current, target);
                self.residuals.swap(current, target);
                self.gains.swap(current, target);
                self.filter_smoother_ratios.swap(current, target);
            }
        }

        self
    }

    pub fn at(&self, epoch: Epoch) -> Option<ODRecord<StateType, EstType, MsrSize>> {
        if let Ok(index) = self
            .estimates
            .binary_search_by(|est| est.epoch().cmp(&epoch))
        {
            Some(ODRecord {
                estimate: self.estimates[index].clone(),
                residual: self.residuals[index].clone(),
                gain: self.gains[index].clone(),
                filter_smoother_ratio: self.filter_smoother_ratios[index].clone(),
            })
        } else {
            None
        }
    }
}
