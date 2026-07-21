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
use super::{MeasurementType, measurement::Measurement};
use core::fmt;
use hifitime::prelude::{Duration, Epoch};
use indexmap::{IndexMap, IndexSet};
use log::{info, warn};
use std::ops::Bound::{self, Excluded, Included, Unbounded};
use std::ops::{Add, AddAssign, RangeBounds};

mod io_ccsds_tdm;
mod io_parquet;

#[cfg(feature = "python")]
use pyo3::prelude::*;
#[cfg(feature = "python")]
mod python;

/// Tracking data storing all of measurements as a B-Tree.
/// It inherently does NOT support multiple concurrent measurements from several trackers.
///
/// # Measurement Moduli, e.g. range modulus
///
/// In the case of ranging, and possibly other data types, a code is used to measure the range to the spacecraft. The length of this code
/// determines the ambiguity resolution, as per equation 9 in section 2.2.2.2 of the JPL DESCANSO, document 214, _Pseudo-Noise and Regenerative Ranging_.
/// For example, using the JPL Range Code and a frequency range clock of 1 MHz, the range ambiguity is 75,660 km. In other words,
/// as soon as the spacecraft is at a range of 75,660 + 1 km the JPL Range Code will report the vehicle to be at a range of 1 km.
/// This is simply because the range code overlaps with itself, effectively loosing track of its own reference:
/// it's due to the phase shift of the signal "lapping" the original signal length.
///
/// ```text
///             (Spacecraft)
///             ^
///             |    Actual Distance = 75,661 km
///             |
/// 0 km                                         75,660 km (Wrap-Around)
/// |-----------------------------------------------|
///   When the "code length" is exceeded,
///   measurements wrap back to 0.
///
/// So effectively:
///     Observed code range = Actual range (mod 75,660 km)
///     75,661 km → 1 km
///
/// ```
///
/// Nyx can only resolve the range ambiguity if the tracking data specifies a modulus for this specific measurement type.
/// For example, in the case of the JPL Range Code and a 1 MHz range clock, the ambiguity interval is 75,660 km.
///
/// The measurement used in the Orbit Determination Process then becomes the following, where `//` represents the [Euclidian division](https://doc.rust-lang.org/std/primitive.f64.html#method.div_euclid).
///
/// ```text
/// k = computed_obs // ambiguity_interval
/// real_obs = measured_obs + k * modulus
/// ```
///
/// Reference: JPL DESCANSO, document 214, _Pseudo-Noise and Regenerative Ranging_.
///
#[derive(Clone, Default)]
#[cfg_attr(feature = "python", pyclass(from_py_object))]
pub struct TrackingDataArc {
    /// All measurements in this data arc
    pub measurements: Vec<Measurement>,
    /// Source file if loaded from a file or saved to a file.
    pub source: Option<String>,
    /// Optionally provide a map of modulos (e.g. the RANGE_MODULO of CCSDS TDM).
    pub moduli: Option<IndexMap<MeasurementType, f64>>,
    /// Reject all of the measurements, useful for debugging passes.
    pub force_reject: bool,
}

#[cfg_attr(feature = "python", pymethods)]
impl TrackingDataArc {
    /// Sort these measurements by epoch
    pub fn sort(&mut self) {
        self.measurements.sort_unstable_by(|a, b| {
            a.epoch
                .cmp(&b.epoch)
                .then_with(|| a.tracker.cmp(&b.tracker))
        });

        // Coalesce adjacent duplicate elements in exactly O(K) time.
        // dedup_by passes pointers to `(next_element, kept_element)`.
        // If the closure returns true, `next_element` is physically dropped.
        self.measurements.dedup_by(|next, kept| {
            if next.epoch == kept.epoch && next.tracker == kept.tracker {
                // The tracker and epoch are identical. Drain the sub-observables
                // from the redundant 'next' measurement and merge them into the 'kept' one.
                kept.data.extend(next.data.drain(..));

                // If either partial record was manually flagged as rejected,
                // the combined radiometric record must retain that suspicion.
                kept.rejected |= next.rejected;

                // Return true to destroy the redundant parent struct.
                true
            } else {
                // Elements differ structurally. Keep both.
                false
            }
        });
    }
    /// Returns the start epoch of this tracking arc
    pub fn start_epoch(&self) -> Option<Epoch> {
        self.measurements.first().map(|msr| msr.epoch)
    }

    /// Returns the end epoch of this tracking arc
    pub fn end_epoch(&self) -> Option<Epoch> {
        self.measurements.last().map(|msr| msr.epoch)
    }

    /// Returns the duration this tracking arc
    pub fn duration(&self) -> Option<Duration> {
        match self.start_epoch() {
            Some(start) => self.end_epoch().map(|end| end - start),
            None => None,
        }
    }

    /// Returns the number of measurements in this data arc
    pub fn len(&self) -> usize {
        self.measurements.len()
    }

    /// Returns whether this arc has no measurements.
    pub fn is_empty(&self) -> bool {
        self.measurements.is_empty()
    }

    /// Returns the minimum duration between two subsequent measurements.
    pub fn min_duration_sep(&self) -> Option<Duration> {
        if self.is_empty() {
            None
        } else {
            let mut min_sep = Duration::MAX;
            let mut prev_epoch = self.start_epoch().unwrap();
            for msr in self.measurements.iter().skip(1) {
                let epoch = msr.epoch;
                let this_sep = epoch - prev_epoch;
                min_sep = min_sep.min(this_sep);
                prev_epoch = epoch;
            }
            Some(min_sep)
        }
    }
    /// Set (or overwrites) the modulus of the provided measurement type.
    pub fn set_moduli(&mut self, msr_type: MeasurementType, modulus: f64) {
        if modulus.is_nan() || modulus.abs() < f64::EPSILON {
            warn!("cannot set modulus for {msr_type:?} to {modulus}");
            return;
        }
        if self.moduli.is_none() {
            self.moduli = Some(IndexMap::new());
        }

        self.moduli.as_mut().unwrap().insert(msr_type, modulus);
    }

    /// Applies the moduli to each measurement, if defined.
    pub fn apply_moduli(&mut self) {
        if let Some(moduli) = &self.moduli {
            for msr in &mut self.measurements {
                for (msr_type, modulus) in moduli {
                    if let Some(msr_value) = msr.data.get_mut(msr_type) {
                        *msr_value %= *modulus;
                    }
                }
            }
        }
    }

    /// Downsamples the tracking data to a lower frequency using a simple moving average low-pass filter followed by decimation,
    /// returning new `TrackingDataArc` with downsampled measurements.
    ///
    /// It provides a computationally efficient approach to reduce the sampling rate while mitigating aliasing effects.
    ///
    /// # Algorithm
    ///
    /// 1. A simple moving average filter is applied as a low-pass filter.
    /// 2. Decimation is performed by selecting every Nth sample after filtering.
    ///
    /// # Advantages
    ///
    /// - Computationally efficient, suitable for large datasets common in spaceflight applications.
    /// - Provides basic anti-aliasing, crucial for preserving signal integrity in orbit determination and tracking.
    /// - Maintains phase information, important for accurate timing in spacecraft state estimation.
    ///
    /// # Limitations
    ///
    /// - The frequency response is not as sharp as more sophisticated filters (e.g., FIR, IIR).
    /// - May not provide optimal stopband attenuation for high-precision applications.
    ///
    /// ## Considerations for Spaceflight Applications
    ///
    /// - Suitable for initial data reduction in ground station tracking pipelines.
    /// - Adequate for many orbit determination and tracking tasks where computational speed is prioritized.
    /// - For high-precision applications (e.g., interplanetary navigation), consider using more advanced filtering techniques.
    ///
    /// :type target_step: Duration
    /// :rtype: TrackingDataArc
    pub fn downsample(&self, target_step: Duration) -> Self {
        if self.is_empty() {
            return self.clone();
        }
        let current_step = self.min_duration_sep().unwrap();

        if current_step >= target_step {
            warn!(
                "cannot downsample tracking data from {current_step} to {target_step} (that would be upsampling)"
            );
            return self.clone();
        }

        let current_hz = 1.0 / current_step.to_seconds();
        let target_hz = 1.0 / target_step.to_seconds();

        // Simple moving average as low-pass filter
        let window_size = (current_hz / target_hz).round() as usize;

        info!(
            "downsampling tracking data from {current_step} ({current_hz:.6} Hz) to {target_step} ({target_hz:.6} Hz) (N = {window_size})"
        );

        let mut result = TrackingDataArc {
            source: self.source.clone(),
            ..Default::default()
        };

        let measurements: Vec<_> = self.measurements.iter().collect();

        for (i, msr) in measurements.iter().enumerate().step_by(window_size) {
            let epoch = msr.epoch;
            let start = i.saturating_sub(window_size / 2);
            let end = (i + window_size / 2 + 1).min(measurements.len());
            let window = &measurements[start..end];

            let mut filtered_measurement = Measurement {
                tracker: window[0].tracker.clone(),
                epoch,
                data: IndexMap::new(),
                rejected: false,
            };

            // Apply moving average filter for each measurement type
            for mtype in self.unique_types() {
                let sum: f64 = window.iter().filter_map(|m| m.data.get(&mtype)).sum();
                let count = window
                    .iter()
                    .filter(|m| m.data.contains_key(&mtype))
                    .count();

                if count > 0 {
                    filtered_measurement.data.insert(mtype, sum / count as f64);
                }
            }

            result.measurements.push(filtered_measurement);
        }
        result.sort();
        result
    }

    /// Splits a long tracking data arc into smaller chunks, each up to `max_duration` long.
    pub fn chunk(&self, max_duration: Duration) -> Vec<TrackingDataArc> {
        let mut chunks = Vec::new();
        if self.is_empty() || max_duration <= Duration::ZERO {
            return chunks;
        }

        let mut start_idx = 0;
        let total_measurements = self.measurements.len();

        while start_idx < total_measurements {
            let chunk_start_epoch = self.measurements[start_idx].epoch;
            let chunk_end_time = chunk_start_epoch + max_duration;

            // Isolate the remaining, unprocessed portion of the vector
            let remaining = &self.measurements[start_idx..];

            // Perform a binary search on the remaining slice to find the first
            // index that strictly exceeds the chunk_end_time.
            let offset = remaining.partition_point(|msr| msr.epoch <= chunk_end_time);

            let end_idx = start_idx + offset;

            // Extract and clone ONLY the measurements belonging to this chunk.
            // This drops the memory complexity from O(K * N) to strictly O(N).
            let chunk_measurements = self.measurements[start_idx..end_idx].to_vec();

            chunks.push(TrackingDataArc {
                measurements: chunk_measurements,
                source: self.source.clone(),
                moduli: self.moduli.clone(),
                force_reject: self.force_reject,
            });

            // Advance the window to the exact start of the next chunk
            start_idx = end_idx;
        }

        chunks
    }
}

impl TrackingDataArc {
    /// Helper method to resolve bounds into slice indices via binary search.
    fn resolve_bounds<R: RangeBounds<Epoch>>(&self, bound: R) -> (usize, usize) {
        // Find the lower bound index via O(log N) binary search
        let start_idx = match bound.start_bound() {
            Bound::Included(&epoch) => self.measurements.partition_point(|m| m.epoch < epoch),
            Bound::Excluded(&epoch) => self.measurements.partition_point(|m| m.epoch <= epoch),
            Bound::Unbounded => 0,
        };

        // Find the upper bound index via O(log N) binary search
        let end_idx = match bound.end_bound() {
            Bound::Included(&epoch) => self.measurements.partition_point(|m| m.epoch <= epoch),
            Bound::Excluded(&epoch) => self.measurements.partition_point(|m| m.epoch < epoch),
            Bound::Unbounded => self.measurements.len(),
        };

        (start_idx, end_idx)
    }

    /// Returns the unique list of aliases in this tracking data arc
    pub fn unique_aliases(&self) -> IndexSet<String> {
        self.unique().0
    }

    /// Returns the unique measurement types in this tracking data arc
    pub fn unique_types(&self) -> IndexSet<MeasurementType> {
        self.unique().1
    }

    /// Returns the unique trackers and unique measurement types in this data arc
    pub fn unique(&self) -> (IndexSet<String>, IndexSet<MeasurementType>) {
        let mut aliases = IndexSet::new();
        let mut types = IndexSet::new();
        for msr in &self.measurements {
            aliases.insert(msr.tracker.clone());
            for k in msr.data.keys() {
                types.insert(*k);
            }
        }
        (aliases, types)
    }

    /// Returns a new tracking arc that only contains measurements that fall within the given epoch range.
    ///
    /// Executes in O(N) time strictly due to memory shifting, requiring zero new allocations.
    pub fn filter_by_epoch<R: RangeBounds<Epoch>>(mut self, bound: R) -> Self {
        let (start_idx, end_idx) = self.resolve_bounds(bound);

        // Handle disjoint bounds or out-of-range queries
        if start_idx >= end_idx || start_idx >= self.measurements.len() {
            self.measurements.clear();
            return self;
        }

        // In-place memory reduction
        // Truncate the tail first. This drops trailing measurements without shifting.
        self.measurements.truncate(end_idx);

        // Drain the head. This removes preceding measurements and shifts the
        // remaining valid data leftward to index 0 in a single memory move.
        self.measurements.drain(0..start_idx);

        // Note that the order is preserved, so we don't need to sort again.

        // Clear unused memory
        self.measurements.shrink_to_fit();

        self
    }

    /// Returns a new tracking arc that only contains measurements that fall within the given offset from the first epoch.
    /// For example, a bound of 30.minutes()..90.minutes() will only read measurements from the start of the arc + 30 minutes until start + 90 minutes.
    pub fn filter_by_offset<R: RangeBounds<Duration>>(self, bound: R) -> Self {
        if self.is_empty() {
            return self;
        }
        // Rebuild an epoch bound.
        let start = match bound.start_bound() {
            Unbounded => self.start_epoch().unwrap(),
            Included(offset) | Excluded(offset) => self.start_epoch().unwrap() + *offset,
        };

        let end = match bound.end_bound() {
            Unbounded => self.end_epoch().unwrap(),
            Included(offset) | Excluded(offset) => self.start_epoch().unwrap() + *offset,
        };

        self.filter_by_epoch(start..end)
    }

    /// Returns a new tracking arc that only contains measurements from the desired tracker.
    pub fn filter_by_tracker(mut self, tracker: String) -> Self {
        self.measurements = self
            .measurements
            .iter()
            .filter_map(|msr| {
                if msr.tracker == tracker {
                    Some(msr.clone())
                } else {
                    None
                }
            })
            .collect::<Vec<Measurement>>();
        self
    }

    /// Returns a new tracking arc that only contains measurements of the provided type.
    pub fn filter_by_measurement_type(mut self, included_type: MeasurementType) -> Self {
        self.measurements.retain_mut(|msr| {
            msr.data.retain(|msr_type, _| *msr_type == included_type);
            !msr.data.is_empty()
        });
        self
    }

    /// Returns a new tracking arc that contains measurements from all trackers except the one provided
    pub fn exclude_tracker(mut self, excluded_tracker: String) -> Self {
        self.measurements = self
            .measurements
            .iter()
            .filter_map(|msr| {
                if msr.tracker != excluded_tracker {
                    Some(msr.clone())
                } else {
                    None
                }
            })
            .collect::<Vec<Measurement>>();
        self
    }

    /// Returns a new tracking arc that excludes measurements within the given epoch range.
    ///
    /// Executes an in-place O(N) memory shift with zero heap allocations.
    pub fn exclude_by_epoch<R: RangeBounds<Epoch>>(mut self, bound: R) -> Self {
        let (start_idx, end_idx) = self.resolve_bounds(bound);

        if start_idx < end_idx && start_idx < self.measurements.len() {
            // Drain removes the specified range and shifts all subsequent elements
            // leftward to fill the gap. The extracted elements are immediately dropped.
            self.measurements.drain(start_idx..end_idx);
        }

        self
    }

    /// Returns a new tracking arc that contains measurements from all trackers except the one provided
    pub fn exclude_measurement_type(mut self, excluded_type: MeasurementType) -> Self {
        self.measurements = self
            .measurements
            .iter_mut()
            .map(|msr| {
                msr.data.retain(|msr_type, _| *msr_type != excluded_type);
                msr.clone()
            })
            .collect::<Vec<Measurement>>();
        self
    }

    /// Marks measurements within the given epoch range as rejected.
    ///
    /// Operates in O(log N) for bound resolution and O(K) for iteration, where K is the slice length.
    pub fn reject_by_epoch<R: RangeBounds<Epoch>>(mut self, bound: R) -> Self {
        let (start_idx, end_idx) = self.resolve_bounds(bound);

        if start_idx < end_idx && start_idx < self.measurements.len() {
            for msr in &mut self.measurements[start_idx..end_idx] {
                msr.rejected = true;
            }
        }
        self
    }

    /// Marks measurements from the provided tracker as rejected.
    /// Requires an O(N) scan. The parameter is downgraded to &str to prevent heap allocations.
    pub fn reject_by_tracker(mut self, tracker: &str) -> Self {
        for msr in &mut self.measurements {
            if msr.tracker == tracker {
                msr.rejected = true;
            }
        }
        self
    }

    pub fn resid_vs_ref_check(mut self) -> Self {
        self.force_reject = true;
        self
    }
}

impl fmt::Display for TrackingDataArc {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.is_empty() {
            write!(f, "Empty tracking arc")
        } else {
            let start = self.start_epoch().unwrap();
            let end = self.end_epoch().unwrap();
            let src = match &self.source {
                Some(src) => format!(" (source: {src})"),
                None => String::new(),
            };
            write!(
                f,
                "Tracking arc with {} measurements of type {:?} over {} (from {start} to {end}) with trackers {:?}{src}",
                self.len(),
                self.unique_types(),
                end - start,
                self.unique_aliases()
            )
        }
    }
}

impl fmt::Debug for TrackingDataArc {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{self} @ {self:p}")
    }
}

impl PartialEq for TrackingDataArc {
    fn eq(&self, other: &Self) -> bool {
        self.measurements == other.measurements
    }
}

impl Add for TrackingDataArc {
    type Output = Self;

    fn add(mut self, rhs: Self) -> Self::Output {
        self.force_reject = false;
        self.measurements.extend(rhs.measurements);
        self.sort();

        self.force_reject = false;
        self
    }
}

impl AddAssign for TrackingDataArc {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}
