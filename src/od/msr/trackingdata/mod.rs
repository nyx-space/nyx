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
use super::{measurement::Measurement, MeasurementType};
use core::fmt;
use hifitime::prelude::{Duration, Epoch};
use indexmap::{IndexMap, IndexSet};
use std::collections::BTreeMap;
use std::ops::Bound::{Excluded, Included, Unbounded};
use std::ops::RangeBounds;

mod io_ccsds_tdm;
mod io_parquet;

/// Tracking data storing all of measurements as a B-Tree.
/// It inherently does NOT support multiple concurrent measurements from several trackers.
#[derive(Clone, Default)]
pub struct TrackingDataArc {
    /// All measurements in this data arc
    pub measurements: BTreeMap<Epoch, Measurement>, // BUG: Consider a map of tracking to epoch!
    /// Source file if loaded from a file or saved to a file.
    pub source: Option<String>,
}

impl TrackingDataArc {
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
        for msr in self.measurements.values() {
            aliases.insert(msr.tracker.clone());
            for k in msr.data.keys() {
                types.insert(*k);
            }
        }
        (aliases, types)
    }

    /// Returns the start epoch of this tracking arc
    pub fn start_epoch(&self) -> Option<Epoch> {
        self.measurements.first_key_value().map(|(k, _)| *k)
    }

    /// Returns the end epoch of this tracking arc
    pub fn end_epoch(&self) -> Option<Epoch> {
        self.measurements.last_key_value().map(|(k, _)| *k)
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
    /// This is important to correctly set up the propagator and not miss any measurement.
    pub fn min_duration_sep(&self) -> Option<Duration> {
        if self.is_empty() {
            None
        } else {
            let mut min_sep = Duration::MAX;
            let mut prev_epoch = self.start_epoch().unwrap();
            for (epoch, _) in self.measurements.iter().skip(1) {
                let this_sep = *epoch - prev_epoch;
                min_sep = min_sep.min(this_sep);
                prev_epoch = *epoch;
            }
            Some(min_sep)
        }
    }

    /// Returns a new tracking arc that only contains measurements that fall within the given epoch range.
    pub fn filter_by_epoch<R: RangeBounds<Epoch>>(mut self, bound: R) -> Self {
        self.measurements = self
            .measurements
            .range(bound)
            .map(|(epoch, msr)| (*epoch, msr.clone()))
            .collect::<BTreeMap<Epoch, Measurement>>();
        self
    }

    /// Returns a new tracking arc that only contains measurements that fall within the given offset from the first epoch
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
            Included(offset) | Excluded(offset) => self.end_epoch().unwrap() - *offset,
        };

        self.filter_by_epoch(start..end)
    }

    /// Returns a new tracking arc that only contains measurements from the desired tracker.
    pub fn filter_by_tracker(mut self, tracker: String) -> Self {
        self.measurements = self
            .measurements
            .iter()
            .filter_map(|(epoch, msr)| {
                if msr.tracker == tracker {
                    Some((*epoch, msr.clone()))
                } else {
                    None
                }
            })
            .collect::<BTreeMap<Epoch, Measurement>>();
        self
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
    pub fn downsample(self, target_step: Duration) -> Self {
        if self.is_empty() {
            return self;
        }
        let current_step = self.min_duration_sep().unwrap();

        if current_step >= target_step {
            warn!("cannot downsample tracking data from {current_step} to {target_step} (that would be upsampling)");
            return self;
        }

        let current_hz = 1.0 / current_step.to_seconds();
        let target_hz = 1.0 / target_step.to_seconds();

        // Simple moving average as low-pass filter
        let window_size = (current_hz / target_hz).round() as usize;

        info!("downsampling tracking data from {current_step} ({current_hz:.6} Hz) to {target_step} ({target_hz:.6} Hz) (N = {window_size})");

        let mut result = TrackingDataArc {
            source: self.source.clone(),
            ..Default::default()
        };

        let measurements: Vec<_> = self.measurements.iter().collect();

        for (i, (epoch, _)) in measurements.iter().enumerate().step_by(window_size) {
            let start = if i >= window_size / 2 {
                i - window_size / 2
            } else {
                0
            };
            let end = (i + window_size / 2 + 1).min(measurements.len());
            let window = &measurements[start..end];

            let mut filtered_measurement = Measurement {
                tracker: window[0].1.tracker.clone(),
                epoch: **epoch,
                data: IndexMap::new(),
            };

            // Apply moving average filter for each measurement type
            for mtype in self.unique_types() {
                let sum: f64 = window.iter().filter_map(|(_, m)| m.data.get(&mtype)).sum();
                let count = window
                    .iter()
                    .filter(|(_, m)| m.data.contains_key(&mtype))
                    .count();

                if count > 0 {
                    filtered_measurement.data.insert(mtype, sum / count as f64);
                }
            }

            result.measurements.insert(**epoch, filtered_measurement);
        }
        result
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
