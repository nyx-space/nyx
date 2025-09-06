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
use std::collections::btree_map::Entry;
use std::collections::BTreeMap;
use std::ops::Bound::{self, Excluded, Included, Unbounded};
use std::ops::{Add, AddAssign, RangeBounds};

mod io_ccsds_tdm;
mod io_parquet;

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
pub struct TrackingDataArc {
    /// All measurements in this data arc
    pub measurements: BTreeMap<Epoch, Measurement>, // BUG: Consider a map of tracking to epoch!
    /// Source file if loaded from a file or saved to a file.
    pub source: Option<String>,
    /// Optionally provide a map of modulos (e.g. the RANGE_MODULO of CCSDS TDM).
    pub moduli: Option<IndexMap<MeasurementType, f64>>,
    /// Reject all of the measurements, useful for debugging passes.
    pub force_reject: bool,
}

impl TrackingDataArc {
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
            for msr in self.measurements.values_mut() {
                for (msr_type, modulus) in moduli {
                    if let Some(msr_value) = msr.data.get_mut(msr_type) {
                        *msr_value %= *modulus;
                    }
                }
            }
        }
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

    /// Returns a new tracking arc that only contains measurements of the provided type.
    pub fn filter_by_measurement_type(mut self, included_type: MeasurementType) -> Self {
        self.measurements.retain(|_epoch, msr| {
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
            .filter_map(|(epoch, msr)| {
                if msr.tracker != excluded_tracker {
                    Some((*epoch, msr.clone()))
                } else {
                    None
                }
            })
            .collect::<BTreeMap<Epoch, Measurement>>();
        self
    }

    /// Returns a new tracking arc that excludes measurements within the given epoch range.
    pub fn exclude_by_epoch<R: RangeBounds<Epoch>>(mut self, bound: R) -> Self {
        info!(
            "Excluding measurements from {:?} to {:?}",
            bound.start_bound(),
            bound.end_bound()
        );

        let mut new_measurements = BTreeMap::new();

        // Include entries before the excluded range
        new_measurements.extend(
            self.measurements
                .range((Bound::Unbounded, bound.start_bound()))
                .map(|(epoch, msr)| (*epoch, msr.clone())),
        );

        // Include entries after the excluded range
        new_measurements.extend(
            self.measurements
                .range((bound.end_bound(), Bound::Unbounded))
                .map(|(epoch, msr)| (*epoch, msr.clone())),
        );

        self.measurements = new_measurements;
        self
    }

    /// Returns a new tracking arc that contains measurements from all trackers except the one provided
    pub fn exclude_measurement_type(mut self, excluded_type: MeasurementType) -> Self {
        self.measurements = self
            .measurements
            .iter_mut()
            .map(|(epoch, msr)| {
                msr.data.retain(|msr_type, _| *msr_type != excluded_type);

                (*epoch, msr.clone())
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
            let start = i.saturating_sub(window_size / 2);
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
        for (epoch, msr) in rhs.measurements {
            if let Entry::Vacant(e) = self.measurements.entry(epoch) {
                e.insert(msr);
            } else {
                error!("merging tracking data with overlapping epoch is not supported");
            }
        }

        self
    }
}

impl AddAssign for TrackingDataArc {
    fn add_assign(&mut self, rhs: Self) {
        *self = self.clone() + rhs;
    }
}
