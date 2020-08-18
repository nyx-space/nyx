use crate::celestia::Frame;
use crate::dimensions::allocator::Allocator;
use crate::dimensions::{DefaultAllocator, DimName, MatrixMN, VectorN};
use crate::time::Epoch;

use std::fmt;

#[derive(Clone)]
pub struct SNC<A: DimName>
where
    DefaultAllocator: Allocator<f64, A> + Allocator<f64, A, A>,
{
    /// Time at which this SNC starts to become applicable
    pub start_time: Option<Epoch>,
    /// Specify the frame of this SNC -- CURRENTLY UNIMPLEMENTED
    pub frame: Option<Frame>,
    /// Enables state noise compensation (process noise) only be applied if the time between measurements is less than the disable_time_s amount in seconds
    pub disable_time_s: f64,
    // Stores the initial epoch when the SNC is requested, needed for decay. Kalman filter will edit this automatically.
    pub init_epoch: Option<Epoch>,
    diag: VectorN<f64, A>,
    decay_diag: Option<Vec<f64>>,
    // Stores the previous epoch of the SNC request, needed for disable time
    pub prev_epoch: Option<Epoch>,
}

impl<A> fmt::Debug for SNC<A>
where
    A: DimName,
    DefaultAllocator: Allocator<f64, A> + Allocator<f64, A, A>,
{
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if let Some(decay) = &self.decay_diag {
            let mut fmt_cov = Vec::with_capacity(A::dim());
            for i in 0..A::dim() {
                fmt_cov.push(format!("{:.3e} × exp(- {:.3e} t)", self.diag[i], decay[i]));
            }
            write!(
                f,
                "SNC: diag({}) {}",
                fmt_cov.join(","),
                if let Some(start) = self.start_time {
                    format!("starting at {}", start.as_gregorian_utc_str())
                } else {
                    "".to_string()
                }
            )
        } else {
            let mut fmt_cov = Vec::with_capacity(A::dim());
            for i in 0..A::dim() {
                fmt_cov.push(format!("{:.3e}", self.diag[i]));
            }
            write!(
                f,
                "SNC: diag({}) {}",
                fmt_cov.join(","),
                if let Some(start) = self.start_time {
                    format!("starting at {}", start.as_gregorian_utc_str())
                } else {
                    "".to_string()
                }
            )
        }
    }
}

impl<A: DimName> SNC<A>
where
    DefaultAllocator: Allocator<f64, A> + Allocator<f64, A, A>,
{
    /// Initialize a state noise compensation structure from the diagonal values
    pub fn from_diagonal(disable_time_s: f64, values: &[f64]) -> Self {
        assert_eq!(
            values.len(),
            A::dim(),
            "Not enough values for the size of the SNC matrix"
        );

        let mut diag = VectorN::zeros();
        for (i, v) in values.iter().enumerate() {
            diag[i] = *v;
        }

        Self {
            diag,
            disable_time_s,
            start_time: None,
            frame: None,
            decay_diag: None,
            init_epoch: None,
            prev_epoch: None,
        }
    }

    /// Initialize an SNC with a time at which it should start
    pub fn with_start_time(disable_time_s: f64, values: &[f64], start_time: Epoch) -> Self {
        let mut me = Self::from_diagonal(disable_time_s, values);
        me.start_time = Some(start_time);
        me
    }

    /// Initialize an exponentially decaying SNC with initial SNC and decay constants.
    /// Decay constants in seconds since start of the tracking pass.
    pub fn with_decay(disable_time_s: f64, initial_snc: &[f64], decay_constants_s: &[f64]) -> Self {
        assert_eq!(
            decay_constants_s.len(),
            A::dim(),
            "Not enough decay constants for the size of the SNC matrix"
        );

        let mut me = Self::from_diagonal(disable_time_s, initial_snc);
        me.decay_diag = Some(decay_constants_s.to_vec());
        me
    }

    /// Returns the SNC matrix (_not_ incl. Gamma matrix approximation) at the provided Epoch.
    /// May be None if:
    ///  1. Start time of this matrix is _after_ epoch
    ///  2. Time between epoch and previous epoch (set in the Kalman filter!) is longer than disabling time
    pub fn to_matrix(&self, epoch: Epoch) -> Option<MatrixMN<f64, A, A>> {
        if let Some(start_time) = self.start_time {
            if start_time > epoch {
                // This SNC applies only later
                return None;
            }
        }

        // Check the disable time, and return no SNC if the previous SNC was computed too long ago
        if let Some(prev_epoch) = self.prev_epoch {
            if epoch - prev_epoch > self.disable_time_s {
                return None;
            }
        }
        // Build a static matrix
        let mut snc = MatrixMN::<f64, A, A>::zeros();
        for i in 0..self.diag.nrows() {
            snc[(i, i)] = self.diag[i];
        }

        if let Some(decay) = &self.decay_diag {
            // Let's apply the decay to the diagonals
            let total_delta_t = epoch - self.init_epoch.unwrap();
            for i in 0..self.diag.nrows() {
                snc[(i, i)] *= (-decay[i] * total_delta_t).exp();
            }
        }

        Some(snc)
    }
}
