use super::rand_pcg::Pcg64Mcg;
use crate::od::measurement::GroundStation;
use crate::time::{Duration, Epoch};
use std::default::Default;
// Should be a scheduler
/// A structure defining how a ground station is configured for a given tracking pass.
#[derive(Debug, Clone)]
pub struct GSSetup {
    /// The ground station itself
    pub gs: GroundStation,
    /// Availability of the ground station
    pub availability: Availability,
    /// Cadence of the station
    pub cadence: Cadence,
    /// Measurement integration method
    pub integration: Integration,
    /// The seed used to add noise to the measurements
    pub seed: u64,
}

/// Availability of the ground station
#[derive(Debug, Clone)]
pub enum Availability {
    /// Ground station is always available.
    Always,
    /// Ground station is never available (useful for checking impact of unavailability).
    Never,
    /// Ground station is only available between those epochs (included).
    Specific { from: Epoch, until: Epoch },
}

impl Default for Availability {
    fn default() -> Self {
        Self::Always
    }
}

/// Cadence of ground station.
#[derive(Debug, Clone)]
pub enum Cadence {
    /// When tracking and in visibility, the ground station will always be able to generate measurements.
    AlwaysOn,
    /// The station will be turned on (and generate measurements) for the first duration starting from the availability and then off.
    /// For example, (2*Unit::Hour, 10*Unit::Hour) here means that the station is will generate measurements for the first
    /// two hours after the start of its availability and then shut off for the following 10 hours.
    FromAvailability { on: Duration, off: Duration },
    /// Same as FromAvailability but the start/end times are computed from when the spacecraft is in visibility.
    FromTrackingStart { on: Duration, off: Duration },
}

impl Default for Cadence {
    fn default() -> Self {
        Self::AlwaysOn
    }
}

/// Integration specifies whether the measurements are integrated over time and if so at what frequency.
/// The measurements are then averaged (as is done with NASA DSN).
#[derive(Debug, Clone)]
pub enum Integration {
    /// There is no integration, the measurement is assumed instantaneous.
    Instantaneous,
    /// The measurements are integrated over the specified duration and at the specified sampling frequency.
    TimeAndSampling {
        duration: Duration,
        sampling_hz: f64,
        stamping: IntegrationTimestamp,
    },
}

impl Default for Integration {
    fn default() -> Self {
        Self::Instantaneous
    }
}

/// IntegrationTimestamp specifies where the measurement is time stamped if it's integrated over some time.
#[derive(Debug, Clone)]
pub enum IntegrationTimestamp {
    /// Measurement is time stamped at the start of the measurement integration.
    Start,
    /// Measurement is time stamped at the end of the measurement integration.
    End,
    /// Measurement is time stamped at the middle of the measurement integration.
    /// This _will_ lead to oddly specific measurement times if the integration time is not a multiple of two.
    Middle,
}

impl Default for IntegrationTimestamp {
    fn default() -> Self {
        Self::Middle
    }
}
