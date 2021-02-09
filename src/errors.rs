use std::error::Error;
use std::fmt;

#[derive(Clone, PartialEq, Debug)]
pub enum NyxError {
    /// STM is singular, check the automatic differentiation function and file a bug
    SingularStateTransitionMatrix,
    /// Fuel exhausted error, try running without fuel depletion, and then adding it.
    FuelExhausted,
    /// Propagation event not triggered withinin the max propagation time
    ConditionNeverTriggered,
    /// Propagation event not hit enough times (requested, found).
    UnsufficientTriggers(usize, usize),
    /// Maximum iterations reached, value corresponds to the number of iterations used
    MaxIterReached(usize),
    /// The STM was not updated prior to requesting a filter update
    StateTransitionMatrixNotUpdated,
    /// The operation was expecting the state to have an STM, but it isn't present.
    StateTransitionMatrixUnset,
    /// The sensitivity matrix was not updated prior to requesting a filter measurement update
    SensitivityNotUpdated,
    /// Kalman gain is singular, file a bug if this is encountered
    SingularKalmanGain,
    /// Covariance is singular
    SingularCovarianceMatrix,
    TargetsTooClose,
    LambertNotReasonablePhi,
    LambertMultiRevNotSupported,
    /// Returns this error if the partials for this model are not defined, thereby preventing the computation of the STM
    PartialsUndefined,
    LoadingError(String),
    ObjectNotFound(String),
    NoInterpolationData(String),
    InvalidInterpolationData(String),
    OutOfInterpolationWindow(String),
    NoStateData(String),
    DisjointFrameOrientations(String, String),
    /// When there is a controller but there isn't any thruster available
    CtrlExistsButNoThrusterAvail,
    /// The control vector returned by a controller must be a unit vector. Use the throttle() function to specify the amount.
    CtrlNotAUnitVector(f64),
    /// The control throttle range must be between 0.0 and 1.0 (both included) as it represents a percentage.
    CtrlThrottleRangeErr(f64),
    /// An objective based analysis or control was attempted, but no objective was defined.
    NoObjectiveDefined,
    /// Some custom error for new dynamics
    CustomError(String),
}

impl fmt::Display for NyxError {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            Self::StateTransitionMatrixNotUpdated => {
                write!(f, "STM was not updated prior to time or measurement update")
            }
            Self::SensitivityNotUpdated => write!(
                f,
                "The measurement matrix H_tilde was not updated prior to measurement update"
            ),
            Self::SingularKalmanGain => write!(
                f,
                "Gain could not be computed because H*P_bar*H + R is singular"
            ),
            Self::SingularStateTransitionMatrix => write!(
                f,
                "STM is singular, propagation or smoothing cannot proceed"
            ),
            Self::SingularCovarianceMatrix => {
                write!(f, "Covariance is singular, smoothing cannot proceed")
            }
            Self::FuelExhausted => write!(
                f,
                "Spacecraft fuel exhausted, disable fuel depletion and place maneuvers"
            ),
            Self::ConditionNeverTriggered => write!(
                f,
                "Try increasing the search space, i.e. increase the maximum propagation time"
            ),
            Self::TargetsTooClose => write!(f, "Lambert too close: Δν ~=0 and A ~=0"),
            Self::LambertMultiRevNotSupported => {
                write!(f, "Use the Izzo algorithm for multi-rev transfers")
            }
            Self::LambertNotReasonablePhi => {
                write!(f, "No reasonable phi found to connect both radii")
            }
            _ => write!(f, "{:?}", self),
        }
    }
}

impl Error for NyxError {}
