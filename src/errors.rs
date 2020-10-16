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
    /// The sensitivity matrix was not updated prior to requesting a filter measurement update
    SensitivityNotUpdated,
    /// Kalman gain is singular, file a bug if this is encountered
    SingularKalmanGain,
    /// Covariance is singular
    SingularCovarianceMatrix,
    /// Some custom error for new dynamics
    CustomError(String),
}

impl fmt::Display for NyxError {
    // Prints the Keplerian orbital elements with units
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match *self {
            NyxError::StateTransitionMatrixNotUpdated => {
                write!(f, "STM was not updated prior to time or measurement update")
            }
            NyxError::SensitivityNotUpdated => write!(
                f,
                "The measurement matrix H_tilde was not updated prior to measurement update"
            ),
            NyxError::SingularKalmanGain => write!(
                f,
                "Gain could not be computed because H*P_bar*H + R is singular"
            ),
            NyxError::SingularStateTransitionMatrix => write!(
                f,
                "STM is singular, propagation or smoothing cannot proceed"
            ),
            NyxError::SingularCovarianceMatrix => {
                write!(f, "Covariance is singular, smoothing cannot proceed")
            }
            NyxError::FuelExhausted => write!(
                f,
                "Spacecraft fuel exhausted, disable fuel depletion and place maneuvers"
            ),
            NyxError::ConditionNeverTriggered => write!(
                f,
                "Try increasing the search space, i.e. increase the maximum propagation time"
            ),
            _ => write!(f, "{:?}", self),
        }
    }
}

impl Error for NyxError {}
