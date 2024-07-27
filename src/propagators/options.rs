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

use std::fmt;

use crate::time::{Duration, Unit};

use super::{ErrorCtrl, RSSCartesianStep};
use typed_builder::TypedBuilder;

/// PropOpts stores the integrator options, including the minimum and maximum step sizes, and the
/// max error size.
///
/// Note that different step sizes and max errors are only used for adaptive
/// methods. To use a fixed step integrator, initialize the options using `with_fixed_step`, and
/// use whichever adaptive step integrator is desired.  For example, initializing an RK45 with
/// fixed step options will lead to an RK4 being used instead of an RK45.
#[derive(Clone, Copy, Debug, TypedBuilder)]
#[builder(doc)]
pub struct PropOpts<E: ErrorCtrl> {
    #[builder(default_code = "60.0 * Unit::Second")]
    pub init_step: Duration,
    #[builder(default_code = "0.001 * Unit::Second")]
    pub min_step: Duration,
    #[builder(default_code = "2700.0 * Unit::Second")]
    pub max_step: Duration,
    #[builder(default = 1e-12)]
    pub tolerance: f64,
    #[builder(default = 50)]
    pub attempts: u8,
    #[builder(default = false)]
    pub fixed_step: bool,
    pub error_ctrl: E,
}

impl<E: ErrorCtrl> PropOpts<E> {
    /// `with_adaptive_step` initializes an `PropOpts` such that the integrator is used with an
    ///  adaptive step size. The number of attempts is currently fixed to 50 (as in GMAT).
    pub fn with_adaptive_step(
        min_step: Duration,
        max_step: Duration,
        tolerance: f64,
        error_ctrl: E,
    ) -> Self {
        PropOpts {
            init_step: max_step,
            min_step,
            max_step,
            tolerance,
            attempts: 50,
            fixed_step: false,
            error_ctrl,
        }
    }

    pub fn with_adaptive_step_s(
        min_step: f64,
        max_step: f64,
        tolerance: f64,
        error_ctrl: E,
    ) -> Self {
        Self::with_adaptive_step(
            min_step * Unit::Second,
            max_step * Unit::Second,
            tolerance,
            error_ctrl,
        )
    }

    /// Returns a string with the information about these options
    pub fn info(&self) -> String {
        format!("{self}")
    }

    /// Set the maximum step size and sets the initial step to that value if currently greater
    pub fn set_max_step(&mut self, max_step: Duration) {
        if self.init_step > max_step {
            self.init_step = max_step;
        }
        self.max_step = max_step;
    }

    /// Set the minimum step size and sets the initial step to that value if currently smaller
    pub fn set_min_step(&mut self, min_step: Duration) {
        if self.init_step < min_step {
            self.init_step = min_step;
        }
        self.min_step = min_step;
    }
}

impl<E: ErrorCtrl> fmt::Display for PropOpts<E> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        if self.fixed_step {
            write!(f, "fixed step: {:e}", self.min_step,)
        } else {
            write!(
                f,
                "min_step: {:e}, max_step: {:e}, tol: {:e}, attempts: {}",
                self.min_step, self.max_step, self.tolerance, self.attempts,
            )
        }
    }
}

impl PropOpts<RSSCartesianStep> {
    /// `with_fixed_step` initializes an `PropOpts` such that the integrator is used with a fixed
    ///  step size.
    pub fn with_fixed_step(step: Duration) -> Self {
        PropOpts {
            init_step: step,
            min_step: step,
            max_step: step,
            tolerance: 0.0,
            fixed_step: true,
            attempts: 0,
            error_ctrl: RSSCartesianStep {},
        }
    }

    pub fn with_fixed_step_s(step: f64) -> Self {
        Self::with_fixed_step(step * Unit::Second)
    }

    /// Returns the default options with a specific tolerance.
    #[allow(clippy::field_reassign_with_default)]
    pub fn with_tolerance(tolerance: f64) -> Self {
        let mut opts = Self::default();
        opts.tolerance = tolerance;
        opts
    }

    /// Creates a propagator with the provided max step, and sets the initial step to that value as well.
    #[allow(clippy::field_reassign_with_default)]
    pub fn with_max_step(max_step: Duration) -> Self {
        let mut opts = Self::default();
        opts.set_max_step(max_step);
        opts
    }
}

impl Default for PropOpts<RSSCartesianStep> {
    /// `default` returns the same default options as GMAT.
    fn default() -> PropOpts<RSSCartesianStep> {
        PropOpts {
            init_step: 60.0 * Unit::Second,
            min_step: 0.001 * Unit::Second,
            max_step: 2700.0 * Unit::Second,
            tolerance: 1e-12,
            attempts: 50,
            fixed_step: false,
            error_ctrl: RSSCartesianStep {},
        }
    }
}

#[test]
fn test_options() {
    use super::error_ctrl::RSSStep;

    let opts = PropOpts::with_fixed_step_s(1e-1);
    assert_eq!(opts.min_step, 1e-1 * Unit::Second);
    assert_eq!(opts.max_step, 1e-1 * Unit::Second);
    assert!(opts.tolerance.abs() < f64::EPSILON);
    assert!(opts.fixed_step);

    let opts = PropOpts::with_adaptive_step_s(1e-2, 10.0, 1e-12, RSSStep {});
    assert_eq!(opts.min_step, 1e-2 * Unit::Second);
    assert_eq!(opts.max_step, 10.0 * Unit::Second);
    assert!((opts.tolerance - 1e-12).abs() < f64::EPSILON);
    assert!(!opts.fixed_step);

    let opts: PropOpts<RSSCartesianStep> = Default::default();
    assert_eq!(opts.init_step, 60.0 * Unit::Second);
    assert_eq!(opts.min_step, 0.001 * Unit::Second);
    assert_eq!(opts.max_step, 2700.0 * Unit::Second);
    assert!((opts.tolerance - 1e-12).abs() < f64::EPSILON);
    assert_eq!(opts.attempts, 50);
    assert!(!opts.fixed_step);

    let opts = PropOpts::with_max_step(1.0 * Unit::Second);
    assert_eq!(opts.init_step, 1.0 * Unit::Second);
    assert_eq!(opts.min_step, 0.001 * Unit::Second);
    assert_eq!(opts.max_step, 1.0 * Unit::Second);
    assert!((opts.tolerance - 1e-12).abs() < f64::EPSILON);
    assert_eq!(opts.attempts, 50);
    assert!(!opts.fixed_step);
}
