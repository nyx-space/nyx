/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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

extern crate meval;

use self::meval::{Context, Expr};
use crate::log::error;
use crate::na::Matrix3;
use crate::time::Epoch;
use crate::utils::{r1, r2, r3};
use std::cmp::PartialEq;
use std::collections::HashMap;
use std::fmt;
use std::io::{Error as IoError, ErrorKind as IoErrorKind};
use std::str::FromStr;

pub trait ParentRotation: Send + Sync + fmt::Debug {
    fn dcm_to_parent(&self, datetime: Epoch) -> Option<Matrix3<f64>>;
}

#[derive(Debug)]
pub struct NoRotation;
impl ParentRotation for NoRotation {
    fn dcm_to_parent(&self, _: Epoch) -> Option<Matrix3<f64>> {
        Some(Matrix3::identity())
    }
}

/// Defines an Euler rotation, angle must be in radians
#[derive(Clone, Copy, Debug)]
pub enum EulerRotation {
    R1(f64),
    R2(f64),
    R3(f64),
}

impl EulerRotation {
    pub fn r1_from_degrees(angle_deg: f64) -> Self {
        Self::R1(angle_deg.to_radians())
    }
    pub fn r2_from_degrees(angle_deg: f64) -> Self {
        Self::R2(angle_deg.to_radians())
    }
    pub fn r3_from_degrees(angle_deg: f64) -> Self {
        Self::R3(angle_deg.to_radians())
    }
    /// Get the DCM from this Euler rotation
    pub fn dcm(&self) -> Matrix3<f64> {
        match *self {
            Self::R1(angle) => r1(angle),
            Self::R2(angle) => r2(angle),
            Self::R3(angle) => r3(angle),
        }
    }
}

impl ParentRotation for EulerRotation {
    fn dcm_to_parent(&self, _: Epoch) -> Option<Matrix3<f64>> {
        Some(self.dcm())
    }
}

/// A fixed three-axis Euler rotation
#[derive(Debug)]
pub struct Euler3Axis {
    /// The first rotation (e.g. R3)
    pub first: EulerRotation,
    /// The second rotation (e.g. R1)
    pub second: EulerRotation,
    /// The third and final rotation (e.g. R3, to complete a 3-1-1 rotation)
    pub third: EulerRotation,
}

impl ParentRotation for Euler3Axis {
    fn dcm_to_parent(&self, _: Epoch) -> Option<Matrix3<f64>> {
        Some(self.third.dcm() * self.second.dcm() * self.first.dcm())
    }
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum AngleUnit {
    Degrees,
    Radians,
}

impl FromStr for AngleUnit {
    type Err = IoError;
    fn from_str(s: &str) -> Result<Self, Self::Err> {
        if s.to_lowercase().starts_with("deg") {
            Ok(AngleUnit::Degrees)
        } else if s.to_lowercase().starts_with("rad") {
            Ok(AngleUnit::Radians)
        } else {
            Err(IoError::new(
                IoErrorKind::InvalidData,
                format!("unknown angles `{}`", s),
            ))
        }
    }
}

/// A time varying three-axis Euler rotation
#[derive(Clone)]
pub struct Euler3AxisDt
where
    Self: Send + Sync,
{
    pub base_context: HashMap<String, String>,
    pub rot_order: [(EulerRotation, Expr); 3],
    pub unit: AngleUnit,
    pub is_ra_dec_w: bool,
}

impl Euler3AxisDt {
    /// Initialize a new time varying transformation.
    /// Reserved keywords in the context are "T" for centuries past 2000 Jan 1 12h TBD
    /// epoch (JDE 2451545.0), and "d" for days since that epoch.
    fn new(
        first_rot: (EulerRotation, Expr),
        second_rot: (EulerRotation, Expr),
        third_rot: (EulerRotation, Expr),
        context: HashMap<String, String>,
        unit: AngleUnit,
        is_ra_dec_w: bool,
    ) -> Self {
        let rot_order = [first_rot, second_rot, third_rot];
        Self {
            base_context: context,
            rot_order,
            unit,
            is_ra_dec_w,
        }
    }

    /// Specify how to compute this frame from the provided Euler angles and their time varying expressions.
    /// Note that these angles define how to go from THIS frame TO the PARENT frame (e.g. Sun fixed to ICRF).
    pub fn from_euler_angles(
        first_rot: (EulerRotation, Expr),
        second_rot: (EulerRotation, Expr),
        third_rot: (EulerRotation, Expr),
        context: HashMap<String, String>,
        unit: AngleUnit,
    ) -> Self {
        Self::new(first_rot, second_rot, third_rot, context, unit, false)
    }

    /// A time varying Right ascension, Declination, and W frame
    /// Conversion TO parent frame (e.g. Sun body to ICRF) defined as:
    /// R3(-(alpha-90 deg)) * R1(delta - 90 deg) * R3(-W)
    /// Where alpha is the right ascension and delta the declination
    pub fn from_ra_dec_w(
        alpha_right_asc: Expr,
        delta_declin: Expr,
        w_twist: Expr,
        context: HashMap<String, String>,
        unit: AngleUnit,
    ) -> Self {
        // We initialize this backward on purpose, cf. https://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/rotation.html#Working%20with%20RA,%20Dec%20and%20Twist
        Self::new(
            (EulerRotation::R3(0.0), alpha_right_asc),
            (EulerRotation::R1(0.0), delta_declin),
            (EulerRotation::R3(0.0), w_twist),
            context,
            unit,
            true,
        )
    }
}

impl fmt::Debug for Euler3AxisDt {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{:?}-{:?}-{:?}",
            self.rot_order[0].0, self.rot_order[1].0, self.rot_order[2].0
        )
    }
}

impl ParentRotation for Euler3AxisDt {
    fn dcm_to_parent(&self, datetime: Epoch) -> Option<Matrix3<f64>> {
        let days_d = datetime.as_tdb_days_since_j2000();
        let centuries_t = datetime.as_tdb_centuries_since_j2000();
        // Let's create a new context, add the time variables, and compute the rotation's context
        // Note: this will be resolved with Chi/racr files.
        let mut ctx = Context::default();
        ctx.var("d", days_d);
        ctx.var("T", centuries_t);
        for (var, expr_str) in &self.base_context {
            let as_expr: Expr = expr_str.parse().unwrap();
            ctx.var(
                var.to_owned(),
                as_expr.eval_with_context(&ctx).unwrap_or_else(|_| {
                    panic!("Could not evaluate variable `{}` as `{}`", var, expr_str)
                }),
            );
        }
        let mut dcm = Matrix3::identity();
        for (angle_no, (rot, expr)) in self.rot_order.iter().rev().enumerate() {
            // Compute the correct angle
            match expr.eval_with_context(&ctx) {
                Ok(eval_angle) => {
                    // Convert the angle to radians if needed
                    let mut angle = if self.unit == AngleUnit::Degrees {
                        eval_angle.to_radians()
                    } else {
                        eval_angle
                    };
                    // Apply the angle transformation if needed
                    if self.is_ra_dec_w {
                        if angle_no == 1 {
                            angle = std::f64::consts::FRAC_PI_2 - angle;
                        } else if angle_no == 2 {
                            angle += std::f64::consts::FRAC_PI_2;
                        }
                    }
                    let rot_with_angl = match rot {
                        EulerRotation::R1(_) => EulerRotation::R1(angle),
                        EulerRotation::R2(_) => EulerRotation::R2(angle),
                        EulerRotation::R3(_) => EulerRotation::R3(angle),
                    };
                    dcm *= rot_with_angl.dcm();
                }
                Err(e) => {
                    error!("{}", e);
                    // Stop here if something when wrong
                    return None;
                }
            }
        }
        Some(dcm)
    }
}

#[test]
fn test_angle_unit_deser() {
    use std::str::FromStr;
    assert_eq!(AngleUnit::from_str("DeGrEes").unwrap(), AngleUnit::Degrees);
    assert_eq!(AngleUnit::from_str("RaDiaNs").unwrap(), AngleUnit::Radians);
    assert!(AngleUnit::from_str("Gradian").is_err());
}
