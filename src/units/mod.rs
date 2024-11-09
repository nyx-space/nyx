pub mod mass;

type Base = f64;

macro_rules! unit {
    ($me:ident, $abbrev:literal) => {
        impl $me {
            pub fn new(base: $crate::units::Base) -> Self {
                $me(base)
            }

            fn inner(&self) -> Base {
                self.0
            }
        }

        impl ::std::fmt::Debug for $me {
            fn fmt(
                &self,
                f: &mut ::std::fmt::Formatter<'_>,
            ) -> ::std::result::Result<(), ::std::fmt::Error> {
                ::std::write!(f, "{}({})", ::std::stringify!($me), self.0)
            }
        }

        impl ::std::fmt::Display for $me {
            fn fmt(
                &self,
                f: &mut ::std::fmt::Formatter<'_>,
            ) -> ::std::result::Result<(), ::std::fmt::Error> {
                ::std::write!(f, "{} {}", self.0, $abbrev)
            }
        }

        impl ::std::ops::Add for $me {
            type Output = Self;

            fn add(self, rhs: Self) -> Self::Output {
                $me(self.0.add(rhs.0))
            }
        }

        impl ::std::ops::AddAssign for $me {
            fn add_assign(&mut self, rhs: Self) {
                (&mut self.0).add_assign(rhs.0);
            }
        }

        impl ::std::ops::Sub for $me {
            type Output = Self;

            fn sub(self, rhs: Self) -> Self::Output {
                $me(self.0.sub(rhs.0))
            }
        }

        impl ::std::ops::SubAssign for $me {
            fn sub_assign(&mut self, rhs: Self) {
                (&mut self.0).sub_assign(rhs.0);
            }
        }

        impl ::std::ops::Mul for $me {
            type Output = Self;

            fn mul(self, rhs: Self) -> Self::Output {
                $me(self.0.mul(rhs.0))
            }
        }

        impl ::std::ops::MulAssign for $me {
            fn mul_assign(&mut self, rhs: Self) {
                (&mut self.0).mul_assign(rhs.0);
            }
        }

        impl ::std::ops::Div for $me {
            type Output = Self;

            fn div(self, rhs: Self) -> Self::Output {
                $me(self.0.div(rhs.0))
            }
        }

        impl ::std::ops::DivAssign for $me {
            fn div_assign(&mut self, rhs: Self) {
                (&mut self.0).div_assign(rhs.0);
            }
        }
    };
}

macro_rules! linear_conversion {
    ($factor:literal, $from:ident, $into:ident) => {
        impl std::convert::From<$from> for $into {
            fn from(value: $from) -> Self {
                let converted = $factor * value.inner();
                $into::new(converted)
            }
        }

        impl ::std::ops::Add<$from> for $into {
            type Output = Self;

            fn add(self, rhs: $from) -> Self::Output {
                let rhs: $into = rhs.into();
                self.add(rhs)
            }
        }

        impl ::std::ops::AddAssign<$from> for $into {
            fn add_assign(&mut self, rhs: $from) {
                let rhs: $into = rhs.into();
                self.add_assign(rhs);
            }
        }

        impl ::std::ops::Sub<$from> for $into {
            type Output = Self;

            fn sub(self, rhs: $from) -> Self::Output {
                let rhs: $into = rhs.into();
                self.sub(rhs)
            }
        }

        impl ::std::ops::SubAssign<$from> for $into {
            fn sub_assign(&mut self, rhs: $from) {
                let rhs: $into = rhs.into();
                self.sub_assign(rhs);
            }
        }

        impl ::std::ops::Mul<$from> for $into {
            type Output = Self;

            fn mul(self, rhs: $from) -> Self::Output {
                let rhs: $into = rhs.into();
                self.mul(rhs)
            }
        }

        impl ::std::ops::MulAssign<$from> for $into {
            fn mul_assign(&mut self, rhs: $from) {
                let rhs: $into = rhs.into();
                self.mul_assign(rhs)
            }
        }

        impl ::std::ops::Div<$from> for $into {
            type Output = Self;

            fn div(self, rhs: $from) -> Self::Output {
                let rhs: $into = rhs.into();
                self.div(rhs)
            }
        }

        impl ::std::ops::DivAssign<$from> for $into {
            fn div_assign(&mut self, rhs: $from) {
                let rhs: $into = rhs.into();
                self.div_assign(rhs)
            }
        }
    };
}

pub(crate) use linear_conversion;
pub(crate) use unit;
