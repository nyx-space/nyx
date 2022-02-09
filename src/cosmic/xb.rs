/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

#[derive(Clone, PartialEq, prost::Message)]
pub struct Metadata {
    #[prost(message, optional, tag = "1")]
    pub version: ::core::option::Option<metadata::Version>,
    #[prost(string, tag = "2")]
    pub publisher: prost::alloc::string::String,
    #[prost(message, optional, tag = "3")]
    pub date: ::core::option::Option<metadata::CeDate>,
    #[prost(string, tag = "4")]
    pub file_version: prost::alloc::string::String,
    #[prost(bool, tag = "5")]
    pub proprietary: bool,
    #[prost(string, tag = "6")]
    pub comments: prost::alloc::string::String,
}
/// Nested message and enum types in `Metadata`.
pub mod metadata {
    #[derive(Clone, PartialEq, prost::Message)]
    pub struct CeDate {
        #[prost(uint32, tag = "1")]
        pub year: u32,
        #[prost(uint32, tag = "2")]
        pub month: u32,
        #[prost(uint32, tag = "3")]
        pub day: u32,
    }
    #[derive(Clone, PartialEq, prost::Message)]
    pub struct Version {
        #[prost(uint32, tag = "1")]
        pub major: u32,
        #[prost(uint32, tag = "2")]
        pub minor: u32,
        #[prost(uint32, tag = "3")]
        pub patch: u32,
    }
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Frame {
    /// Center object of this frame, may be a celestial body, a
    /// spacecraft, a ground station, an instrument, etc.
    /// Note that this corresponds to the full path of the object, e.g.
    /// "/EarthMoonBarycenter/Earth/CAPSTONE/Camera1".
    #[prost(string, tag = "1")]
    pub center: prost::alloc::string::String,
    /// Orientation used by this frame (e.g. "J2000", "IAU Earth", etc.)
    #[prost(string, tag = "2")]
    pub orientation: prost::alloc::string::String,
}
/// A constant value, may be used to specify celestial object properties (e.g.
/// GM), or spacecraft properties (e.g. mass).
#[derive(Clone, PartialEq, prost::Message)]
pub struct Constant {
    #[prost(double, tag = "1")]
    pub value: f64,
    #[prost(enumeration = "Unit", tag = "2")]
    pub unit: i32,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Vector {
    #[prost(double, tag = "1")]
    pub x: f64,
    #[prost(double, tag = "2")]
    pub y: f64,
    #[prost(double, tag = "3")]
    pub z: f64,
    #[prost(enumeration = "Unit", tag = "4")]
    pub unit: i32,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Vector4 {
    #[prost(double, tag = "1")]
    pub x: f64,
    #[prost(double, tag = "2")]
    pub y: f64,
    #[prost(double, tag = "3")]
    pub z: f64,
    #[prost(double, tag = "4")]
    pub w: f64,
}
/// Used to store coefficients to interpolate a three dimensional vector.
#[derive(Clone, PartialEq, prost::Message)]
pub struct VectorCoefficients {
    #[prost(double, repeated, tag = "1")]
    pub x: prost::alloc::vec::Vec<f64>,
    #[prost(double, repeated, tag = "2")]
    pub y: prost::alloc::vec::Vec<f64>,
    #[prost(double, repeated, tag = "3")]
    pub z: prost::alloc::vec::Vec<f64>,
}
/// Defines an Epoch. For example, DE438 uses the ET system with an DaysJDE
/// representation.
#[derive(Clone, PartialEq, prost::Message)]
pub struct Epoch {
    #[prost(enumeration = "TimeSystem", tag = "1")]
    pub ts: i32,
    #[prost(enumeration = "TimeRepr", tag = "2")]
    pub repr: i32,
    /// Days since the epoch of the repr used. Recommendation: use
    #[prost(int32, tag = "3")]
    pub days: i32,
    /// days for any `seconds` value greater than 86400 seconds.
    ///
    /// Seconds are stored in a double to trivial support
    #[prost(double, tag = "4")]
    pub seconds: f64,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Equation {
    #[prost(string, tag = "1")]
    pub expression: prost::alloc::string::String,
    #[prost(enumeration = "Unit", tag = "2")]
    pub unit: i32,
    #[prost(map = "string, message", tag = "3")]
    pub constants: ::std::collections::HashMap<prost::alloc::string::String, Constant>,
    #[prost(map = "string, message", tag = "4")]
    pub context: ::std::collections::HashMap<prost::alloc::string::String, Equation>,
}
/// Grouping all units together allow setting a generic unit to a Vector.
/// Note that we try to group the enum numbering in sections. This provides some
/// lee-way for extending the specs to other units, while maintaining some
/// numbering scheme.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, prost::Enumeration)]
#[repr(i32)]
#[allow(clippy::enum_variant_names)]
pub enum Unit {
    Dimensionless = 0,
    /// Distance units
    Au = 1,
    Km = 2,
    M = 3,
    Cm = 4,
    Mm = 5,
    /// Velocity units
    KmS = 10,
    MS = 11,
    CmS = 12,
    MmS = 13,
    /// Acceleration units
    KmS2 = 20,
    MS2 = 21,
    CmS2 = 22,
    MmS2 = 23,
    /// Angular units
    Deg = 30,
    Rad = 31,
    /// Angular velocity units
    DegS = 32,
    RadS = 33,
    /// Angular acceleration units
    DegS2 = 34,
    RadS2 = 35,
    /// Frequency units
    Hz = 40,
    KHz = 41,
    MHz = 42,
    Ghz = 43,
    /// Duration units -- anything above days is poorly defined, so we ignore that
    Ns = 50,
    Us = 51,
    Ms = 52,
    S = 53,
    Min = 54,
    H = 55,
    Days = 56,
    /// Other units
    CustomUnit = 60,
    /// Used for gravitational parameters
    Km3S2 = 61,
    N = 62,
}
/// Allows specifying the time system used
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, prost::Enumeration)]
#[repr(i32)]
pub enum TimeSystem {
    /// Temps Atomique International
    Tai = 0,
    /// Ephemeris Time (very slightly different from IERS TDB)
    Et = 1,
    /// Terrestrial Time
    Tt = 2,
    /// Universal Coordinated Time (assumes ALL leap seconds)
    Utc = 3,
    /// Barycentric Dynamical Time (IERS TDB)
    Tdb = 4,
    CustomTs = 5,
}
/// Allow specifying the time representation used
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, prost::Enumeration)]
#[repr(i32)]
#[allow(clippy::enum_variant_names)]
pub enum TimeRepr {
    /// Seconds past 01 January 1900 at MIDNIGHT (or TAI Epoch)
    SecondsJ1900 = 0,
    /// Days of 86400.0 seconds past of the Epoch of the time system
    DaysJ1900 = 1,
    /// Seconds past 01 January 2000 at NOON (NOT midnight, hence
    SecondsJ2k = 2,
    /// different notation)
    ///
    /// Days of 86400.0 seconds past 01 January 2000 at NOON (NOT
    DaysJ2k = 3,
    /// midnight, hence different notation)
    ///
    /// Days past the Modified Julian Date Epoch defined as 17
    DaysMjd = 4,
    /// November 1858 at zero hours.
    ///
    /// Days past the Julian Date Epoch.
    DaysJde = 5,
    CustomTimeRepr = 6,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, prost::Enumeration)]
#[repr(i32)]
pub enum InterpType {
    /// Supported interpolations.
    Chebyshev = 0,
    Hermite = 1,
    Polynomial = 2,
    Lagrange = 3,
    Fourier = 4,
    /// NOTE: Requires additional communication between provider and implementer.
    CustomInterpolationType = 5,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Orientation {
    /// Frame of this orientation, with center and orientation
    #[prost(message, optional, tag = "1")]
    pub frame: ::core::option::Option<Frame>,
    /// A start time of all of the states, discrete or interpolated.
    #[prost(message, optional, tag = "2")]
    pub start_epoch: ::core::option::Option<Epoch>,
    #[prost(message, optional, tag = "3")]
    pub records: ::core::option::Option<AttitudeRegistry>,
    #[prost(message, optional, tag = "4")]
    pub interpolator: ::core::option::Option<AttitudeInterp>,
    /// A map of constant name to constant value and unit.
    #[prost(map = "string, message", tag = "5")]
    pub constants: ::std::collections::HashMap<prost::alloc::string::String, Constant>,
    /// Child orientations
    #[prost(message, repeated, tag = "6")]
    pub children: prost::alloc::vec::Vec<Orientation>,
}
/// An AttitudeRegistry contains a quickly searchable registry attitude records.
///
/// An AttitudeRegistry should be used to communicate historical data, e.g.
/// instrument X had attitude A at time T with angular velocity W and covariance
/// C. It may also be used to communicate a future desired attitude state, e.g.
/// instrument X should have attitude A at future time T with angular velocity W
/// and within accuracy covariance C all while tracking (or maintaining) the same
/// attitude / angular velocity / both.
#[derive(Clone, PartialEq, prost::Message)]
pub struct AttitudeRegistry {
    /// A pre-sorted list of all of the times (in seconds) available in the map of
    /// states. These time entries are seconds past the start_epoch provided in the
    /// higher struct. Perform a binary search in this index to retrieve time key
    /// for the desired time. In other words, search for the closest time to the
    /// desired time, and retrive the Attitude for this time. Check the repr enum
    /// to understand the attitude representation. If it isn't set, check the
    /// comments or discuss with the publisher of the file.
    ///
    /// NOTE: This provides an O(log(n)) + O(1) access to all of the states. The
    /// O(log(n)) corresponds to the binary search in the index, which then leads
    /// to an O(1) access.
    ///
    /// NOTE: Limitations of protobufs require this index to be an integer.
    /// NOTE: For better platform support, these reference times are limited to 32
    /// bits.
    #[prost(uint32, repeated, tag = "1")]
    pub time_index: prost::alloc::vec::Vec<u32>,
    /// A map associating each time (in seconds) from the index with an Attitude.
    #[prost(map = "uint32, message", tag = "2")]
    pub states: ::std::collections::HashMap<u32, Attitude>,
}
/// Nested message and enum types in `AttitudeRegistry`.
pub mod attitude_registry {
    #[derive(Clone, PartialEq, prost::Message)]
    pub struct AttitudeState {
        /// Absolute epoch
        #[prost(message, optional, tag = "1")]
        pub epoch: ::core::option::Option<super::Epoch>,
        /// The attitude itself
        #[prost(message, optional, tag = "2")]
        pub attitude: ::core::option::Option<super::Attitude>,
        /// The angular velocity, may be unspecified.
        #[prost(message, optional, tag = "3")]
        pub velocity: ::core::option::Option<super::Vector>,
        /// A disambiguation flag specifying whether the velocity is defined.
        /// Defaults to false (i.e vector is indeed defined and should be used).
        #[prost(bool, tag = "4")]
        pub velocity_is_zero: bool,
        /// The covariance array should form a matrix whose diagonal size is either 3
        /// or 4 elements of attitude (depending on its representation) plus 3
        /// elements corresponding to the angular velocity.
        #[prost(double, repeated, tag = "5")]
        pub covariance: prost::alloc::vec::Vec<f64>,
        /// The covariance exponent specifies an optional exponent for all of the
        /// components of the covariance. This enables storing high precision
        /// covariance while not losing precision of floating point values.
        #[prost(double, tag = "6")]
        pub covariance_exponent: f64,
        /// Tracking information to specify whether the spacecraft should be tracking
        /// the quaternion, the angular velocity, or both.
        #[prost(enumeration = "Tracking", tag = "7")]
        pub track: i32,
        /// An optional map of parameters associated to this state
        #[prost(map = "string, message", tag = "8")]
        pub constants: ::std::collections::HashMap<prost::alloc::string::String, super::Constant>,
    }
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, prost::Enumeration)]
    #[repr(i32)]
    pub enum Tracking {
        /// Defaults to None, which should also be used for
        /// historical data (to communicate the attitude of an instrument when a
        /// measurement was taken for example).
        None = 0,
        TrkAttitude = 1,
        TrkVelocity = 2,
        Both = 3,
    }
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Attitude {
    /// The attitude representation (defaults to unspecified with None as a
    /// representation).
    #[prost(enumeration = "AttitudeRepr", tag = "1")]
    pub repr: i32,
    // Note: although attitude and equations can BOTH be set, only one _should_ be
    // set. If both are set, only the attitude will be used, not the equations.
    // A limitation of protobufs prevents the structure from enforcing either/or
    // in the case of a repeated field.
    /// Exact attitude as an array whose size depends on the
    /// representation (e.g. 4 for a quaternion, 3 for an MRP).
    #[prost(double, repeated, tag = "2")]
    pub attitude: prost::alloc::vec::Vec<f64>,
    /// A list of equations to compute the attitude whose default context include
    /// `T` as centuries pas J2000 TDB and `d` as days past J2000 TDB.
    #[prost(message, repeated, tag = "3")]
    pub equations: prost::alloc::vec::Vec<Equation>,
}
/// An AttitudeInterp contains a set of a continuous attitude information.
///
/// It should be used to communicate a specific attitude and angular velocity a
/// spacecraft should maintain between two different times. It may be used to
/// store attitude estimates for historical data as well.
#[derive(Clone, PartialEq, prost::Message)]
pub struct AttitudeInterp {
    /// Type of interpolation used
    #[prost(enumeration = "InterpType", tag = "1")]
    pub itype: i32,
    /// The attitude representation (defaults to unspecified with None as a
    /// representation).
    #[prost(enumeration = "AttitudeRepr", tag = "2")]
    pub repr: i32,
    /// Degree of the interpolation used for computing the attitude (e.g. Piecewise
    /// Linear would have a degree 1, but a Hermite interpolation would usually
    /// have 2*nval - 1 where nval corresponds to the number of states used to
    /// compute the interpolation coefficients).
    #[prost(uint32, tag = "3")]
    pub attitude_degree: u32,
    /// Degree of the interpolation used for computing the velocity. Only used if
    /// the interpolation includes the velocity data.
    #[prost(uint32, tag = "4")]
    pub angular_velocity_degree: u32,
    /// A pre-sorted list of all of the times (in seconds) available in the map of
    /// interpolated states. These time entries are seconds past the
    /// start_mod_julian dates (which is in days). Perform a binary search in this
    /// index to retrieve time key for the desired time. In other words, search for
    /// the closest time to the desired time, retrive the InterpState for this
    /// time, build the interpolation functions, and finally apply these at the
    /// desired time.
    ///
    /// NOTE: This provides an O(log(n)) + O(1) access to all of the states. The
    /// O(log(n)) corresponds to the binary search in the index, which then leads
    /// to an O(1) access. NOTE: Only variable size (or "unequally spaced") windows
    /// are supported. Attitude information is usually provided for short-enough
    /// periods of time that equally spaced interpolations do not add significant
    /// value. One may always use these states as equally spaced windows. NOTE:
    /// Limitations of protobufs require this index to be an integer. NOTE: For
    /// better platform support, these reference times are limited to 32 bits.
    #[prost(uint32, repeated, tag = "5")]
    pub time_index_s: prost::alloc::vec::Vec<u32>,
    /// A map associating each time (in seconds) from the index with an
    /// InterpState.
    #[prost(map = "uint32, message", tag = "6")]
    pub states: ::std::collections::HashMap<u32, attitude_interp::InterpState>,
    /// The minimum value of the time input to the interpolation function (usually
    /// -1 for Chebyshev and Hermite interpolations)
    #[prost(double, tag = "7")]
    pub time_normalization_min: f64,
    /// The maxmum value of the time input to the interpolation function (usually
    /// +1 for Chebyshev and Hermite interpolations)
    #[prost(double, tag = "8")]
    pub time_normalization_max: f64,
}
/// Nested message and enum types in `AttitudeInterp`.
pub mod attitude_interp {
    /// These coefficients should be computed through an interpolation where the
    /// time data is aligned between 0 and 1, unless noted otherwise.
    #[derive(Clone, PartialEq, prost::Message)]
    pub struct InterpState {
        /// Relative time in seconds compared to the indexed time.
        #[prost(double, tag = "1")]
        pub time_offset_s: f64,
        /// Duration in seconds for which these states are valid.
        #[prost(float, tag = "2")]
        pub window_duration: f32,
        /// Unit of the window duration. If unset, assume seconds!
        #[prost(enumeration = "super::Unit", tag = "3")]
        pub time_unit: i32,
        /// All angular velocity coefficients for this time offset.
        #[prost(message, optional, tag = "6")]
        pub angular_velocity: ::core::option::Option<super::VectorCoefficients>,
        #[prost(oneof = "interp_state::Attitude", tags = "4, 5")]
        pub attitude: ::core::option::Option<interp_state::Attitude>,
    }
    /// Nested message and enum types in `InterpState`.
    pub mod interp_state {
        #[derive(Clone, PartialEq, prost::Oneof)]
        pub enum Attitude {
            /// Attitude interpolation can be as quaternion
            #[prost(message, tag = "4")]
            Dim4(super::super::QuaternionCoefficients),
            /// Or as dimension 3 (CRP, MRP, Euler angles, etc.)
            #[prost(message, tag = "5")]
            Dim3(super::super::VectorCoefficients),
        }
    }
}
/// Used to store quaternion coefficients to interpolate.
#[derive(Clone, PartialEq, prost::Message)]
pub struct QuaternionCoefficients {
    #[prost(double, repeated, tag = "1")]
    pub q0: prost::alloc::vec::Vec<f64>,
    #[prost(double, repeated, tag = "2")]
    pub q1: prost::alloc::vec::Vec<f64>,
    #[prost(double, repeated, tag = "3")]
    pub q2: prost::alloc::vec::Vec<f64>,
    #[prost(double, repeated, tag = "4")]
    pub q3: prost::alloc::vec::Vec<f64>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, prost::Enumeration)]
#[repr(i32)]
pub enum AttitudeRepr {
    NotApplicable = 0,
    /// Quaternion where the leading item is the scalar
    QuaternionLead = 1,
    /// Quaternion where the trailing item is the scalar
    QuaternionTrail = 2,
    /// Euler representations
    Euler121 = 3,
    Euler123 = 4,
    Euler131 = 5,
    Euler132 = 6,
    Euler212 = 7,
    Euler213 = 8,
    Euler231 = 9,
    Euler232 = 10,
    Euler312 = 11,
    Euler313 = 12,
    Euler321 = 13,
    Euler323 = 14,
    /// Classical Rodriguez parameters
    Crp = 15,
    /// Modified Rodriguez parameters
    Mrp = 16,
    /// Right ascension / declination / prime meridian location (W)
    RaDecW = 17,
    /// Representation should be specified in the file or ICD
    CustomAttRepr = 18,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Ephemeris {
    /// Nmae of this ephemeris
    #[prost(string, tag = "1")]
    pub name: prost::alloc::string::String,
    /// Name of the orientation frame
    #[prost(string, tag = "2")]
    pub orientation: prost::alloc::string::String,
    #[prost(message, optional, tag = "3")]
    pub start_epoch: ::core::option::Option<Epoch>,
    /// A registry of states without interpolation.
    #[prost(message, optional, tag = "4")]
    pub records: ::core::option::Option<EphemRegistry>,
    #[prost(message, optional, tag = "5")]
    pub interpolator: ::core::option::Option<EphemInterp>,
    /// A map of constant name to constant value and unit.
    #[prost(map = "string, message", tag = "6")]
    pub constants: ::std::collections::HashMap<prost::alloc::string::String, Constant>,
    /// A map of gravity fields available
    #[prost(map = "string, message", tag = "7")]
    pub harmonics: ::std::collections::HashMap<prost::alloc::string::String, GravityField>,
    /// Child ephemerides, ordered from closest to the parent to furthest
    #[prost(message, repeated, tag = "8")]
    pub children: prost::alloc::vec::Vec<Ephemeris>,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct EphemRegistry {
    #[prost(message, repeated, tag = "1")]
    pub states: prost::alloc::vec::Vec<State>,
    #[prost(enumeration = "Unit", tag = "2")]
    pub distance_unit: i32,
    #[prost(enumeration = "Unit", tag = "3")]
    pub velocity_unit: i32,
    /// Frame of these states, useful for reusability
    #[prost(message, optional, tag = "4")]
    pub frame: ::core::option::Option<Frame>,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct State {
    /// Absolute epoch
    #[prost(message, optional, tag = "1")]
    pub epoch: ::core::option::Option<Epoch>,
    /// The position may be unset, check the is_zero flag of the Vector.
    #[prost(message, optional, tag = "2")]
    pub position: ::core::option::Option<Vector>,
    /// The velocity may be unset, check the is_zero flag of the Vector.
    #[prost(message, optional, tag = "3")]
    pub velocity: ::core::option::Option<Vector>,
    /// If the covariance is composed entirely of zeros, it is not informative,
    /// and therefore can be assumed to be unset.
    #[prost(message, optional, tag = "4")]
    pub covariance: ::core::option::Option<state::Covariance>,
    /// The covariance exponent specifies an optional exponent for all of the
    /// components of the covariance. This enables storing high precision
    /// covariance while not losing precision of floating point values.
    #[prost(double, tag = "5")]
    pub covariance_exponent: f64,
    /// A map of constants associated with this state
    #[prost(map = "string, message", tag = "6")]
    pub constants: ::std::collections::HashMap<prost::alloc::string::String, Constant>,
}
/// Nested message and enum types in `State`.
pub mod state {
    #[derive(Clone, PartialEq, prost::Message)]
    pub struct Covariance {
        /// The Covariance of the object is based on the CCSDS OEM Data format.
        ///
        /// Covariance matrix [1,1]
        #[prost(double, tag = "1")]
        pub cx_x: f64,
        /// Covariance matrix [2,1]
        #[prost(double, tag = "2")]
        pub cy_x: f64,
        /// Covariance matrix [2,2]
        #[prost(double, tag = "3")]
        pub cy_y: f64,
        /// Covariance matrix [3,1]
        #[prost(double, tag = "4")]
        pub cz_x: f64,
        /// Covariance matrix [3,2]
        #[prost(double, tag = "5")]
        pub cz_y: f64,
        /// Covariance matrix [3,3]
        #[prost(double, tag = "6")]
        pub cz_z: f64,
        /// Covariance matrix [4,1]
        #[prost(double, tag = "7")]
        pub cx_dot_x: f64,
        /// Covariance matrix [4,2]
        #[prost(double, tag = "8")]
        pub cx_dot_y: f64,
        /// Covariance matrix [4,3]
        #[prost(double, tag = "9")]
        pub cx_dot_z: f64,
        /// Covariance matrix [4,4]
        #[prost(double, tag = "10")]
        pub cx_dot_x_dot: f64,
        /// Covariance matrix [5,1]
        #[prost(double, tag = "11")]
        pub cy_dot_x: f64,
        /// Covariance matrix [5,2]
        #[prost(double, tag = "12")]
        pub cy_dot_y: f64,
        /// Covariance matrix [5,3]
        #[prost(double, tag = "13")]
        pub cy_dot_z: f64,
        /// Covariance matrix [5,4]
        #[prost(double, tag = "14")]
        pub cy_dot_x_dot: f64,
        /// Covariance matrix [5,5]
        #[prost(double, tag = "15")]
        pub cy_dot_y_dot: f64,
        /// Covariance matrix [6,1]
        #[prost(double, tag = "16")]
        pub cz_dot_x: f64,
        /// Covariance matrix [6,2]
        #[prost(double, tag = "17")]
        pub cz_dot_y: f64,
        /// Covariance matrix [6,3]
        #[prost(double, tag = "18")]
        pub cz_dot_z: f64,
        /// Covariance matrix [6,4]
        #[prost(double, tag = "19")]
        pub cz_dot_x_dot: f64,
        /// Covariance matrix [6,5]
        #[prost(double, tag = "20")]
        pub cz_dot_y_dot: f64,
        /// Covariance matrix [6,6]
        #[prost(double, tag = "21")]
        pub cz_dot_z_dot: f64,
    }
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct EphemInterp {
    /// Type of interpolation used
    #[prost(enumeration = "InterpType", tag = "1")]
    pub itype: i32,
    /// Degree of the interpolation used for computing the position (e.g. Piecewise
    /// Linear would have a degree 1, but a Hermite interpolation would usually
    /// have 2*nval - 1 where nval corresponds to the number of states used to
    /// compute the interpolation coefficients).
    #[prost(uint32, tag = "2")]
    pub position_degree: u32,
    /// Degree of the interpolation used for computing the velocity. Only used if
    /// the interpolation includes the velocity data.
    #[prost(uint32, tag = "3")]
    pub velocity_degree: u32,
    #[prost(enumeration = "Unit", tag = "6")]
    pub distance_unit: i32,
    #[prost(enumeration = "Unit", tag = "7")]
    pub velocity_unit: i32,
    #[prost(oneof = "ephem_interp::StateData", tags = "4, 5")]
    pub state_data: ::core::option::Option<ephem_interp::StateData>,
}
/// Nested message and enum types in `EphemInterp`.
pub mod ephem_interp {
    #[derive(Clone, PartialEq, prost::Oneof)]
    pub enum StateData {
        #[prost(message, tag = "4")]
        EqualStates(super::EqualStepStates),
        #[prost(message, tag = "5")]
        VarwindowStates(super::VarWindowStates),
    }
}
/// EqualStepStates provides an O(1) access to all of the states.
/// To access state of time T, get the index as
/// floor((t_in_mod_julian - epoch_mod_julian)/window_duration) .
/// This state's coefficients can be cached for quicker continuous computation.
/// Note that we store the position and velocity outside of a message for smaller
/// serialized structure (two contiguous lists of structures).
#[derive(Clone, PartialEq, prost::Message)]
pub struct EqualStepStates {
    /// Fixed window duration for all of the states (unit specified below, or
    /// assume days)
    #[prost(double, tag = "1")]
    pub window_duration: f64,
    /// Unit of the window duration (is typically days)
    #[prost(enumeration = "Unit", tag = "2")]
    pub window_duration_unit: i32,
    /// All position coefficients for this time offset.
    #[prost(message, repeated, tag = "3")]
    pub position: prost::alloc::vec::Vec<VectorCoefficients>,
    /// All velocity coefficients for this time offset. Optional, but if used, it
    /// **must** be of the same length as the list of position coefficients.
    #[prost(message, repeated, tag = "4")]
    pub velocity: prost::alloc::vec::Vec<VectorCoefficients>,
}
/// VarWindowStates provides an O(log(n)) + O(1) access to all of the states. The
/// O(log(n))
/// corresponds to the binary search in the index, which then leads to an O(1)
/// access.
#[derive(Clone, PartialEq, prost::Message)]
pub struct VarWindowStates {
    /// A pre-sorted list of all of the times (in seconds) available in the map of
    /// interpolated states. These time entries are seconds past the start_epoch
    /// (defined in the parent Ephemeris object). Perform a binary search in this
    /// index to retrieve time key for the desired time. In other words, search for
    /// the closest time to the desired time, retrive the InterpState for this
    /// time, build the interpolation functions, and finally apply these at the
    /// desired time. NOTE: Limitations of protobufs require this index to be an
    /// integer. NOTE: For better platform support, these reference times are
    /// limited to 32 bits.
    #[prost(uint32, repeated, tag = "1")]
    pub time_index_s: prost::alloc::vec::Vec<u32>,
    /// A map associating each time (in seconds) from the index with a InterpState.
    #[prost(map = "uint32, message", tag = "2")]
    pub interp_states: ::std::collections::HashMap<u32, var_window_states::InterpState>,
    /// The minimum value of the time input to the interpolation function (usually
    /// -1 for Chebyshev and Hermite interpolations)
    #[prost(double, tag = "3")]
    pub time_normalization_min: f64,
    /// The maxmum value of the time input to the interpolation function (usually
    /// +1 for Chebyshev and Hermite interpolations)
    #[prost(double, tag = "4")]
    pub time_normalization_max: f64,
}
/// Nested message and enum types in `VarWindowStates`.
pub mod var_window_states {
    /// These coefficients should be computed through an interpolation where the
    /// time data is aligned between 0 and 1, unless noted otherwise.
    #[derive(Clone, PartialEq, prost::Message)]
    pub struct InterpState {
        /// Relative time compared to the indexed time (in seconds). Add this to the
        /// indexed time to compute the relative epoch of the start of this window
        #[prost(double, tag = "1")]
        pub time_offset_s: f64,
        /// Duration of the window: used to compute the inverse interpolant
        #[prost(double, tag = "2")]
        pub window_duration: f64,
        /// Unit of the window duration. If unset, assume seconds!
        #[prost(enumeration = "super::Unit", tag = "3")]
        pub time_unit: i32,
        /// All position coefficients for this time offset.
        #[prost(message, optional, tag = "4")]
        pub position: ::core::option::Option<super::VectorCoefficients>,
        /// All velocity coefficients for this time offset.
        #[prost(message, optional, tag = "5")]
        pub velocity: ::core::option::Option<super::VectorCoefficients>,
    }
}
/// Structure of the gravity field is based on information from
/// https://earth-info.nga.mil/GandG/wgs84/gravitymod/egm2008/README_WGS84_2.pdf
#[derive(Clone, PartialEq, prost::Message)]
pub struct GravityField {
    /// Frame in which this gravity field applies, with center and orientation
    #[prost(message, optional, tag = "1")]
    pub frame: ::core::option::Option<Frame>,
    /// For example, the Earth gravity fields typically start at degree 2
    #[prost(uint32, tag = "2")]
    pub min_degree: u32,
    /// For example, as 70x70 gravity field would set the max_degree to 70
    #[prost(uint32, tag = "3")]
    pub max_degree: u32,
    /// The order always starts at zero, hence there is no min_order
    #[prost(uint32, tag = "4")]
    pub max_order: u32,
    /// The spherical harmonic coefficients. This array must be rebuilt as a
    /// matrix using the degree and order defined above.
    #[prost(double, repeated, tag = "5")]
    pub c_nm: prost::alloc::vec::Vec<f64>,
    #[prost(double, repeated, tag = "6")]
    pub s_nm: prost::alloc::vec::Vec<f64>,
    /// The associated error standard deviations. This array must be rebuilt as a
    /// matrix using the degree and order defined above.
    #[prost(double, repeated, tag = "7")]
    pub sigma_c_nm: prost::alloc::vec::Vec<f64>,
    #[prost(double, repeated, tag = "8")]
    pub sigma_s_nm: prost::alloc::vec::Vec<f64>,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Instrument {
    /// Frame should contain the full path to this instrument
    #[prost(message, optional, tag = "1")]
    pub frame: ::core::option::Option<Frame>,
    /// The position offset from the parent frame
    #[prost(message, optional, tag = "2")]
    pub position_offset: ::core::option::Option<Vector>,
    /// The orientation offset from the parent frame
    #[prost(message, optional, tag = "3")]
    pub orientation_offset: ::core::option::Option<Attitude>,
    #[prost(map = "string, message", tag = "4")]
    pub constants: ::std::collections::HashMap<prost::alloc::string::String, Constant>,
    /// Child instruments
    #[prost(message, repeated, tag = "5")]
    pub children: prost::alloc::vec::Vec<Instrument>,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct NavigationObject {
    /// Identifier of the object
    #[prost(string, tag = "1")]
    pub object: prost::alloc::string::String,
    /// Observer frame which include its center and orientation (e.g. "DSS
    /// 65" and "IAU Earth").
    #[prost(message, optional, tag = "2")]
    pub observer: ::core::option::Option<Frame>,
    /// A list of tracking passes for this object
    #[prost(message, repeated, tag = "3")]
    pub passes: prost::alloc::vec::Vec<TrackingPass>,
    /// A map of constants of that object applicable to all tracking passes.
    #[prost(map = "string, message", tag = "5")]
    pub constants: ::std::collections::HashMap<prost::alloc::string::String, Constant>,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct TrackingPass {
    /// Unique identifier of the navigation data
    #[prost(string, tag = "1")]
    pub pass_uid: prost::alloc::string::String,
    /// Start epoch of the tracking pass
    #[prost(message, optional, tag = "2")]
    pub start_epoch: ::core::option::Option<Epoch>,
    /// End epoch of the tracking pass
    #[prost(message, optional, tag = "3")]
    pub end_epoch: ::core::option::Option<Epoch>,
    /// Actual observations
    #[prost(message, repeated, tag = "4")]
    pub obs: prost::alloc::vec::Vec<Observation>,
    /// An optional set of solutions
    #[prost(message, optional, tag = "5")]
    pub sol: ::core::option::Option<NavSolutions>,
    /// An optional map of parameter name to parameter value and unit.
    #[prost(map = "string, message", tag = "6")]
    pub constants: ::std::collections::HashMap<prost::alloc::string::String, Constant>,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Observation {
    /// Absolute epoch
    #[prost(message, optional, tag = "1")]
    pub epoch: ::core::option::Option<Epoch>,
    /// The observation itself
    #[prost(double, tag = "2")]
    pub obs: f64,
    /// Kind of observation
    #[prost(enumeration = "ObsKind", tag = "3")]
    pub kind: i32,
    /// Unit of the observation
    #[prost(enumeration = "Unit", tag = "4")]
    pub unit: i32,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct NavSolutions {
    /// A registry of states without interpolation.
    #[prost(message, optional, tag = "1")]
    pub state_reg: ::core::option::Option<EphemRegistry>,
    /// Stores the residuals. If the filter processes several observations at the
    /// same time (e.g. Range and Range-Rate), then it is recommended that two
    /// entries we added to this list both with the same epoch.
    #[prost(message, repeated, tag = "2")]
    pub residuals: prost::alloc::vec::Vec<Residual>,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Residual {
    /// Absolute epoch
    #[prost(message, optional, tag = "1")]
    pub epoch: ::core::option::Option<Epoch>,
    /// The prefit residual itself
    #[prost(double, tag = "2")]
    pub prefit: f64,
    /// The postfit residual
    #[prost(double, tag = "3")]
    pub postfit: f64,
    /// Kind of observation which resulted in this set of residuals
    #[prost(enumeration = "ObsKind", tag = "4")]
    pub obskind: i32,
    /// Unit of the residual
    #[prost(enumeration = "Unit", tag = "5")]
    pub unit: i32,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, prost::Enumeration)]
#[repr(i32)]
pub enum ObsKind {
    /// In distance units
    Range = 0,
    /// In velocity units
    RangeRate = 1,
    /// In frequency units
    Doppler = 2,
    /// Dimensionless
    RangeUnits = 3,
    /// In angle units
    Angles = 4,
    /// In angular rate units
    AngluarRate = 5,
    /// Valid only for residual
    None = 10,
    Custom = 11,
}
#[derive(Clone, PartialEq, prost::Message)]
pub struct Xb {
    #[prost(message, optional, tag = "1")]
    pub meta: ::core::option::Option<Metadata>,
    #[prost(map = "string, message", tag = "2")]
    pub constants: ::std::collections::HashMap<prost::alloc::string::String, Constant>,
    #[prost(message, optional, tag = "3")]
    pub ephemeris_root: ::core::option::Option<Ephemeris>,
    #[prost(message, optional, tag = "4")]
    pub orientation_root: ::core::option::Option<Orientation>,
    #[prost(message, repeated, tag = "5")]
    pub navigation_objects: prost::alloc::vec::Vec<NavigationObject>,
    #[prost(message, optional, tag = "6")]
    pub instrument_root: ::core::option::Option<Instrument>,
}
