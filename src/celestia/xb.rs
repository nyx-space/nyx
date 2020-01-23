#[derive(Clone, PartialEq, Message)]
pub struct Metadata {
    #[prost(message, optional, tag = "1")]
    pub starxb: ::std::option::Option<metadata::Version>,
    #[prost(string, tag = "2")]
    pub publisher: std::string::String,
    #[prost(message, optional, tag = "3")]
    pub date: ::std::option::Option<metadata::CeDate>,
    #[prost(string, tag = "4")]
    pub file_version: std::string::String,
    #[prost(bool, tag = "5")]
    pub proprietary: bool,
    #[prost(string, tag = "6")]
    pub comments: std::string::String,
}
pub mod metadata {
    #[derive(Clone, PartialEq, Message)]
    pub struct CeDate {
        #[prost(uint32, tag = "1")]
        pub year: u32,
        #[prost(uint32, tag = "2")]
        pub month: u32,
        #[prost(uint32, tag = "3")]
        pub day: u32,
    }
    #[derive(Clone, PartialEq, Message)]
    pub struct Version {
        #[prost(uint32, tag = "1")]
        pub major: u32,
        #[prost(uint32, tag = "2")]
        pub minor: u32,
        #[prost(uint32, tag = "3")]
        pub patch: u32,
    }
}
#[derive(Clone, PartialEq, Message)]
pub struct Parameter {
    /// A parameter value, may be used to specify celestial object properties (e.g. GM), or spacecraft
    /// properties (e.g. mass).
    #[prost(double, tag = "1")]
    pub value: f64,
    #[prost(enumeration = "Unit", tag = "2")]
    pub unit: i32,
}
#[derive(Clone, PartialEq, Message)]
pub struct Identifier {
    /// All objects are given an Identifier. Identifiers must be have either a non-zero number
    /// and/or a non-empty string name. The number may be the NAIF_ID, or another identifier
    /// which is understood by the publisher and user of the *XB file. The name may be used to store
    /// the SPACEWARN identifier, or another string understood by both the publisher and the user.
    #[prost(sint32, tag = "1")]
    pub number: i32,
    #[prost(string, tag = "2")]
    pub name: std::string::String,
}
#[derive(Clone, PartialEq, Message)]
pub struct Vector {
    /// The components of the vector, which should be ignored if is_zero is set to true.
    #[prost(double, tag = "1")]
    pub x: f64,
    #[prost(double, tag = "2")]
    pub y: f64,
    #[prost(double, tag = "3")]
    pub z: f64,
    /// A disambiguation flag specifying whether this vector is defined. Defaults to false (i.e vector
    /// is indeed defined and should be used).
    #[prost(bool, tag = "4")]
    pub is_zero: bool,
    #[prost(enumeration = "Unit", tag = "5")]
    pub unit: i32,
}
/// Used to store coefficients to interpolate a three dimensional vector.
#[derive(Clone, PartialEq, Message)]
pub struct VectorCoefficients {
    #[prost(double, repeated, tag = "1")]
    pub x: ::std::vec::Vec<f64>,
    #[prost(double, repeated, tag = "2")]
    pub y: ::std::vec::Vec<f64>,
    #[prost(double, repeated, tag = "3")]
    pub z: ::std::vec::Vec<f64>,
}
/// Defines an Epoch. For example, DE438 uses the ET system with an DaysJDE representation.
#[derive(Clone, PartialEq, Message)]
pub struct Epoch {
    #[prost(enumeration = "TimeSystem", tag = "1")]
    pub ts: i32,
    #[prost(enumeration = "TimeRepr", tag = "2")]
    pub repr: i32,
    #[prost(double, tag = "3")]
    pub value: f64,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
#[repr(i32)]
pub enum Unit {
    Undefined = 0,
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
    /// Angular velocity units
    DegS = 20,
    RadS = 21,
    /// Acceleration units
    KmS2 = 30,
    MS2 = 31,
    CmS2 = 32,
    MmS2 = 33,
    /// Angular acceleration units
    DegS2 = 40,
    RadS2 = 41,
    /// Other units
    ///
    /// Used for graviational parameters
    Km3S2 = 90,
    N = 91,
    CustomUnit = 100,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
#[repr(i32)]
pub enum AttitudeRepr {
    /// Default is that the Attitude state does not have any attitude
    None = 0,
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
    /// Representation should be specified in the file or ICD
    CustomAttRepr = 17,
}
/// Allows specifying the time system used
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
#[repr(i32)]
pub enum TimeSystem {
    /// Temps Atomique International
    Tai = 0,
    /// Ephemeris Time or Barycentric Dynamical Time (IERS TDB)
    Et = 1,
    /// Terrestrial Time
    Tt = 2,
    /// Universal Coordinated Time (assumes ALL leap seconds)
    Utc = 3,
    CustomTs = 4,
}
/// Allow specifying the time representation used
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
#[repr(i32)]
pub enum TimeRepr {
    /// Seconds past 01 January 1900 at MIDNIGHT (or TAI Epoch)
    SecondsJ1900 = 0,
    /// Days of 86400.0 seconds past of the Epoch of the time system
    DaysJ1900 = 1,
    /// Seconds past 01 January 2000 at NOON (NOT midnight, hence different notation)
    SecondsJ2k = 2,
    /// Days of 86400.0 seconds past 01 January 2000 at NOON (NOT midnight, hence different notation)
    DaysJ2k = 3,
    /// Days past the Modified Julian Date Epoch defined as 17 November 1858 at zero hours.
    DaysMjd = 4,
    /// Days past the Julian Date Epoch.
    DaysJde = 5,
    CustomTimeRepr = 6,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
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
#[derive(Clone, PartialEq, Message)]
pub struct EphemerisContainer {
    #[prost(message, optional, tag = "1")]
    pub meta: ::std::option::Option<Metadata>,
    #[prost(map = "string, message", tag = "2")]
    pub parameters: ::std::collections::HashMap<std::string::String, Parameter>,
    #[prost(message, repeated, tag = "3")]
    pub ephemerides: ::std::vec::Vec<Ephemeris>,
}
#[derive(Clone, PartialEq, Message)]
pub struct Ephemeris {
    /// Unique identifier of the object
    #[prost(message, optional, tag = "1")]
    pub id: ::std::option::Option<Identifier>,
    /// Unique identifier of the reference frame
    #[prost(message, optional, tag = "2")]
    pub ref_frame: ::std::option::Option<Identifier>,
    #[prost(message, optional, tag = "3")]
    pub start_epoch: ::std::option::Option<Epoch>,
    /// An optional registry of states without interpolation.
    #[prost(message, optional, tag = "4")]
    pub records: ::std::option::Option<EphemRegistry>,
    #[prost(message, optional, tag = "5")]
    pub interpolator: ::std::option::Option<EphemInterp>,
    /// An optional map of parameter name to parameter value and unit.
    #[prost(map = "string, message", tag = "6")]
    pub parameters: ::std::collections::HashMap<std::string::String, Parameter>,
    #[prost(enumeration = "Unit", tag = "7")]
    pub distance_unit: i32,
    #[prost(enumeration = "Unit", tag = "8")]
    pub velocity_unit: i32,
    #[prost(enumeration = "Unit", tag = "9")]
    pub acceleration_unit: i32,
}
#[derive(Clone, PartialEq, Message)]
pub struct EphemRegistry {
    /// A pre-sorted list of all of the times (in seconds) available in the map of states.
    /// These time entries are seconds past the start_epoch provided in the higher struct.
    /// Perform a binary search in this index to retrieve time key for the desired time.
    /// In other words, search for the closest time to the desired time, and retrive the State
    /// for this time.
    ///
    /// NOTE: This provides an O(log(n)) + O(1) access to all of the states. The O(log(n))
    /// corresponds to the binary search in the index, which then leads to an O(1) access.
    ///
    /// NOTE: Limitations of protobufs require this index to be an integer.
    /// NOTE: For better platform support, these reference times are limited to 32 bits.
    #[prost(uint32, repeated, tag = "1")]
    pub time_index: ::std::vec::Vec<u32>,
    /// A map associating each time (in seconds) from the index with a State.
    #[prost(map = "uint32, message", tag = "2")]
    pub states: ::std::collections::HashMap<u32, ephem_registry::State>,
}
pub mod ephem_registry {
    #[derive(Clone, PartialEq, Message)]
    pub struct State {
        /// Absolute epoch
        #[prost(message, optional, tag = "1")]
        pub epoch: ::std::option::Option<super::Epoch>,
        /// The position may be unset, check the is_zero flag of the Vector.
        #[prost(message, optional, tag = "2")]
        pub position: ::std::option::Option<super::Vector>,
        /// The velocity may be unset, check the is_zero flag of the Vector.
        #[prost(message, optional, tag = "3")]
        pub velocity: ::std::option::Option<super::Vector>,
        /// If the covariance is composed entirely of zeros, it is not informative, and therefore
        /// can be assumed to be unset.
        #[prost(message, optional, tag = "4")]
        pub covariance: ::std::option::Option<state::Covariance>,
        /// The covariance exponent specifies an optional exponent for all of the components of the
        /// covariance. This enables storing high precision covariance while not losing precision of
        /// floating point values.
        #[prost(double, tag = "5")]
        pub covariance_exponent: f64,
        /// An optional map of parameters associated to this state
        #[prost(map = "string, message", tag = "6")]
        pub parameters: ::std::collections::HashMap<std::string::String, super::Parameter>,
    }
    pub mod state {
        #[derive(Clone, PartialEq, Message)]
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
}
#[derive(Clone, PartialEq, Message)]
pub struct EphemInterp {
    /// Type of interpolation used
    #[prost(enumeration = "InterpType", tag = "1")]
    pub itype: i32,
    /// Degree of the interpolation used for computing the position (e.g. Piecewise Linear would have
    /// a degree 1, but a Hermite interpolation would usually have 2*nval - 1 where nval corresponds
    /// to the number of states used to compute the interpolation coefficients).
    #[prost(uint32, tag = "2")]
    pub position_degree: u32,
    /// Degree of the interpolation used for computing the velocity. Only used if the interpolation
    /// includes the velocity data.
    #[prost(uint32, tag = "3")]
    pub velocity_degree: u32,
    #[prost(oneof = "ephem_interp::StateData", tags = "4, 5")]
    pub state_data: ::std::option::Option<ephem_interp::StateData>,
}
pub mod ephem_interp {
    #[derive(Clone, PartialEq, Oneof)]
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
/// Note that we store the position and velocity outside of a message for smaller serialized
/// structure (two contiguous lists of structures).
#[derive(Clone, PartialEq, Message)]
pub struct EqualStepStates {
    /// Fixed window duration for all of the states (in days)
    #[prost(double, tag = "1")]
    pub window_duration: f64,
    /// All position coefficients for this time offset.
    #[prost(message, repeated, tag = "2")]
    pub position: ::std::vec::Vec<VectorCoefficients>,
    /// All velocity coefficients for this time offset. Optional, but if used, it **must** be of the
    /// same length as the list of position coefficients.
    #[prost(message, repeated, tag = "3")]
    pub velocity: ::std::vec::Vec<VectorCoefficients>,
}
/// VarWindowStates provides an O(log(n)) + O(1) access to all of the states. The O(log(n))
/// corresponds to the binary search in the index, which then leads to an O(1) access.
#[derive(Clone, PartialEq, Message)]
pub struct VarWindowStates {
    /// A pre-sorted list of all of the times (in seconds) available in the map of interpolated states.
    /// These time entries are seconds past the epoch (defined in the higher message).
    /// Perform a binary search in this index to retrieve time key for the desired time.
    /// In other words, search for the closest time to the desired time, retrive the InterpState
    /// for this time, build the interpolation functions, and finally apply these at the desired time.
    /// NOTE: Limitations of protobufs require this index to be an integer.
    /// NOTE: For better platform support, these reference times are limited to 32 bits.
    #[prost(uint32, repeated, tag = "1")]
    pub time_index: ::std::vec::Vec<u32>,
    /// A map associating each time (in seconds) from the index with a InterpState.
    #[prost(map = "uint32, message", tag = "2")]
    pub interp_states: ::std::collections::HashMap<u32, var_window_states::InterpState>,
}
pub mod var_window_states {
    /// These coefficients should be computed through an interpolation where the time data is aligned
    /// between 0 and 1 for Hermite interpolatinos, unless noted otherwise in the EXB meta data.
    #[derive(Clone, PartialEq, Message)]
    pub struct InterpState {
        /// Relative time in seconds compared to the indexed time.
        #[prost(float, tag = "1")]
        pub time_offset: f32,
        /// Duration in seconds for which these states are valid.
        #[prost(float, tag = "2")]
        pub window_duration: f32,
        /// All position coefficients for this time offset.
        #[prost(message, optional, tag = "3")]
        pub position: ::std::option::Option<super::VectorCoefficients>,
        /// All velocity coefficients for this time offset. (optional)
        #[prost(message, optional, tag = "4")]
        pub velocity: ::std::option::Option<super::VectorCoefficients>,
    }
}
#[derive(Clone, PartialEq, Message)]
pub struct AttitudeContainer {
    #[prost(message, optional, tag = "1")]
    pub meta: ::std::option::Option<Metadata>,
    #[prost(map = "string, message", tag = "2")]
    pub parameters: ::std::collections::HashMap<std::string::String, Parameter>,
    #[prost(message, repeated, tag = "3")]
    pub kinematics: ::std::vec::Vec<Kinematic>,
}
#[derive(Clone, PartialEq, Message)]
pub struct Kinematic {
    /// Unique identifier of the object
    #[prost(message, optional, tag = "1")]
    pub id: ::std::option::Option<Identifier>,
    /// Unique identifier of the reference frame
    #[prost(message, optional, tag = "2")]
    pub ref_frame: ::std::option::Option<Identifier>,
    /// A start time of all of the states, discrete or interpolated.
    #[prost(message, optional, tag = "3")]
    pub start_epoch: ::std::option::Option<Epoch>,
    #[prost(message, optional, tag = "4")]
    pub records: ::std::option::Option<AttitudeRegistry>,
    #[prost(message, optional, tag = "5")]
    pub interpolator: ::std::option::Option<AttitudeInterp>,
    /// An optional map of parameter name to parameter value and unit.
    #[prost(map = "string, message", tag = "6")]
    pub parameters: ::std::collections::HashMap<std::string::String, Parameter>,
    #[prost(enumeration = "Unit", tag = "7")]
    pub angular_velocity_unit: i32,
    #[prost(enumeration = "Unit", tag = "8")]
    pub angular_acceleration_unit: i32,
}
/// An AttitudeRegistry contains a quickly searchable registry attitude records.
///
/// An AttitudeRegistry should be used to communicate historical data, e.g. instrument X had attitude A at
/// time T with angular velocity W and covariance C.
/// It may also be used to communicate a future desired attitude state, e.g. instrument
/// X should have attitude A at future time T with angular velocity W and within accuracy covariance
/// C all while tracking (or maintaining) the same attitude / angular velocity / both.
#[derive(Clone, PartialEq, Message)]
pub struct AttitudeRegistry {
    /// A pre-sorted list of all of the times (in seconds) available in the map of states.
    /// These time entries are seconds past the start_epoch provided in the higher struct.
    /// Perform a binary search in this index to retrieve time key for the desired time.
    /// In other words, search for the closest time to the desired time, and retrive the Attitude
    /// for this time.
    /// Check the repr enum to understand the attitude representation. If it isn't set, check the
    /// comments or discuss with the publisher of the file.
    ///
    /// NOTE: This provides an O(log(n)) + O(1) access to all of the states. The O(log(n))
    /// corresponds to the binary search in the index, which then leads to an O(1) access.
    ///
    /// NOTE: Limitations of protobufs require this index to be an integer.
    /// NOTE: For better platform support, these reference times are limited to 32 bits.
    #[prost(uint32, repeated, tag = "1")]
    pub time_index: ::std::vec::Vec<u32>,
    /// A map associating each time (in seconds) from the index with an Attitude.
    #[prost(map = "uint32, message", tag = "2")]
    pub states: ::std::collections::HashMap<u32, attitude_registry::Attitude>,
}
pub mod attitude_registry {
    #[derive(Clone, PartialEq, Message)]
    pub struct Attitude {
        /// Absolute epoch
        #[prost(message, optional, tag = "1")]
        pub epoch: ::std::option::Option<super::Epoch>,
        /// The attitude representation (defaults to unspecified with None as a representation).
        #[prost(enumeration = "super::AttitudeRepr", tag = "2")]
        pub repr: i32,
        /// The exact attitude at this time. If no attitude is provided, this array will be empty.
        /// Attitude representation is left to the freedom of the implementer.
        /// Refer to common.proto for the recommended representations.
        #[prost(double, repeated, tag = "3")]
        pub attitude: ::std::vec::Vec<f64>,
        /// The angular velocity, may be unspecified.
        #[prost(message, optional, tag = "4")]
        pub velocity: ::std::option::Option<super::Vector>,
        /// If the covariance is composed entirely of zeros, it is not informative, and therefore
        /// can be assumed to be unset.
        #[prost(message, optional, tag = "5")]
        pub covariance: ::std::option::Option<attitude::Covariance>,
        /// The covariance exponent specifies an optional exponent for all of the components of the
        /// covariance. This enables storing high precision covariance while not losing precision of
        /// floating point values.
        #[prost(double, tag = "6")]
        pub covariance_exponent: f64,
        /// Tracking information to specify whether the spacecraft should be tracking the
        /// quaternion, the angular velocity, or both.
        #[prost(enumeration = "Tracking", tag = "7")]
        pub track: i32,
        /// An optional map of parameters associated to this state
        #[prost(map = "string, message", tag = "8")]
        pub parameters: ::std::collections::HashMap<std::string::String, super::Parameter>,
    }
    pub mod attitude {
        #[derive(Clone, PartialEq, Message)]
        pub struct Covariance {
            /// The items starting with a `a` correspond to the attitude components.
            /// The items starting with a `w` correspond to the angular velocity vector omega.
            ///
            /// Covariance matrix [1,1]
            #[prost(double, tag = "1")]
            pub ca0_a0: f64,
            /// Covariance matrix [2,1]
            #[prost(double, tag = "2")]
            pub ca1_a0: f64,
            /// Covariance matrix [2,2]
            #[prost(double, tag = "3")]
            pub ca1_a1: f64,
            /// Covariance matrix [3,1]
            #[prost(double, tag = "4")]
            pub ca2_a0: f64,
            /// Covariance matrix [3,2]
            #[prost(double, tag = "5")]
            pub ca2_a1: f64,
            /// Covariance matrix [3,3]
            #[prost(double, tag = "6")]
            pub ca2_a2: f64,
            /// Covariance matrix [4,1]
            #[prost(double, tag = "7")]
            pub ca3_a0: f64,
            /// Covariance matrix [4,2]
            #[prost(double, tag = "8")]
            pub ca3_a1: f64,
            /// Covariance matrix [4,3]
            #[prost(double, tag = "9")]
            pub ca3_a2: f64,
            /// Covariance matrix [4,4]
            #[prost(double, tag = "10")]
            pub ca3_a3: f64,
            /// Covariance matrix [5,1]
            #[prost(double, tag = "11")]
            pub cwx_a0: f64,
            /// Covariance matrix [5,2]
            #[prost(double, tag = "12")]
            pub cwx_a1: f64,
            /// Covariance matrix [5,3]
            #[prost(double, tag = "13")]
            pub cwx_a2: f64,
            /// Covariance matrix [5,4]
            #[prost(double, tag = "14")]
            pub cwx_a3: f64,
            /// Covariance matrix [5,5]
            #[prost(double, tag = "15")]
            pub cwx_wx: f64,
            /// Covariance matrix [6,1]
            #[prost(double, tag = "16")]
            pub cwy_a0: f64,
            /// Covariance matrix [6,2]
            #[prost(double, tag = "17")]
            pub cwy_a1: f64,
            /// Covariance matrix [6,3]
            #[prost(double, tag = "18")]
            pub cwy_a2: f64,
            /// Covariance matrix [6,4]
            #[prost(double, tag = "19")]
            pub cwy_a3: f64,
            /// Covariance matrix [6,5]
            #[prost(double, tag = "20")]
            pub cwy_wx: f64,
            /// Covariance matrix [6,5]
            #[prost(double, tag = "21")]
            pub cwy_wy: f64,
            /// Covariance matrix [7,1]
            #[prost(double, tag = "22")]
            pub cwz_a0: f64,
            /// Covariance matrix [7,2]
            #[prost(double, tag = "23")]
            pub cwz_a1: f64,
            /// Covariance matrix [7,3]
            #[prost(double, tag = "24")]
            pub cwz_a2: f64,
            /// Covariance matrix [7,4]
            #[prost(double, tag = "25")]
            pub cwz_a3: f64,
            /// Covariance matrix [7,5]
            #[prost(double, tag = "26")]
            pub cwz_wx: f64,
            /// Covariance matrix [7,6]
            #[prost(double, tag = "27")]
            pub cwz_wy: f64,
            /// Covariance matrix [7,7]
            #[prost(double, tag = "28")]
            pub cwz_wz: f64,
        }
    }
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
    #[repr(i32)]
    pub enum Tracking {
        /// Defaults to None, which should also be used for
        /// historical data (to communicate the attitude of an instrument when a measurement was taken
        /// for example).
        None = 0,
        TrkAttitude = 1,
        TrkVelocity = 2,
        Both = 3,
    }
}
/// An AttitudeInterp contains a set of a continuous attitude information.
///
/// It should be used to communicate a specific attitude and angular velocity a spacecraft should
/// maintain between two different times. It may be used to store attitude estimates for historical
/// data as well.
#[derive(Clone, PartialEq, Message)]
pub struct AttitudeInterp {
    /// Type of interpolation used
    #[prost(enumeration = "InterpType", tag = "1")]
    pub itype: i32,
    /// The attitude representation (defaults to unspecified with None as a representation).
    #[prost(enumeration = "AttitudeRepr", tag = "2")]
    pub repr: i32,
    /// Degree of the interpolation used for computing the attitude (e.g. Piecewise Linear would have
    /// a degree 1, but a Hermite interpolation would usually have 2*nval - 1 where nval corresponds
    /// to the number of states used to compute the interpolation coefficients).
    #[prost(uint32, tag = "3")]
    pub attitude_degree: u32,
    /// Degree of the interpolation used for computing the velocity. Only used if the interpolation
    /// includes the velocity data.
    #[prost(uint32, tag = "4")]
    pub angular_velocity_degree: u32,
    /// A pre-sorted list of all of the times (in seconds) available in the map of interpolated states.
    /// These time entries are seconds past the start_mod_julian dates (which is in days).
    /// Perform a binary search in this index to retrieve time key for the desired time.
    /// In other words, search for the closest time to the desired time, retrive the InterpState
    /// for this time, build the interpolation functions, and finally apply these at the desired time.
    ///
    /// NOTE: This provides an O(log(n)) + O(1) access to all of the states. The O(log(n))
    /// corresponds to the binary search in the index, which then leads to an O(1) access.
    /// NOTE: Only variable size (or "unequally spaced") windows are supported. Attitude information
    /// is usually provided for short-enough periods of time that equally spaced interpolations do not
    /// add significant value. One may always use these states as equally spaced windows.
    /// NOTE: Limitations of protobufs require this index to be an integer.
    /// NOTE: For better platform support, these reference times are limited to 32 bits.
    #[prost(uint32, repeated, tag = "5")]
    pub time_index: ::std::vec::Vec<u32>,
    /// A map associating each time (in seconds) from the index with an InterpState.
    #[prost(map = "uint32, message", tag = "6")]
    pub states: ::std::collections::HashMap<u32, attitude_interp::InterpState>,
}
pub mod attitude_interp {
    /// These coefficients are used to represent an attitude as a set of interpolation coefficients.
    #[derive(Clone, PartialEq, Message)]
    pub struct AttitudeCoefficients {
        #[prost(double, repeated, tag = "1")]
        pub a0: ::std::vec::Vec<f64>,
        #[prost(double, repeated, tag = "2")]
        pub a1: ::std::vec::Vec<f64>,
        #[prost(double, repeated, tag = "3")]
        pub a2: ::std::vec::Vec<f64>,
        #[prost(double, repeated, tag = "4")]
        pub a3: ::std::vec::Vec<f64>,
    }
    #[derive(Clone, PartialEq, Message)]
    pub struct InterpState {
        /// Relative time in seconds compared to the indexed time.
        #[prost(float, tag = "1")]
        pub time_offset: f32,
        /// Duration in seconds for which these states are valid.
        #[prost(float, tag = "2")]
        pub window_duration: f32,
        /// All attitude coefficients for this time offset.
        #[prost(message, optional, tag = "3")]
        pub attitude: ::std::option::Option<AttitudeCoefficients>,
        /// All angular velocity coefficients for this time offset. (optional)
        #[prost(message, optional, tag = "4")]
        pub angular_velocity: ::std::option::Option<super::VectorCoefficients>,
    }
}
#[derive(Clone, PartialEq, Message)]
pub struct FrameContainer {
    #[prost(message, optional, tag = "1")]
    pub meta: ::std::option::Option<Metadata>,
    #[prost(map = "string, message", tag = "2")]
    pub parameters: ::std::collections::HashMap<std::string::String, Parameter>,
    #[prost(message, repeated, tag = "3")]
    pub frames: ::std::vec::Vec<Frame>,
}
#[derive(Clone, PartialEq, Message)]
pub struct Frame {
    /// Unique identifier of the frame
    #[prost(message, optional, tag = "1")]
    pub id: ::std::option::Option<frame::Identifier>,
    /// Unique identifier of the *XB for the translation component (EXB).
    #[prost(message, optional, tag = "2")]
    pub exb_id: ::std::option::Option<frame::Identifier>,
    /// Unique identifier of the *XB for the rotation component (AXB).
    #[prost(message, optional, tag = "3")]
    pub axb_id: ::std::option::Option<frame::Identifier>,
    /// An optional map of parameter name to parameter value and unit.
    #[prost(map = "string, message", tag = "4")]
    pub parameters: ::std::collections::HashMap<std::string::String, Parameter>,
}
pub mod frame {
    #[derive(Clone, PartialEq, Message)]
    pub struct Identifier {
        /// All objects are given an Identifier. Identifiers must be have either a non-zero number
        /// and/or a non-empty string name. The number may be the NAIF_ID, or another identifier
        /// which is understood by the publisher and user of the *XB file. The name may be used to store
        /// the SPACEWARN identifier, or another string understood by both the publisher and the user.
        #[prost(sint32, tag = "1")]
        pub number: i32,
        #[prost(string, tag = "2")]
        pub name: std::string::String,
    }
}
