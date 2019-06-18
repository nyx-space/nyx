#[derive(Clone, PartialEq, Message)]
pub struct AttitudeContainer {
    #[prost(message, optional, tag = "1")]
    pub meta: ::std::option::Option<Metadata>,
    /// Offset of time references `T` wrt to the Modified Julian Date of "1858 Nov 17 zero hours".
    /// The computation must be as such: `MJD = start_mod_julian + mod_julian_tai_offset`. Hence, if
    /// the reference of the calendar is _after_ the MJD reference, mod_julian_tai_offset is a
    /// positive number. For example, if the time references are in Julian Dates, then
    /// mod_julian_tai_offset = -2400000.50.
    #[prost(double, tag = "2")]
    pub mod_julian_tai_offset: f64,
    #[prost(map = "string, message", tag = "3")]
    pub parameters: ::std::collections::HashMap<String, Parameter>,
    #[prost(message, repeated, tag = "4")]
    pub kinematics: ::std::vec::Vec<Kinematic>,
}
#[derive(Clone, PartialEq, Message)]
pub struct Metadata {
    #[prost(message, optional, tag = "1")]
    pub starxb: ::std::option::Option<metadata::Version>,
    #[prost(string, tag = "2")]
    pub publisher: String,
    #[prost(message, optional, tag = "3")]
    pub date: ::std::option::Option<metadata::CeDate>,
    #[prost(string, tag = "4")]
    pub file_version: String,
    #[prost(bool, tag = "5")]
    pub proprietary: bool,
    #[prost(string, tag = "6")]
    pub comments: String,
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
    #[prost(string, tag = "2")]
    pub unit: String,
}
#[derive(Clone, PartialEq, Message)]
pub struct Kinematic {
    /// Unique identifier of the object
    #[prost(message, optional, tag = "1")]
    pub id: ::std::option::Option<kinematic::Identifier>,
    /// Unique identifier of the reference frame
    #[prost(message, optional, tag = "2")]
    pub ref_frame: ::std::option::Option<kinematic::Identifier>,
    /// A start time of all of the states, discrete or interpolated.
    #[prost(double, tag = "3")]
    pub start_mod_julian: f64,
    #[prost(message, optional, tag = "4")]
    pub schedule: ::std::option::Option<StateMap>,
    #[prost(message, optional, tag = "5")]
    pub interpolator: ::std::option::Option<Interpolation>,
    /// An optional map of parameter name to parameter value and unit.
    #[prost(map = "string, message", tag = "6")]
    pub frame_parameters: ::std::collections::HashMap<String, Parameter>,
    #[prost(string, tag = "7")]
    pub angular_velocity_unit: String,
}
pub mod kinematic {
    #[derive(Clone, PartialEq, Message)]
    pub struct Identifier {
        /// All objects are given an Identifier. Identifiers must be have either a non-zero number
        /// and/or a non-empty string name. The number may be the NAIF_ID, or another identifier
        /// which is understood by the publisher and user of the *XB file. The name may be used to store
        /// the SPACEWARN identifier, or another string understood by both the publisher and the user.
        #[prost(sint32, tag = "1")]
        pub number: i32,
        #[prost(string, tag = "2")]
        pub name: String,
    }
}
/// A StateMap contains a map of times to discrete attitude states.
///
/// A StateMap should be used to communicate historical data, e.g. instrument X had attitude A at
/// time T with angular velocity W and covariance C.
/// It may also be used to communicate a future desired attitude state, e.g. instrument
/// X should have attitude A at future time T with angular velocity W and within accuracy covariance
/// C all while tracking (or maintaining) the same attitude / angular velocity / both.
#[derive(Clone, PartialEq, Message)]
pub struct StateMap {
    /// A pre-sorted list of all of the times (in seconds) available in the map of states.
    /// These time entries are seconds past the start_mod_julian dates (which is in days).
    /// Perform a binary search in this index to retrieve time key for the desired time.
    /// In other words, search for the closest time to the desired time, and retrive the State
    /// for this time.
    /// Check the repr enum to understand the attitude representation. If it isn't set, check the
    /// comments or discuss with the publisher of the file.
    /// Check if the has_attitude flag is set prior to interpreting the attitude struct.
    /// Check if the has_velocity flag is set prior to interpreting the angular_velocity vector.
    ///
    /// NOTE: This provides an O(log(n)) + O(1) access to all of the states. The O(log(n))
    /// corresponds to the binary search in the index, which then leads to an O(1) access.
    /// NOTE: Only variable size (or "unequally spaced") windows are supported. Attitude information
    /// is usually provided for short periods of time so equally spaced interpolations do not
    /// add significant value. One may always use these states as equally spaced windows.
    ///
    /// NOTE: Limitations of protobufs require this index to be an integer.
    /// NOTE: For better platform support, these reference times are limited to 32 bits.
    #[prost(uint32, repeated, tag = "1")]
    pub time_index: ::std::vec::Vec<u32>,
    /// A map associating each time (in seconds) from the index with a State.
    #[prost(map = "uint32, message", tag = "2")]
    pub states: ::std::collections::HashMap<u32, state_map::State>,
}
pub mod state_map {
    #[derive(Clone, PartialEq, Message)]
    pub struct State {
        /// Relative time in seconds compared to the indexed time.
        #[prost(float, tag = "1")]
        pub time_offset: f32,
        /// The attitude representation
        #[prost(enumeration = "state::AttitudeRepresentation", tag = "2")]
        pub repr: i32,
        /// A disambiguation flag specifying whether this discrete state has an attitude.
        /// The struct defaults to zeros which may be a valid attitude. If this flag is set to
        /// false, then the covariance for the attitude will also be zeros.
        #[prost(bool, tag = "3")]
        pub has_attitude: bool,
        /// The exact attitude at this time. If no attitude is provided, this struct will
        /// only contain zeros, which is an invalid quaternion (the norm of a quaternion is always one).
        #[prost(message, optional, tag = "4")]
        pub attitude: ::std::option::Option<state::Attitude>,
        /// A disambiguation flag specifying whether this discrete state has an angular velocity.
        /// The vector defaults to zeros which would be a valid angular velocity. If this flag is set to
        /// false, then the covariance for the angular velocity will also be zeros.
        #[prost(bool, tag = "5")]
        pub has_velocity: bool,
        /// The exact angular velocity at this time. Data in this struct should be ignored if the value
        /// of has_velocity is set to false.
        #[prost(message, optional, tag = "6")]
        pub angular_velocity: ::std::option::Option<state::Vector>,
        /// The exact covariance at this time. Data in this struct should be ignored if the value
        /// of has_velocity is set to false.
        #[prost(message, optional, tag = "7")]
        pub covariance: ::std::option::Option<state::Covariance>,
        /// The covariance exponent specifies an optional exponent for all of the components of the
        /// covariance. This enables storing high precision covariance while not losing precision of
        /// floating point values.
        #[prost(double, tag = "8")]
        pub covariance_exponent: f64,
        /// Optional tracking information to specify whether the spacecraft should be tracking the
        /// quaternion, the angular velocity, or both.
        #[prost(enumeration = "Tracking", tag = "9")]
        pub track: i32,
    }
    pub mod state {
        /// Attitude representation is left to the freedom of the implementer. However, it must fit
        /// within four dimensions.
        /// If using a shorter than 4-dimensional representation of attitude (like MRPs), it is
        /// recommended to specify in the comments section which component is always fixed to zero. A
        /// further recommendation is to set the fourth component to zero in those cases.
        #[derive(Clone, PartialEq, Message)]
        pub struct Attitude {
            #[prost(double, tag = "1")]
            pub a0: f64,
            #[prost(double, tag = "2")]
            pub a1: f64,
            #[prost(double, tag = "3")]
            pub a2: f64,
            #[prost(double, tag = "4")]
            pub a3: f64,
        }
        #[derive(Clone, PartialEq, Message)]
        pub struct Vector {
            #[prost(double, tag = "1")]
            pub x: f64,
            #[prost(double, tag = "2")]
            pub y: f64,
            #[prost(double, tag = "3")]
            pub z: f64,
        }
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
        #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
        pub enum AttitudeRepresentation {
            Quaternion = 0,
            Custom = 1,
        }
    }
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
    pub enum Tracking {
        /// Defaults to None, which should also be used for
        /// historical data (to communicate the attitude of an instrument when a measurement was taken
        /// for example).
        None = 0,
        Attitude = 1,
        Velocity = 2,
        Both = 3,
    }
}
/// An Interpolation contains a set of a continuous attitude information.
///
/// It should be used to communicate a specific attitude and angular velocity a spacecraft should
/// maintain between two different times. It may be used to store attitude estimates for historical
/// data as well.
#[derive(Clone, PartialEq, Message)]
pub struct Interpolation {
    /// Type of interpolation used
    #[prost(enumeration = "interpolation::InterpType", tag = "1")]
    pub itype: i32,
    /// Degree of the interpolation used for computing the quaternion (e.g. Piecewise Linear would have
    /// a degree 1, but a Hermite interpolation would usually have 2*nval - 1 where nval corresponds
    /// to the number of states used to compute the interpolation coefficients).
    #[prost(uint32, tag = "2")]
    pub quaternion_degree: u32,
    /// Degree of the interpolation used for computing the velocity. Only used if the interpolation
    /// includes the velocity data.
    #[prost(uint32, tag = "3")]
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
    #[prost(uint32, repeated, tag = "4")]
    pub time_index: ::std::vec::Vec<u32>,
    /// A map associating each time (in seconds) from the index with an InterpState.
    #[prost(map = "uint32, message", tag = "5")]
    pub states: ::std::collections::HashMap<u32, interpolation::InterpState>,
}
pub mod interpolation {
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
    /// These coefficients are used to represent an angular velocity vector as a set of interpolation
    /// coefficients.
    #[derive(Clone, PartialEq, Message)]
    pub struct VectorCoefficients {
        #[prost(double, repeated, tag = "1")]
        pub x: ::std::vec::Vec<f64>,
        #[prost(double, repeated, tag = "2")]
        pub y: ::std::vec::Vec<f64>,
        #[prost(double, repeated, tag = "3")]
        pub z: ::std::vec::Vec<f64>,
    }
    #[derive(Clone, PartialEq, Message)]
    pub struct InterpState {
        /// Relative time in seconds compared to the indexed time.
        #[prost(float, tag = "1")]
        pub time_offset: f32,
        /// Duration in seconds for which these states are valid.
        #[prost(float, tag = "2")]
        pub window_duration: f32,
        /// All quaternion coefficients for this time offset.
        #[prost(message, optional, tag = "3")]
        pub quaternion: ::std::option::Option<AttitudeCoefficients>,
        /// All angular velocity coefficients for this time offset. (optional)
        #[prost(message, optional, tag = "4")]
        pub angular_velocity: ::std::option::Option<VectorCoefficients>,
    }
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
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
}
