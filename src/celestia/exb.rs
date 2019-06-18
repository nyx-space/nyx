#[derive(Clone, PartialEq, Message)]
pub struct EphemerisContainer {
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
    pub ephemerides: ::std::vec::Vec<Ephemeris>,
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
pub struct Ephemeris {
    /// Unique identifier of the object
    #[prost(message, optional, tag = "1")]
    pub id: ::std::option::Option<ephemeris::Identifier>,
    /// Unique identifier of the reference frame
    #[prost(message, optional, tag = "2")]
    pub ref_frame: ::std::option::Option<ephemeris::Identifier>,
    /// An optional list of known states. This may be used by an implementer to verify their
    /// computation, or to transmit future propagated states with an optional covariance.
    #[prost(message, repeated, tag = "3")]
    pub states: ::std::vec::Vec<State>,
    #[prost(message, optional, tag = "4")]
    pub interpolator: ::std::option::Option<Interpolation>,
    /// An optional map of parameter name to parameter value and unit.
    #[prost(map = "string, message", tag = "5")]
    pub ephem_parameters: ::std::collections::HashMap<String, Parameter>,
    #[prost(enumeration = "ephemeris::UnitDistance", tag = "6")]
    pub distance_unit: i32,
    #[prost(enumeration = "ephemeris::UnitVelocity", tag = "7")]
    pub velocity_unit: i32,
}
pub mod ephemeris {
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
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
    pub enum UnitDistance {
        Au = 0,
        Km = 1,
        M = 2,
        CustomDistanceUnit = 3,
    }
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
    pub enum UnitVelocity {
        KmS = 0,
        MS = 1,
        CustomVelocityUnit = 2,
    }
}
#[derive(Clone, PartialEq, Message)]
pub struct State {
    #[prost(double, tag = "1")]
    pub mod_julian: f64,
    #[prost(message, optional, tag = "2")]
    pub position: ::std::option::Option<state::Vector>,
    #[prost(message, optional, tag = "3")]
    pub velocity: ::std::option::Option<state::Vector>,
    #[prost(message, optional, tag = "4")]
    pub covariance: ::std::option::Option<state::Covariance>,
    /// The covariance exponent specifies an optional exponent for all of the components of the
    /// covariance. This enables storing high precision covariance while not losing precision of
    /// floating point values.
    #[prost(double, tag = "5")]
    pub covariance_exponent: f64,
}
pub mod state {
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
#[derive(Clone, PartialEq, Message)]
pub struct Interpolation {
    /// Type of interpolation used
    #[prost(enumeration = "interpolation::InterpType", tag = "1")]
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
    /// A start time of all of the states.
    #[prost(double, tag = "4")]
    pub start_mod_julian: f64,
    #[prost(oneof = "interpolation::StateData", tags = "5, 6")]
    pub state_data: ::std::option::Option<interpolation::StateData>,
}
pub mod interpolation {
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
    pub enum InterpType {
        /// Supported interpolations.
        Chebyshev = 0,
        Hermite = 1,
        Polynomial = 2,
        Lagrange = 3,
        /// NOTE: Requires additional communication between provider and implementer.
        CustomInterpolationType = 4,
    }
    #[derive(Clone, Oneof, PartialEq)]
    pub enum StateData {
        #[prost(message, tag = "5")]
        EqualStates(super::EqualStepStates),
        #[prost(message, tag = "6")]
        VarwindowStates(super::VarWindowStates),
    }
}
///
/// These coefficients should be computed through an interpolation where the time data is aligned
/// between 0 and 1 for Hermite interpolatinos, unless noted otherwise in the EXB meta data.
#[derive(Clone, PartialEq, Message)]
pub struct Coefficients {
    #[prost(double, repeated, tag = "1")]
    pub x: ::std::vec::Vec<f64>,
    #[prost(double, repeated, tag = "2")]
    pub y: ::std::vec::Vec<f64>,
    #[prost(double, repeated, tag = "3")]
    pub z: ::std::vec::Vec<f64>,
}
/// EqualStepStates provides an O(1) access to all of the states.
/// To access state of time T, get the index as
/// floor((t_in_mod_julian - start_mod_julian)/window_duration) .
/// This state's coefficients can be cached for quicker continuous computation.
/// Note that we store the position and velocity outside of a message for smaller serialized
/// structure (two contiguous lists of structures).
#[derive(Clone, PartialEq, Message)]
pub struct EqualStepStates {
    /// Fixed window duration for all of the states
    #[prost(double, tag = "1")]
    pub window_duration: f64,
    /// Unit of the window duration
    #[prost(enumeration = "equal_step_states::UnitTime", tag = "2")]
    pub window_duration_unit: i32,
    /// All position coefficients for this time offset.
    #[prost(message, repeated, tag = "3")]
    pub position: ::std::vec::Vec<Coefficients>,
    /// All velocity coefficients for this time offset. Optional, but if used, it **must** be of the
    /// same length as the list of position coefficients.
    #[prost(message, repeated, tag = "4")]
    pub velocity: ::std::vec::Vec<Coefficients>,
}
pub mod equal_step_states {
    #[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, PartialOrd, Ord, Enumeration)]
    pub enum UnitTime {
        Days = 0,
        Seconds = 1,
        CustomTimeUnit = 2,
    }
}
/// VarWindowStates provides an O(log(n)) + O(1) access to all of the states. The O(log(n))
/// corresponds to the binary search in the index, which then leads to an O(1) access.
#[derive(Clone, PartialEq, Message)]
pub struct VarWindowStates {
    /// A pre-sorted list of all of the times (in seconds) available in the map of interpolated states.
    /// These time entries are seconds paste the start_mod_julian dates (which is in days).
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
        pub position: ::std::option::Option<super::Coefficients>,
        /// All velocity coefficients for this time offset. (optional)
        #[prost(message, optional, tag = "4")]
        pub velocity: ::std::option::Option<super::Coefficients>,
    }
}
