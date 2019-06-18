#[derive(Clone, PartialEq, Message)]
pub struct FrameContainer {
    #[prost(message, optional, tag = "1")]
    pub meta: ::std::option::Option<Metadata>,
    #[prost(map = "string, message", tag = "2")]
    pub parameters: ::std::collections::HashMap<String, Parameter>,
    #[prost(message, repeated, tag = "3")]
    pub frames: ::std::vec::Vec<Frame>,
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
    /// A parameter value, may be used to specify time system properties (e.g. number of leap seconds).
    #[prost(double, tag = "1")]
    pub value: f64,
    #[prost(string, tag = "2")]
    pub unit: String,
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
    pub frame_parameters: ::std::collections::HashMap<String, Parameter>,
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
        pub name: String,
    }
}
