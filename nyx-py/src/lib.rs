use anise::almanac::metaload::{MetaAlmanac, MetaFile};
use anise::almanac::Almanac;
use anise::analysis::prelude::{
    find_arc_intersections, Condition, Event, EventArc, EventDetails, EventEdge, OrbitalElement,
    Plane, VisibilityArc,
};
use anise::analysis::python::{
    PyFrameSpec, PyOrthogonalFrame, PyScalarExpr, PyStateSpec, PyVectorExpr,
};
use anise::analysis::report::PyReportScalars;
use anise::astro::orbit::Orbit;
use anise::astro::Aberration;
use anise::astro::{AzElRange, Location, Occultation, TerrainMask};
use anise::ephemerides::ephemeris::{Covariance, Ephemeris, EphemerisRecord, LocalFrame};
use anise::frames::Frame;
use anise::frames::FrameUid;
use anise::naif::daf::DafDataType;
use anise::structure::dataset::location_dhall::PyLocationDataSet;
use anise::structure::dataset::location_dhall::{LocationDhallSet, LocationDhallSetEntry};
use anise::structure::instrument::{FovShape, Instrument};
use anise::structure::planetocentric::ellipsoid::Ellipsoid;
use anise::structure::spacecraft::{DragData, Mass, SRPData};
use hifitime::leap_seconds::*;
use hifitime::python::*;
use hifitime::ut1::*;
use hifitime::*;
use nyx_space::dynamics::guidance::Thruster;
use nyx_space::dynamics::sequence::SpacecraftSequence;
use nyx_space::mc::{MvnSpacecraft, StateDispersion};
use nyx_space::md::StateParameter;
use nyx_space::{
    cosmic::GuidanceMode,
    od::msr::{Measurement, MeasurementType, TrackingDataArc},
    Spacecraft,
};

use pyo3::{prelude::*, wrap_pymodule};

mod constants;
mod utils;

#[pymodule]
#[pyo3(name = "_nyx")]
fn nyx(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
    m.add_wrapped(wrap_pymodule!(pyanise))?;
    m.add_wrapped(wrap_pymodule!(orbit_determination))?;
    m.add_wrapped(wrap_pymodule!(monte_carlo))?;
    m.add_wrapped(wrap_pymodule!(time))?;

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__doc__", env!("CARGO_PKG_DESCRIPTION"))?;
    m.add("__author__", env!("CARGO_PKG_AUTHORS"))?;

    // Export everything needed to build a spacecraft
    m.add_class::<Spacecraft>()?;
    m.add_class::<GuidanceMode>()?;
    m.add_class::<Mass>()?;
    m.add_class::<SRPData>()?;
    m.add_class::<DragData>()?;
    m.add_class::<Thruster>()?;
    Ok(())
}

#[pymodule]
fn orbit_determination(_py: Python, sm: &Bound<PyModule>) -> PyResult<()> {
    sm.add_class::<TrackingDataArc>()?;
    sm.add_class::<MeasurementType>()?;
    sm.add_class::<Measurement>()?;

    Ok(())
}

#[pymodule]
fn monte_carlo(_py: Python, sm: &Bound<PyModule>) -> PyResult<()> {
    sm.add_class::<MvnSpacecraft>()?;
    sm.add_class::<StateParameter>()?;
    sm.add_class::<StateDispersion>()?;

    Ok(())
}

#[pymodule]
fn mission_design(_py: Python, sm: &Bound<PyModule>) -> PyResult<()> {
    sm.add_class::<SpacecraftSequence>()?;

    Ok(())
}

/// Reexport hifitime as nyx_space.time
#[pymodule]
fn time(_py: Python, sm: &Bound<PyModule>) -> PyResult<()> {
    sm.add_class::<Epoch>()?;
    sm.add_class::<TimeScale>()?;
    sm.add_class::<TimeSeries>()?;
    sm.add_class::<Duration>()?;
    sm.add_class::<Unit>()?;
    sm.add_class::<LatestLeapSeconds>()?;
    sm.add_class::<LeapSecondsFile>()?;
    sm.add_class::<Ut1Provider>()?;
    sm.add_class::<MonthName>()?;
    sm.add_class::<PyHifitimeError>()?;
    sm.add_class::<PyDurationError>()?;
    sm.add_class::<PyParsingError>()?;
    sm.add_class::<Polynomial>()?;
    sm.add_class::<Weekday>()?;
    Ok(())
}

/// Reexport of ANISE as nyx_space.anise

#[pymodule]
#[pyo3(name = "anise")]
fn pyanise(_py: Python, m: &Bound<PyModule>) -> PyResult<()> {
    m.add_wrapped(wrap_pymodule!(time))?;
    m.add_wrapped(wrap_pymodule!(analysis))?;
    m.add_wrapped(wrap_pymodule!(instrument))?;
    m.add_wrapped(wrap_pymodule!(astro))?;
    m.add_wrapped(wrap_pymodule!(constants::constants))?;
    m.add_wrapped(wrap_pymodule!(utils::utils))?;

    m.add_class::<Almanac>()?;
    m.add_class::<Aberration>()?;
    m.add_class::<MetaAlmanac>()?;
    m.add_class::<MetaFile>()?;
    m.add_class::<LocationDhallSet>()?;
    m.add_class::<LocationDhallSetEntry>()?;
    m.add_class::<PyLocationDataSet>()?;

    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add("__doc__", env!("CARGO_PKG_DESCRIPTION"))?;
    m.add("__author__", env!("CARGO_PKG_AUTHORS"))?;

    Ok(())
}

#[pymodule]
fn analysis(_py: Python, sm: &Bound<PyModule>) -> PyResult<()> {
    sm.add_class::<PyStateSpec>()?;
    sm.add_class::<PyOrthogonalFrame>()?;
    sm.add_class::<PyVectorExpr>()?;
    sm.add_class::<PyScalarExpr>()?;
    sm.add_class::<PyFrameSpec>()?;
    sm.add_class::<OrbitalElement>()?;
    sm.add_class::<Condition>()?;
    sm.add_class::<Plane>()?;
    sm.add_class::<Event>()?;
    sm.add_class::<EventDetails>()?;
    sm.add_class::<EventEdge>()?;
    sm.add_class::<EventArc>()?;
    sm.add_class::<VisibilityArc>()?;
    sm.add_class::<PyReportScalars>()?;
    sm.add_wrapped(wrap_pyfunction!(find_arc_intersections))?;
    Ok(())
}

#[pymodule]
fn instrument(_py: Python, sm: &Bound<PyModule>) -> PyResult<()> {
    sm.add_class::<Instrument>()?;
    sm.add_class::<FovShape>()?;
    sm.add_class::<Ellipsoid>()?;
    Ok(())
}

#[pymodule]
pub(crate) fn astro(_py: Python, sm: &Bound<'_, PyModule>) -> PyResult<()> {
    sm.add_class::<Ellipsoid>()?;
    sm.add_class::<Frame>()?;
    sm.add_class::<FrameUid>()?;
    sm.add_class::<Orbit>()?;
    sm.add_class::<AzElRange>()?;
    sm.add_class::<Occultation>()?;
    sm.add_class::<Location>()?;
    sm.add_class::<TerrainMask>()?;
    sm.add_class::<Ephemeris>()?;
    sm.add_class::<EphemerisRecord>()?;
    sm.add_class::<Covariance>()?;
    sm.add_class::<LocalFrame>()?;
    sm.add_class::<DafDataType>()?;

    Ok(())
}
