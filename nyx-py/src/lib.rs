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
use anise::almanac::Almanac;
use anise::almanac::metaload::{MetaAlmanac, MetaFile};
use anise::analysis::prelude::{
    Condition, Event, EventArc, EventDetails, EventEdge, OrbitalElement, Plane, VisibilityArc,
    find_arc_intersections,
};
use anise::analysis::python::{
    PyFrameSpec, PyOrthogonalFrame, PyScalarExpr, PyStateSpec, PyVectorExpr,
};
use anise::analysis::report::PyReportScalars;
use anise::astro::Aberration;
use anise::astro::orbit::Orbit;
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

use nyx_space::cosmic::eclipse::ShadowModel;
use nyx_space::dynamics::guidance::Thruster;
use nyx_space::dynamics::sequence::{
    AccelModels, Dynamics, ForceModels, PropagatorConfig, SpacecraftSequence,
};
use nyx_space::dynamics::{AtmDensity, Drag, PointMasses, SolarPressure};
use nyx_space::io::gravity::GravityFieldConfig;
use nyx_space::mc::{MvnSpacecraft, StateDispersion};
use nyx_space::md::StateParameter;
use nyx_space::md::trajectory::ExportCfg;
use nyx_space::od::GroundStation;
use nyx_space::od::msr::{Measurement, MeasurementType, TrackingDataArc};
use nyx_space::od::noise::link_specific::{CN0, CarrierFreq, ChipRate, SN0};
use nyx_space::od::noise::{GaussMarkov, StochasticNoise, WhiteNoise};
use nyx_space::od::simulator::{Handoff, PyCadence, Scheduler, Strand, TrkConfig};
use nyx_space::propagators::{IntegratorMethod, IntegratorOptions};
use nyx_space::{Spacecraft, cosmic::GuidanceMode};

use pyo3::{prelude::*, wrap_pymodule};

use crate::py_md::Propagator;

mod constants;
mod py_md;
mod py_od;
mod utils;

#[pymodule]
#[pyo3(name = "_nyx")]
fn nyx(m: &Bound<'_, PyModule>) -> PyResult<()> {
    pyo3_log::init();
    m.add_wrapped(wrap_pymodule!(pyanise))?;
    m.add_wrapped(wrap_pymodule!(mission_design))?;
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
    m.add_class::<ExportCfg>()?;
    Ok(())
}

#[pymodule]
fn orbit_determination(_py: Python, sm: &Bound<PyModule>) -> PyResult<()> {
    sm.add_class::<TrackingDataArc>()?;
    sm.add_class::<MeasurementType>()?;
    sm.add_class::<Measurement>()?;
    sm.add_class::<Location>()?;
    sm.add_class::<GroundStation>()?;
    sm.add_class::<PyCadence>()?;
    sm.add_class::<Handoff>()?;
    sm.add_class::<Strand>()?;
    sm.add_class::<TrkConfig>()?;
    sm.add_class::<Scheduler>()?;
    sm.add_class::<ExportCfg>()?;
    sm.add_class::<FrameUid>()?;
    sm.add_class::<StochasticNoise>()?;
    sm.add_class::<WhiteNoise>()?;
    sm.add_class::<GaussMarkov>()?;
    sm.add_class::<SN0>()?;
    sm.add_class::<CN0>()?;
    sm.add_class::<ChipRate>()?;
    sm.add_class::<CarrierFreq>()?;
    sm.add_class::<py_od::GroundTrackingArcSim>()?;
    sm.add_class::<py_od::PySpacecraftODProcess>()?;
    sm.add_class::<py_od::PySpacecraftODSolution>()?;
    sm.add_class::<py_od::PyKfEstimate>()?;
    sm.add_class::<py_od::PyResidual>()?;

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
    sm.add_class::<PropagatorConfig>()?;
    sm.add_class::<Propagator>()?;
    sm.add_class::<Dynamics>()?;
    sm.add_class::<IntegratorMethod>()?;
    sm.add_class::<IntegratorOptions>()?;
    sm.add_class::<AccelModels>()?;
    sm.add_class::<ForceModels>()?;
    sm.add_class::<SpacecraftSequence>()?;
    sm.add_class::<GravityFieldConfig>()?;
    sm.add_class::<PointMasses>()?;
    sm.add_class::<ExportCfg>()?;
    sm.add_class::<SolarPressure>()?;
    sm.add_class::<Drag>()?;
    sm.add_class::<AtmDensity>()?;
    sm.add_class::<ShadowModel>()?;

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
