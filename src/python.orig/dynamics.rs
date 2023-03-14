use std::sync::Arc;
use pyo3::{prelude::*, types::PyType, types::PyDict};
use crate::dynamics::orbital::OrbitalDynamics as OrbitalDynamicsRs;
use crate::dynamics::solarpressure::SolarPressure as SolarPressureRs;
use crate::dynamics::drag::Drag as DragRs;
use crate::dynamics::spacecraft::SpacecraftDynamics as SpacecraftDynamicsRs;
use crate::cosmic::Bodies as BodiesRs;

// Traits
use crate::dynamics::ForceModel as ForceModelRs;
use crate::dynamics::Dynamics as DynamicsRs;

use crate::python::cosmic::Bodies;
use crate::python::cosmic::Frame;
use crate::python::cosmic::Cosm;

/// nyx_space.dynamics
#[pymodule]
fn dynamics(py: Python, m: &PyModule) -> PyResult<()>
{
    let drag_mod = PyModule::new(py, "drag")?;
    drag(drag_mod)?;
    m.add_submodule(drag_mod)?;
    add_module(py, drag_mod);

    let orbital_mod = PyModule::new(py, "orbital")?;
    orbital(orbital_mod)?;
    m.add_submodule(orbital_mod)?;
    add_module(py, orbital_mod);

    let solarpressure_mod = PyModule::new(py, "solarpressure")?;
    solarpressure(solarpressure_mod)?;
    m.add_submodule(solarpressure_mod)?;
    add_module(py, solarpressure_mod);

    let spacecraft_mod = PyModule::new(py, "spacecraft")?;
    spacecraft(spacecraft_mod)?;
    m.add_submodule(spacecraft_mod)?;
    add_module(py, spacecraft_mod);

    // m.add_submodule(wrap_pymodule!(orbital))?;
    // m.add_submodule(wrap_pymodule!(solarpressure))?;
    // m.add_submodule(wrap_pymodule!(spacecraft))?;
    Ok(())
}
/// Ref: https://github.com/PyO3/pyo3/issues/471#issuecomment-489583516
fn add_module(py: Python, module: &PyModule) {
    py.import("sys")
        .expect("failed to import python sys module")
        .dict()
        .get_item("modules")
        .expect("failed to get python modules dictionary")
        .downcast::<PyDict>()
        .expect("failed to turn sys.modules into a PyDict")
        .set_item(module.name().expect("module missing name"), module)
        .expect("failed to inject module");
}
// #[pymodule]
fn drag(m: &PyModule) -> PyResult<()>
{
    m.add_class::<Drag>()?;
    Ok(())
}

// #[pymodule]
fn orbital(m: &PyModule) -> PyResult<()>
{
    m.add_class::<OrbitalDynamics>()?;
    Ok(())
}
// #[pymodule]
fn solarpressure(m: &PyModule) -> PyResult<()>
{
    m.add_class::<SolarPressure>()?;
    Ok(())
}
// #[pymodule]
fn spacecraft(m: &PyModule) -> PyResult<()>
{
    m.add_class::<SpacecraftDynamics>()?;
    Ok(())
}


#[pyclass]
pub struct Drag
{
    pub inner: Arc<DragRs>
}
#[pymethods]
impl Drag {
    #[classmethod]
    /// Common exponential drag model for the Earth
    pub fn earth_exp(_cls: &PyType, cosm: PyRef<Cosm> ) -> Self {
        Self { 
            inner: DragRs::std_atm1976(cosm.inner.clone())
        }
    }
    #[classmethod]
    /// Drag model which uses the standard atmosphere 1976 model for atmospheric density
    pub fn std_atm1976(_cls: &PyType, cosm: PyRef<Cosm> ) -> Self {
        Self { 
            inner: DragRs::std_atm1976(cosm.inner.clone())
        }
    }
    pub fn force_model(&self) -> ForceModel {
        ForceModel {
            inner: self.inner.clone()
        }
    }
}

#[pyclass]
pub struct SolarPressure
{
    pub inner: Arc<SolarPressureRs>
}
#[pymethods]
impl SolarPressure {
    #[classmethod]
    pub fn default(_cls: &PyType, frame: PyRef<Frame>, cosm: PyRef<Cosm> ) -> Self {
        Self { 
            inner: SolarPressureRs::default(frame.inner.clone(), cosm.inner.clone())
        }
    }
    pub fn force_model(&self) -> ForceModel {
        ForceModel {
            inner: self.inner.clone()
        }
    }
}

#[pyclass]
pub struct ForceModel {
    pub inner: Arc<dyn ForceModelRs>
}
#[pymethods]
impl ForceModel {}

#[pyclass]
pub struct OrbitalDynamics
{
    pub inner: Arc<OrbitalDynamicsRs<'static>>
}
#[pymethods]
impl OrbitalDynamics {
    #[classmethod]
    fn point_masses(_cls: &PyType, bodies: Vec<Bodies>, cosm: PyRef<Cosm>) -> Self {
        let bodies_slice: Vec<BodiesRs> = bodies.iter().map(|b| b.inner).collect();
        Self { inner: OrbitalDynamicsRs::point_masses(&bodies_slice, cosm.inner.clone())}
    }
}

#[pyclass]
pub struct SpacecraftDynamics {
    pub inner: Arc<SpacecraftDynamicsRs<'static>>
}
#[pymethods]
impl SpacecraftDynamics {
    #[classmethod]
    pub fn from_models(
        _cls: &PyType, 
        orbital_dyn: PyRef<OrbitalDynamics>, 
        force_models: Vec<PyRef<ForceModel>>
    ) -> Self
    {
        let force_models_rs: Vec<_> = force_models
                                        .iter()
                                        .map(|f| f.inner.clone())
                                        .collect();
        Self{
            inner: SpacecraftDynamicsRs::from_models(orbital_dyn.inner.clone(), force_models_rs)
        }
    }
}
