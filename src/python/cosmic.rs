use pyo3::prelude::*;
use pyo3::py_run;
use pyo3::types::PyType;

pub use crate::cosmic::Bodies;
use crate::cosmic::Cosm as CosmRs;
use crate::cosmic::Frame as FrameRs;
use crate::cosmic::GuidanceMode;
pub use crate::cosmic::Orbit;
pub use crate::cosmic::{DragConfig, Spacecraft, SrpConfig};
use crate::dynamics::guidance::Thruster;
use std::sync::Arc;

pub(crate) fn register_cosmic(py: Python<'_>, parent_module: &PyModule) -> PyResult<()> {
    let sm = PyModule::new(py, "_nyx_space.cosmic")?;
    sm.add_class::<Cosm>()?;
    sm.add_class::<Bodies>()?;
    sm.add_class::<Frame>()?;
    sm.add_class::<Orbit>()?;
    sm.add_class::<Spacecraft>()?;
    sm.add_class::<SrpConfig>()?;
    sm.add_class::<DragConfig>()?;
    sm.add_class::<Thruster>()?;
    sm.add_class::<GuidanceMode>()?;

    py_run!(py, sm, "import sys; sys.modules['nyx_space.cosmic'] = sm");
    parent_module.add_submodule(sm)?;
    Ok(())
}

#[derive(Clone)]
#[pyclass]
pub struct Frame {
    pub inner: FrameRs,
}

#[pymethods]
impl Frame {
    pub fn is_geoid(&self) -> bool {
        self.inner.is_geoid()
    }
    pub fn is_celestial(&self) -> bool {
        self.inner.is_celestial()
    }
    pub fn ephem_path(&self) -> Vec<usize> {
        self.inner.ephem_path()
    }
    pub fn frame_path(&self) -> Vec<usize> {
        self.inner.frame_path()
    }
    pub fn gm(&self) -> f64 {
        self.inner.gm()
    }
    pub fn equatorial_radius(&self) -> f64 {
        self.inner.equatorial_radius()
    }
    pub fn flattening(&self) -> f64 {
        self.inner.flattening()
    }
    pub fn semi_major_radius(&self) -> f64 {
        self.inner.semi_major_radius()
    }
    pub fn angular_velocity(&self) -> f64 {
        self.inner.angular_velocity()
    }
}

#[pyclass]
pub struct Cosm {
    pub inner: Arc<CosmRs>,
}

#[pymethods]
impl Cosm {
    #[classmethod]
    pub fn from_xb(_cls: &PyType, filename: &str) -> PyResult<Self> {
        Ok(Cosm {
            inner: Arc::new(CosmRs::from_xb(filename)?),
        })
    }
    #[classmethod]
    pub fn try_de438(_cls: &PyType) -> PyResult<Self> {
        Ok(Cosm {
            inner: Arc::new(CosmRs::try_de438()?),
        })
    }
    #[classmethod]
    pub fn de438(_cls: &PyType) -> PyResult<Self> {
        Ok(Cosm {
            inner: CosmRs::de438(),
        })
    }
    #[classmethod]
    pub fn de438_raw(_cls: &PyType) -> PyResult<Self> {
        Ok(Cosm {
            inner: Arc::new(CosmRs::de438_raw()),
        })
    }
    #[classmethod]
    pub fn de438_gmat(_cls: &PyType) -> PyResult<Self> {
        Ok(Cosm {
            inner: CosmRs::de438_gmat(),
        })
    }

    pub fn frame(&self, name: &str) -> PyResult<Frame> {
        Ok(Frame {
            inner: self.inner.try_frame(name)?,
        })
    }

    pub fn frames_get_names(&self) -> PyResult<Vec<String>> {
        Ok(self.inner.frames_get_names())
    }
}
