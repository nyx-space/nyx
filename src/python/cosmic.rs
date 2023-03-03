use pyo3::prelude::*;
use pyo3::types::PyType;

use super::Epoch;
use crate::cosmic::Bodies as BodiesRs;
use crate::cosmic::Cosm as CosmRs;
use crate::cosmic::Frame as FrameRs;
pub use crate::cosmic::Orbit;
use crate::io::orbit::OrbitSerde;
use crate::io::{ConfigError, ConfigRepr, Configurable};
use std::sync::Arc;

#[pyclass]
#[derive(Clone, Copy)]
pub struct Bodies {
    pub inner: BodiesRs,
}
#[allow(non_snake_case)]
#[pymethods]
impl Bodies {
    #[classattr]
    pub fn SSB() -> Self {
        Self {
            inner: BodiesRs::SSB,
        }
    }
    #[classattr]
    pub fn Sun() -> Self {
        Self {
            inner: BodiesRs::Sun,
        }
    }
    #[classattr]
    pub fn MercuryBarycenter() -> Self {
        Self {
            inner: BodiesRs::MercuryBarycenter,
        }
    }
    #[classattr]
    pub fn Mercury() -> Self {
        Self {
            inner: BodiesRs::Mercury,
        }
    }
    #[classattr]
    pub fn VenusBarycenter() -> Self {
        Self {
            inner: BodiesRs::VenusBarycenter,
        }
    }
    #[classattr]
    pub fn Venus() -> Self {
        Self {
            inner: BodiesRs::Venus,
        }
    }
    #[classattr]
    pub fn EarthBarycenter() -> Self {
        Self {
            inner: BodiesRs::EarthBarycenter,
        }
    }
    #[classattr]
    pub fn Earth() -> Self {
        Self {
            inner: BodiesRs::Earth,
        }
    }
    #[classattr]
    pub fn Luna() -> Self {
        Self {
            inner: BodiesRs::Luna,
        }
    }
    #[classattr]
    pub fn MarsBarycenter() -> Self {
        Self {
            inner: BodiesRs::MarsBarycenter,
        }
    }
    #[classattr]
    pub fn JupiterBarycenter() -> Self {
        Self {
            inner: BodiesRs::JupiterBarycenter,
        }
    }
    #[classattr]
    pub fn SaturnBarycenter() -> Self {
        Self {
            inner: BodiesRs::SaturnBarycenter,
        }
    }
    #[classattr]
    pub fn UranusBarycenter() -> Self {
        Self {
            inner: BodiesRs::UranusBarycenter,
        }
    }
    #[classattr]
    pub fn NeptuneBarycenter() -> Self {
        Self {
            inner: BodiesRs::NeptuneBarycenter,
        }
    }
    #[classattr]
    pub fn PlutoBarycenter() -> Self {
        Self {
            inner: BodiesRs::PlutoBarycenter,
        }
    }
}

#[pymethods]
impl Orbit {
    /// Creates a new Orbit in the provided frame at the provided Epoch.
    ///
    /// **Units:** km, km, km, km/s, km/s, km/s
    #[new]
    pub fn from_cartesian(
        x_km: f64,
        y_km: f64,
        z_km: f64,
        vx_km_s: f64,
        vy_km_s: f64,
        vz_km_s: f64,
        epoch: Epoch,
        frame: PyRef<Frame>,
    ) -> Self {
        Self::cartesian(
            x_km,
            y_km,
            z_km,
            vx_km_s,
            vy_km_s,
            vz_km_s,
            epoch,
            frame.inner,
        )
    }

    #[staticmethod]
    fn load_yaml(path: &str) -> Result<Self, ConfigError> {
        let serde = OrbitSerde::load_yaml(path)?;

        let cosm = CosmRs::de438();

        Self::from_config(&serde, cosm)
    }

    #[staticmethod]
    fn load_many_yaml(path: &str) -> Result<Vec<Self>, ConfigError> {
        let orbits = OrbitSerde::load_many_yaml(path)?;

        // Create a new Cosm until ANISE switch
        let cosm = CosmRs::de438();

        let mut selves = Vec::with_capacity(orbits.len());

        for serde in orbits {
            selves.push(Self::from_config(&serde, cosm.clone())?);
        }

        Ok(selves)
    }

    fn __repr__(&self) -> String {
        format!("{self:?}")
    }

    fn __str__(&self) -> String {
        format!("{self}")
    }
}

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
    pub fn angualar_velocity(&self) -> f64 {
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
