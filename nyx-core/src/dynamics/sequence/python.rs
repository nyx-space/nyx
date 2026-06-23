use super::{AccelModels, Dynamics, ForceModels, PropagatorConfig, SpacecraftSequence, Thruster};
#[cfg(feature = "python")]
use crate::dynamics::{Drag, SolarPressure};
use crate::dynamics::sequence::config::TidalBody;
use crate::propagators::{IntegratorMethod, IntegratorOptions};
use crate::{dynamics::PointMasses, io::gravity::GravityFieldConfig};
use pyo3::exceptions::PyException;
use {
    crate::Spacecraft,
    anise::almanac::Almanac,
    pyo3::prelude::*,
    std::{collections::HashMap, sync::Arc},
};

#[cfg(feature = "python")]
#[pymethods]
impl SpacecraftSequence {
    #[new]
    fn py_new() -> Self {
        SpacecraftSequence::default()
    }

    #[classmethod]
    #[pyo3(name = "from_dhall")]
    fn py_from_dhall(_cls: &Bound<'_, pyo3::types::PyType>, dhall_str: &str) -> PyResult<Self> {
        serde_dhall::from_str(dhall_str)
            .parse()
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[classmethod]
    #[pyo3(name = "from_yaml")]
    fn py_from_yaml(_cls: &Bound<'_, pyo3::types::PyType>, yaml_str: &str) -> PyResult<Self> {
        serde_yml::from_str(yaml_str)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[pyo3(name = "setup")]
    fn py_setup(&mut self, py: Python<'_>, almanac: Py<Almanac>) -> PyResult<()> {
        let almanac_ref = almanac.borrow(py);
        self.setup(Arc::new(almanac_ref.clone()))
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    #[pyo3(name = "propagate")]
    fn py_propagate(
        &self,
        py: Python<'_>,
        state: Spacecraft,
        until_phase: Option<String>,
        almanac: Py<Almanac>,
    ) -> PyResult<Vec<(Option<String>, Vec<Spacecraft>)>> {
        let almanac_ref = almanac.borrow(py);
        let trajs = self
            .propagate(state, until_phase, Arc::new(almanac_ref.clone()))
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        let result = trajs
            .into_iter()
            .map(|traj| (traj.name, traj.states))
            .collect();
        Ok(result)
    }

    #[getter]
    fn get_thruster_sets(&self) -> HashMap<String, Thruster> {
        self.thruster_sets.clone()
    }

    fn thruster_set_insert(&mut self, name: String, thruster: Thruster) {
        self.thruster_sets.insert(name, thruster);
    }
    fn thruster_set_remove(&mut self, name: String) -> PyResult<()> {
        if self.thruster_sets.remove(&name).is_none() {
            Err(PyException::new_err(format!("{name} not in thruster set")))
        } else {
            Ok(())
        }
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl AccelModels {
    #[pyo3(signature=(point_masses=None, gravity_field=None, solid_tides=None))]
    #[new]
    fn py_new(
        point_masses: Option<PointMasses>,
        gravity_field: Option<GravityFieldConfig>,
        solid_tides: Option<Vec<TidalBody>>,
    ) -> PyResult<Self> {
        let mut st = [None, None];
        if let Some(tides) = solid_tides {
            if tides.len() > 2 {
                return Err(pyo3::exceptions::PyValueError::new_err(
                    "A maximum of 2 tidal bodies is supported.",
                ));
            }
            for (i, body) in tides.into_iter().enumerate() {
                st[i] = Some(body);
            }
        }
        Ok(Self {
            point_masses,
            gravity_field,
            solid_tides: st,
        })

    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl ForceModels {
    #[pyo3(signature=(solar_pressure=None, drag=None))]
    #[new]
    fn py_new(solar_pressure: Option<SolarPressure>, drag: Option<Drag>) -> Self {
        Self {
            solar_pressure,
            drag,
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl Dynamics {
    #[pyo3(signature=(accel_models=AccelModels::default(), force_models=ForceModels::default()))]
    #[new]
    fn py_new(accel_models: AccelModels, force_models: ForceModels) -> Self {
        Self {
            accel_models,
            force_models,
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl PropagatorConfig {
    #[new]
    fn py_new(dynamics: Dynamics, method: IntegratorMethod, options: IntegratorOptions) -> Self {
        Self {
            dynamics,
            method,
            options,
        }
    }

    fn __str__(&self) -> String {
        format!("{self:?}")
    }

    fn __repr__(&self) -> String {
        format!("{self:?} @ {self:p}")
    }
}

#[cfg(feature = "python")]
#[cfg_attr(feature = "python", pymethods)]
impl TidalBody {
    #[new]
    fn py_new(frame: anise::frames::FrameUid, compute_degree_3: bool) -> Self {
        Self { frame, compute_degree_3 }
    }
}
