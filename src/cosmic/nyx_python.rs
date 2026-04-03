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

use crate::cosmic::GuidanceMode;
use crate::dynamics::guidance::Thruster;
use crate::Spacecraft;
use anise::prelude::Orbit;
use anise::structure::spacecraft::{DragData, Mass, SRPData};
use pyo3::prelude::*;

#[pymethods]
impl Spacecraft {
    #[pyo3(signature=(orbit, mass=None, srp=None, drag=None, thruster=None, mode=None))]
    #[new]
    fn py_new(
        orbit: Orbit,
        mass: Option<Mass>,
        srp: Option<SRPData>,
        drag: Option<DragData>,
        thruster: Option<Thruster>,
        mode: Option<GuidanceMode>,
    ) -> Self {
        let mut sc = Spacecraft {
            orbit,
            thruster,
            ..Default::default()
        };
        if let Some(mass) = mass {
            sc.mass = mass;
        }
        if let Some(srp) = srp {
            sc.srp = srp;
        }
        if let Some(drag) = drag {
            sc.drag = drag;
        }
        if let Some(mode) = mode {
            sc.mode = mode;
        }

        sc
    }

    #[getter]
    fn orbit(&self) -> Orbit {
        self.orbit
    }

    #[getter]
    fn mass(&self) -> Mass {
        self.mass
    }
    #[getter]
    fn srp(&self) -> SRPData {
        self.srp
    }
    #[getter]
    fn drag(&self) -> DragData {
        self.drag
    }
    fn __str__(&self) -> String {
        format!("{self}")
    }

    fn __repr__(&self) -> String {
        format!("{self} @ {self:p}")
    }

    /// Decodes an ASN.1 DER encoded byte array into a Mass object.
    ///
    /// :type data: bytes
    /// :rtype: Mass
    #[classmethod]
    pub fn from_asn1(_cls: &Bound<'_, PyType>, data: &[u8]) -> PyResult<Self> {
        match Self::from_der(data) {
            Ok(obj) => Ok(obj),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 decoding error: {e}"))),
        }
    }

    /// Encodes this Mass object into an ASN.1 DER encoded byte array.
    ///
    /// :rtype: bytes
    pub fn to_asn1<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyBytes>> {
        let mut buf = Vec::new();
        match self.encode_to_vec(&mut buf) {
            Ok(_) => Ok(PyBytes::new(py, &buf)),
            Err(e) => Err(PyValueError::new_err(format!("ASN.1 encoding error: {e}"))),
        }
    }
}
