/*
 * ANISE Toolkit
 * Copyright (C) 2021-onward Christopher Rabotin <christopher.rabotin@gmail.com> et al. (cf. AUTHORS.md)
 * This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at https://mozilla.org/MPL/2.0/.
 *
 * Documentation: https://nyxspace.com/
 */

use std::path::PathBuf;

use anise::naif::kpl::parser::{convert_fk as convert_fk_rs, convert_tpc as convert_tpc_rs};
use anise::structure::dataset::DataSetError;
use pyo3::prelude::*;

#[pymodule]
pub(crate) fn utils(_py: Python, sm: &Bound<PyModule>) -> PyResult<()> {
    sm.add_function(wrap_pyfunction!(convert_fk, sm)?)?;
    sm.add_function(wrap_pyfunction!(convert_tpc, sm)?)?;

    Ok(())
}

/// Converts a KPL/FK file, that defines frame constants like fixed rotations, and frame name to ID mappings into the EulerParameterDataSet equivalent ANISE file.
/// KPL/FK files must be converted into "PCA" (Planetary Constant ANISE) files before being loaded into ANISE.
///
/// :type fk_file_path: str
/// :type anise_output_path: str
/// :type show_comments: bool, optional
/// :type overwrite: bool, optional
/// :rtype: None
#[pyfunction]
#[pyo3(signature = (fk_file_path, anise_output_path, show_comments=None, overwrite=None))]
fn convert_fk(
    fk_file_path: String,
    anise_output_path: String,
    show_comments: Option<bool>,
    overwrite: Option<bool>,
) -> Result<(), DataSetError> {
    let dataset = convert_fk_rs(fk_file_path, show_comments.unwrap_or(false))?;

    dataset.save_as(
        &PathBuf::from(anise_output_path),
        overwrite.unwrap_or(false),
    )?;

    Ok(())
}

/// Converts two KPL/TPC files, one defining the planetary constants as text, and the other defining the gravity parameters, into the PlanetaryDataSet equivalent ANISE file.
/// KPL/TPC files must be converted into "PCA" (Planetary Constant ANISE) files before being loaded into ANISE.
///
/// :type pck_file_path: str
/// :type gm_file_path: str
/// :type anise_output_path: str
/// :type overwrite: bool, optional
/// :rtype: None
#[pyfunction]
#[pyo3(signature = (pck_file_path, gm_file_path, anise_output_path, overwrite=None))]
fn convert_tpc(
    pck_file_path: String,
    gm_file_path: String,
    anise_output_path: String,
    overwrite: Option<bool>,
) -> Result<(), DataSetError> {
    let dataset = convert_tpc_rs(pck_file_path, gm_file_path)?;

    dataset.save_as(
        &PathBuf::from(anise_output_path),
        overwrite.unwrap_or(false),
    )?;

    Ok(())
}
