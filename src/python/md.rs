use pyo3::prelude::*;
use pyo3::types::PyType;

use crate::md::ui::MDProcess as MDProcessRs;
use crate::md::MdHdlr;
use crate::md::OrbitStateOutput;
use std::sync::Arc;

use crate::io::{odp::OdpScenario as OdpScenarioRs, ParsingError};
use crate::python::cosmic::Cosm;
use crate::python::io::ScenarioSerde;
use crate::Spacecraft;

/// nyx_space.md
#[pymodule]
fn md(py: Python, m: &PyModule) -> PyResult<()> {
    m.add_class::<MDProcess>()?;
    Ok(())
}

#[pyclass]
pub struct MDProcess {
    // pub inner: MDProcessRs<'static>
}

#[pymethods]
impl MDProcess {
    #[classmethod]
    pub fn execute_from_scenario(
        _cls: &PyType,
        scenario: PyRef<ScenarioSerde>,
        prop_name: String,
        stm_flag: bool,
        cosm: PyRef<Cosm>,
    ) {
        let (mut md_process, maybe_fmtr) = MDProcessRs::try_from_scenario(
            &scenario.inner,
            prop_name,
            stm_flag,
            Arc::clone(&cosm.inner),
        )
        .unwrap();

        let mut hdlrs: Vec<Box<dyn MdHdlr<Spacecraft>>> = Vec::new();
        if let Some(fmtr) = maybe_fmtr {
            let out = Box::new(OrbitStateOutput::new(fmtr.clone()).unwrap());
            hdlrs.push(out);
        }

        // println!("*******************************{:*<1$}", "", seq_name.len());
        // println!("===> Executing sequence `{}` <===", seq_name);
        // println!("*******************************{:*<1$}", "", seq_name.len());
        md_process.execute_with(hdlrs).unwrap();
    }
    #[classmethod]
    pub fn execute_all_in_scenario(
        _cls: &PyType,
        scenario: PyRef<ScenarioSerde>,
        cosm: PyRef<Cosm>,
    ) {
        let cosm = Arc::clone(&cosm.inner);
        let scenario = &scenario.inner;
        for seq_name in &scenario.sequence {
            match OdpScenarioRs::try_from_scenario(&scenario, seq_name.to_string(), cosm.clone()) {
                Ok(odp) => odp.execute().unwrap(),
                Err(e) => match e {
                    ParsingError::UseMdInstead => {
                        // Build the MDP
                        match MDProcessRs::try_from_scenario(
                            scenario,
                            seq_name.to_string(),
                            false,
                            cosm.clone(),
                        ) {
                            Ok((mut md, maybe_fmtr)) => {
                                let mut hdlrs: Vec<Box<dyn MdHdlr<Spacecraft>>> = Vec::new();
                                if let Some(fmtr) = maybe_fmtr {
                                    let out =
                                        Box::new(OrbitStateOutput::new(fmtr.clone()).unwrap());
                                    hdlrs.push(out);
                                }

                                println!(
                                    "*******************************{:*<1$}",
                                    "",
                                    seq_name.len()
                                );
                                println!("===> Executing sequence `{}` <===", seq_name);
                                println!(
                                    "*******************************{:*<1$}",
                                    "",
                                    seq_name.len()
                                );
                                md.execute_with(hdlrs).unwrap();
                            }
                            Err(e) => {
                                panic!("{:?}", e);
                            }
                        };
                    }
                    _ => {
                        panic!("{:?}", e);
                    }
                },
            }
        }
    }
}
