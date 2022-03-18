/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2022 Christopher Rabotin <christopher.rabotin@gmail.com>

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

extern crate clap;
extern crate config;
extern crate dialoguer;
extern crate glob;
extern crate lazy_static;
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;
extern crate rust_embed;

use clap::{App, Arg};
use config::{Config, File};
use dialoguer::{theme::ColorfulTheme, Select};
use glob::glob;
use lazy_static::lazy_static;
use log::{debug, error, info};
use nyx::cosmic::{Cosm, Xb};
use nyx::io::{odp::OdpScenario, scenario::*, ParsingError};
use nyx::md::ui::MDProcess;
use nyx::md::MdHdlr;
use nyx::md::OrbitStateOutput;
use nyx::Spacecraft;
use rust_embed::RustEmbed;
use std::borrow::Cow;
use std::env::{set_var, var};
use std::sync::Arc;

const LOG_VAR: &str = "NYX_LOG";

lazy_static! {
    static ref COSM: Arc<Cosm> = {
        let de438_buf: Cow<'static, [u8]> = EmbeddedAsset::get("de438s-00-50.xb")
            .expect("Could not find de438s-00-55.xb as asset")
            .data;
        let xb = Xb::from_buffer(&de438_buf).unwrap();
        let mut cosm: Cosm = Cosm::try_from_xb(xb).unwrap();
        cosm.use_gmat_gm();
        debug!("Loaded\n{}\n", cosm);
        Arc::new(cosm)
    };
}

#[derive(RustEmbed)]
#[folder = "data/embed/"]
struct EmbeddedAsset;

fn main() -> Result<(), ParsingError> {
    let app = App::new("nyx")
        .version("1.0.0-alpha.2")
        .author("Chris Rabotin <chris.rabotin@pm.me>")
        .about("Mission design, orbit determination and Monte-Carlo tool.")
        .arg(
            Arg::with_name("SCENARIO")
                .help("Sets the scenario file or path to use")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::with_name("sequence")
                .short("s")
                .long("seq")
                .takes_value(true)
                .value_name("sequence name")
                .help("Specify which sequence to run, program will fail if sequence is not found"),
        )
        .arg(
            Arg::with_name("all")
                .short("a")
                .long("all")
                .takes_value(false)
                .help("Execute all sequences in order"),
        );

    let matches = app.get_matches();

    let mut s = Config::new();

    // Start off by merging in the "default" configuration file
    let scenario_path = matches.value_of("SCENARIO").unwrap();
    if scenario_path.contains('*') {
        s.merge(
            glob(scenario_path)
                .unwrap()
                .map(|path| File::from(path.unwrap()))
                .collect::<Vec<_>>(),
        )
        .expect("Could not load scenario from folder");
    } else {
        s.merge(File::with_name(scenario_path))
            .expect("Could not load scenario from file");
    }

    let exec_all = matches.is_present("all");
    // Try to deserialize the scenario
    let scenario: ScenarioSerde;
    match s.try_into() {
        Ok(s) => scenario = s,
        Err(e) => return Err(ParsingError::LoadingError(e.to_string())),
    };

    if var(LOG_VAR).is_err() {
        set_var(LOG_VAR, "INFO");
    }

    if pretty_env_logger::try_init_custom_env(LOG_VAR).is_err() {
        println!("could not init logger");
    }

    info!("Loaded scenario `{}`", scenario_path);

    // Select the sequence to run (or run all)
    let req_seq_name = if exec_all {
        None
    } else if let Some(seq_name) = matches.value_of("sequence") {
        Some(seq_name.to_string())
    } else {
        // Build the list of sequences
        let sequences = &scenario.sequence;

        if sequences.len() == 1 {
            // Always execute the only sequence available
            Some(sequences[0].clone())
        } else {
            let selection = Select::with_theme(&ColorfulTheme::default())
                .with_prompt("\n\nSelect the sequence to execute")
                .default(0)
                .items(&sequences[..])
                .interact()
                .unwrap();
            Some(sequences[selection].clone())
        }
    };

    // If there is an ODP setup, let's try to decode that
    for seq_name in &scenario.sequence {
        let should_exec = if let Some(req_seq_name) = &req_seq_name {
            seq_name == req_seq_name
        } else {
            true
        };
        if should_exec {
            match OdpScenario::try_from_scenario(&scenario, seq_name.to_string(), (*COSM).clone()) {
                Ok(odp) => odp.execute()?,
                Err(e) => match e {
                    ParsingError::UseMdInstead | ParsingError::SequenceNotFound(_) => {
                        // Build the MDP
                        match MDProcess::try_from_scenario(
                            &scenario,
                            seq_name.to_string(),
                            false,
                            (*COSM).clone(),
                        ) {
                            Ok((mut md, maybe_fmtr)) => {
                                let mut hdlrs: Vec<Box<dyn MdHdlr<Spacecraft>>> = Vec::new();
                                if let Some(fmtr) = maybe_fmtr {
                                    let out =
                                        Box::new(OrbitStateOutput::new(fmtr.clone()).unwrap());
                                    hdlrs.push(out);
                                }

                                info!("*******************************{:*<1$}", "", seq_name.len());
                                info!("===> Executing sequence `{}` <===", seq_name);
                                info!("*******************************{:*<1$}", "", seq_name.len());
                                md.execute_with(hdlrs)?;
                            }
                            Err(e) => {
                                error!("{:?}", e);
                                return Err(e);
                            }
                        };
                    }
                    _ => {
                        error!("{:?}", e);
                        return Err(e);
                    }
                },
            }
        }
    }

    Ok(())
}
