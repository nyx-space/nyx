extern crate clap;
extern crate config;
extern crate dialoguer;
extern crate glob;
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;
extern crate rust_embed;

use clap::{App, Arg};
use config::{Config, File};
use dialoguer::{theme::ColorfulTheme, Select};
use glob::glob;
use log::{error, info};
use nyx::celestia::{load_ephemeris_from_buf, Cosm};
use nyx::io::{odp::OdpScenario, scenario::*, ParsingError};
use nyx::md::ui::{MDProcess, StmStateFlag};
use rust_embed::RustEmbed;
use std::env::{set_var, var};

const LOG_VAR: &str = "NYX_LOG";

#[derive(RustEmbed)]
#[folder = "data/embed/"]
struct EmbeddedAsset;

fn main() -> Result<(), ParsingError> {
    let app = App::new("nyx")
        .version("0.0.20")
        .author(
            "\tChris Rabotin <rabotin@advancedspace.com>\n\tSai Chikine <sai.chikine@advancespace.com>",
        )
        .about("Mission design, orbit determination and Monte-Carlo tool.")
        .arg(
            Arg::with_name("SCENARIO")
                .help("Sets the scenario file or path to use")
                .required(true)
                .index(1),
        ).arg(
            Arg::with_name("sequence")
                .short("s")
                .long("seq")
                .takes_value(true)
                .value_name("sequence name")
                .help("Specify which sequence to run, program will fail if sequence is not found")
        ).arg(
            Arg::with_name("all")
                .short("a")
                .long("all")
                .takes_value(false)
                .help("Execute all sequences in order")
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

    // Load cosm from the embedded
    let de438_buf: Vec<u8> = EmbeddedAsset::get("de438s-00-50.exb")
        .expect("Could not find de438s-00-550.exb as asset")
        .to_vec();
    let cosm = Cosm::try_from_xb(load_ephemeris_from_buf(de438_buf).unwrap()).unwrap();

    // If there is an ODP setup, let's try to decode that

    for seq_name in &scenario.sequence {
        let should_exec = if let Some(req_seq_name) = &req_seq_name {
            seq_name == req_seq_name
        } else {
            true
        };
        if should_exec {
            match OdpScenario::try_from_scenario(&scenario, seq_name.to_string(), &cosm) {
                Ok(odp) => {
                    odp.execute();
                }
                Err(e) => match e {
                    ParsingError::UseMdInstead => {
                        // Build the MDP
                        match MDProcess::try_from_scenario(
                            &scenario,
                            seq_name.to_string(),
                            StmStateFlag::Without(()),
                            &cosm,
                        ) {
                            Ok(md) => {
                                info!("Executing sequence `{}`", seq_name);
                                md.execute();
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
