extern crate clap;
extern crate config;
extern crate dialoguer;
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;
use clap::{App, Arg};
use config::{Config, File};
use dialoguer::{theme::ColorfulTheme, Select};
use log::{error, info};
use nyx::celestia::Cosm;
use nyx::io::{scenario::*, ParsingError};
use nyx::md::ui::MDProcess;
use std::env::{set_var, var};

const LOG_VAR: &str = "NYX_LOG";

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
    s.merge(File::with_name(scenario_path))
        .expect("Could not load scenario from file or directory");

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
    let seq_name = if exec_all {
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

    // Load cosm
    let cosm = Cosm::de438();

    for scen_seq_name in &scenario.sequence {
        if let Some(seq_name) = &seq_name {
            if scen_seq_name == seq_name {
                // Build the MDP
                match MDProcess::try_from_scenario(&scenario, seq_name.to_string(), &cosm) {
                    Ok(mut md) => {
                        info!("Executing sequence `{}`", seq_name);
                        md.execute();
                        break;
                    }
                    Err(e) => {
                        error!("{:?}", e);
                        return Err(e);
                    }
                };
            }
        } else {
            // Build the MDP
            match MDProcess::try_from_scenario(&scenario, scen_seq_name.to_string(), &cosm) {
                Ok(mut md) => {
                    info!("Executing sequence `{}`", scen_seq_name);
                    md.execute();
                    break;
                }
                Err(e) => {
                    error!("{:?}", e);
                    return Err(e);
                }
            };
        }
    }

    Ok(())
}
