extern crate clap;
extern crate config;
extern crate log;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;
use clap::{App, Arg};
use config::{Config, File};
use log::info;
use nyx::celestia::Cosm;
use nyx::io::scenario::*;
use nyx::md::ui::MDProcess;
use std::env::{set_var, var};

const LOG_VAR: &str = "NYX_LOG";

fn main() {
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
        );

    let matches = app.get_matches();

    let mut s = Config::new();

    // Start off by merging in the "default" configuration file
    let scenario_path = matches.value_of("SCENARIO").unwrap();
    s.merge(File::with_name(scenario_path))
        .expect("Could not load scenario from file or directory");

    // Try to deserialize the scenario
    let scenario: ScenarioSerde;
    match s.try_into() {
        Ok(s) => scenario = s,
        Err(e) => panic!("{:?}", e),
    };

    if var(LOG_VAR).is_err() {
        set_var(LOG_VAR, "INFO");
    }

    if pretty_env_logger::try_init_custom_env(LOG_VAR).is_err() {
        println!("could not init logger");
    }

    // Load cosm
    let cosm = Cosm::de438();

    match MDProcess::try_from_scenario(scenario, &cosm) {
        Ok(sequence) => {
            info!("Loaded scenario `{}`", scenario_path);
            for mut item in sequence {
                item.execute();
            }
        }
        Err(e) => println!("ERROR\n{:?}", e), // Here is where we try to load this as an ODP
    };
}
