extern crate clap;
extern crate config;
extern crate nyx_space as nyx;
extern crate pretty_env_logger;
use clap::{App, Arg};
use config::{Config, File};
use nyx::celestia::Cosm;
use nyx::io::scenario::*;
use nyx::md::ui::MDProcess;

fn main() {
    let app = App::new("nyx")
        .version("0.0.20")
        .author(
            "\tChris Rabotin <rabotin@advancedspace.com>\n\tSai Chikine <sai.chikine@advancespace.com>",
        )
        .about("Mission design, orbit determination and Monte-Carlo tool")
        .arg(
            Arg::with_name("SCENARIO")
                .help("Sets the scenario file or path to use")
                .required(true)
                .index(1),
        )
        .arg(
            Arg::with_name("v")
                .short("v")
                .multiple(true)
                .help("Sets the level of verbosity"),
        );

    let matches = app.get_matches();

    let mut s = Config::new();

    // Start off by merging in the "default" configuration file
    s.merge(File::with_name(matches.value_of("SCENARIO").unwrap()))
        .expect("Could not load scenario from file or directory");

    // Try to deserialize the scenario
    let scenario: ScenarioSerde;
    match s.try_into() {
        Ok(s) => scenario = s,
        Err(e) => panic!("{:?}", e),
    };

    if pretty_env_logger::try_init().is_err() {
        println!("could not init env_logger");
    }

    // Load cosm
    let cosm = Cosm::from_xb("./de438s");

    match MDProcess::try_from_scenario(scenario, &cosm) {
        Ok(sequence) => {
            for mut item in sequence {
                item.execute();
            }
        }
        Err(e) => println!("ERROR\n{:?}", e),
    };
}
