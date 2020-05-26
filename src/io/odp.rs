pub use crate::celestia::*;
use crate::io::scenario::ScenarioSerde;
use crate::io::ParsingError;
use crate::md::ui::{MDProcess, StmStateFlag};
use crate::od::ranging::GroundStation;
use crate::od::{Measurement, MeasurementDevice};
use crate::time::SECONDS_PER_DAY;
use std::sync::mpsc::channel;
use std::time::Instant;

pub struct OdpScenario<'a> {
    cosm: &'a Cosm,
}

impl<'a> OdpScenario<'a> {
    pub fn try_from_scenario(
        scenario: &ScenarioSerde,
        seq_name: String,
        cosm: &'a Cosm,
    ) -> Result<Self, ParsingError> {
        if scenario.odp.is_none() {
            return Err(ParsingError::UseMdInstead);
        }

        let odp = &scenario.odp.as_ref().unwrap();
        // Ensure valid ODP scenario
        if scenario.measurements.is_none() {
            return Err(ParsingError::OD(
                "missing `measurements` section".to_string(),
            ));
        } else if scenario.estimate.is_none() {
            return Err(ParsingError::OD("missing `estimate` section".to_string()));
        } else if scenario.stations.is_none() {
            return Err(ParsingError::OD("missing `stations` section".to_string()));
        }
        let all_measurements = scenario.measurements.as_ref().unwrap();
        let all_estimates = scenario.estimate.as_ref().unwrap();
        let all_stations = scenario.stations.as_ref().unwrap();

        if let Some(odp_seq) = odp.get(&seq_name) {
            // Get the measurement generation
            match all_measurements.get(&odp_seq.measurements) {
                None => unimplemented!(),
                Some(ref msr) => {
                    // Get the IAU Earth frame
                    let iau_earth = cosm.frame("IAU Earth");
                    // Build the stations
                    let mut stations = Vec::with_capacity(5);
                    for device in &msr.msr_device {
                        match all_stations.get(device) {
                            None => {
                                return Err(ParsingError::OD(format!(
                                    "station `{}` in sequence `{}` not found",
                                    device, seq_name
                                )))
                            }
                            Some(s) => {
                                let gs = if let Some(base) = &s.from {
                                    match base.to_lowercase().as_str() {
                                        "dss13" => GroundStation::dss13_goldstone(
                                            s.elevation,
                                            s.range_noise,
                                            s.range_rate_noise,
                                            cosm,
                                        ),
                                        "dss34" => GroundStation::dss34_canberra(
                                            s.elevation,
                                            s.range_noise,
                                            s.range_rate_noise,
                                            cosm,
                                        ),
                                        "dss65" => GroundStation::dss65_madrid(
                                            s.elevation,
                                            s.range_noise,
                                            s.range_rate_noise,
                                            cosm,
                                        ),
                                        _ => {
                                            return Err(ParsingError::OD(format!(
                                                "unknown base station `{}`",
                                                base
                                            )))
                                        }
                                    }
                                } else {
                                    let station_name = device.clone();
                                    GroundStation::from_noise_values(
                                        station_name,
                                        s.elevation,
                                        s.latitude.unwrap(),
                                        s.longitude.unwrap(),
                                        s.height.unwrap(),
                                        s.range_noise,
                                        s.range_rate_noise,
                                        iau_earth,
                                        cosm,
                                    )
                                };
                                stations.push(gs);
                            }
                        }
                    }
                    // Generate the measurements.
                    let mut md = MDProcess::try_from_scenario(
                        scenario,
                        msr.propagator.as_ref().unwrap().to_string(),
                        StmStateFlag::Without(()),
                        cosm,
                    )?;

                    let prop_time = md.prop_time_s.unwrap();
                    let mut prop = md.propagator();

                    // Set up the channels
                    let (tx, rx) = channel();
                    prop.tx_chan = Some(&tx);

                    // Run
                    info!(
                        "Generating measurements over {} seconds (~ {:.3} days)",
                        prop_time,
                        prop_time / SECONDS_PER_DAY
                    );
                    let start = Instant::now();
                    prop.until_time_elapsed(prop_time);
                    info!(
                        "Done in {:.3} seconds",
                        (Instant::now() - start).as_secs_f64()
                    );

                    let mut sim_measurements = Vec::with_capacity(10000);
                    while let Ok(rx_state) = rx.try_recv() {
                        for station in stations.iter() {
                            let meas = station.measure(&rx_state).unwrap();
                            if meas.visible() {
                                sim_measurements.push(meas);
                                break; // We know that only one station is in visibility at each time.
                            }
                        }
                    }
                }
            }
        }

        Ok(Self { cosm })
    }
}
