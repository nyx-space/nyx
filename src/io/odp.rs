extern crate csv;

pub use crate::celestia::*;
use crate::dimensions::{Matrix2, Matrix3, Matrix6, Vector2, Vector3, Vector6, U2, U3, U6};
use crate::dynamics::spacecraft::SpacecraftState;
use crate::io::formatter::NavSolutionFormatter;
use crate::io::scenario::ScenarioSerde;
use crate::io::{parse_duration, ParsingError};
use crate::md::ui::{MDProcess, StmStateFlag};
use crate::od::ranging::GroundStation;
use crate::od::ui::*;
use crate::od::{Measurement, MeasurementDevice};
use crate::propagators::PropOpts;
use crate::time::SECONDS_PER_DAY;
use std::sync::mpsc::channel;
use std::time::Instant;

pub struct OdpScenario<'a> {
    truth: MDProcess<'a>,
    nav: MDProcess<'a>,
    ekf_msr_trigger: usize,
    kf: KF<U6, U3, U2, SpacecraftState>,
    stations: Vec<GroundStation<'a>>,
    formatter: Option<NavSolutionFormatter<'a>>,
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
                    // Build the truth MD
                    let md = MDProcess::try_from_scenario(
                        scenario,
                        msr.propagator.as_ref().unwrap().to_string(),
                        StmStateFlag::Without(()),
                        cosm,
                    )?;

                    // Build the initial estimate
                    if all_estimates.get(&odp_seq.initial_estimate).is_none() {
                        return Err(ParsingError::OD(format!(
                            "estimate `{}` not found",
                            &odp_seq.initial_estimate
                        )));
                    }

                    let est_serde = &all_estimates[&odp_seq.initial_estimate];
                    if scenario.state.get(&est_serde.state).is_none() {
                        return Err(ParsingError::OD(format!(
                            "state `{}` not found",
                            &est_serde.state
                        )));
                    }

                    let est_init_state_serde = &scenario.state[&est_serde.state];
                    let state_frame = cosm.frame(est_init_state_serde.frame.as_str());
                    let est_init_state = est_init_state_serde.as_state(state_frame);

                    // Build the covariance
                    let mut cov;
                    if let Some(covar_mat) = &est_serde.covar_mat {
                        cov = Matrix6::from_element(0.0);
                        let max_i = covar_mat.len();
                        let mut max_j = 0;
                        for (i, row) in covar_mat.iter().enumerate() {
                            for (j, item) in row.iter().enumerate() {
                                cov[(i, j)] = *item;
                                if j > max_j {
                                    max_j = j;
                                }
                            }
                        }
                        if max_i != 6 || max_j != 6 {
                            return Err(ParsingError::OD(format!(
                                "initial covariance `{}` does not have 36 items",
                                &odp_seq.initial_estimate
                            )));
                        }
                    } else if let Some(covar_diag) = &est_serde.covar_diag {
                        let mut as_v = Vector6::zeros();
                        for (i, item) in covar_diag.iter().enumerate() {
                            as_v[i] = *item;
                        }
                        cov = Matrix6::from_diagonal(&as_v);
                    } else {
                        return Err(ParsingError::OD(format!(
                            "initial covariance not specified in `{}`",
                            &odp_seq.initial_estimate
                        )));
                    }
                    let mut init_sc_state = md.state();
                    init_sc_state.orbit = est_init_state;
                    let initial_estimate = KfEstimate::from_covar(init_sc_state, cov);
                    let measurement_noise = Matrix2::from_diagonal(&Vector2::new(
                        odp_seq.msr_noise[0],
                        odp_seq.msr_noise[1],
                    ));

                    // Build the filter
                    let kf = match &odp_seq.snc {
                        None => KF::no_snc(initial_estimate, measurement_noise),
                        Some(snc) => {
                            if snc.len() != 3 {
                                return Err(ParsingError::OD(
                                    "SNC must have three components".to_string(),
                                ));
                            }
                            let process_noise =
                                Matrix3::from_diagonal(&Vector3::new(snc[0], snc[1], snc[2]));
                            // Disable SNC if there is more than 120 seconds between two measurements
                            let process_noise_dt = match &odp_seq.snc_disable {
                                None => {
                                    warn!("No SNC disable time specified, assuming 120 seconds");
                                    Some(120.0)
                                }
                                Some(snd_disable_dt) => Some(parse_duration(snd_disable_dt)?),
                            };

                            // And build the filter
                            KF::new(
                                initial_estimate,
                                process_noise,
                                measurement_noise,
                                process_noise_dt,
                            )
                        }
                    };

                    // Build the estimation propagator
                    let estimator = MDProcess::try_from_scenario(
                        scenario,
                        odp_seq.navigation_prop.to_string(),
                        StmStateFlag::With(()),
                        cosm,
                    )?;

                    // Get the formatter
                    let formatter = match &odp_seq.output {
                        Some(output) => match &scenario.output.get(output) {
                            None => {
                                return Err(ParsingError::OD(format!(
                                    "output `{}` not found",
                                    output
                                )))
                            }
                            Some(output) => Some(output.to_nav_sol_formatter(&cosm)),
                        },
                        None => None,
                    };

                    return Ok(Self {
                        truth: md,
                        nav: estimator,
                        kf,
                        ekf_msr_trigger: match &odp_seq.ekf_msr_trigger {
                            Some(val) => *val,
                            None => 100_000,
                        },
                        stations,
                        formatter,
                    });
                }
            }
        }

        Err(ParsingError::OD("sequence not found".to_string()))
    }

    /// Will generate the measurements and run the filter.
    pub fn execute(mut self) {
        // Generate the measurements.
        let prop_time = self.truth.prop_time_s.unwrap();
        let mut prop = self.truth.propagator();
        prop.set_step(60.0, true);

        // Set up the channels
        let (tx, rx) = channel();
        prop.tx_chan = Some(tx);

        // Generate the measurements
        info!(
            "Generating measurements over {} seconds (~ {:.3} days)",
            prop_time,
            prop_time / SECONDS_PER_DAY
        );

        let start = Instant::now();
        prop.until_time_elapsed(prop_time);

        let mut sim_measurements = Vec::with_capacity(10000);
        while let Ok(rx_state) = rx.try_recv() {
            for station in self.stations.iter() {
                let meas = station.measure(&rx_state).unwrap();
                if meas.visible() {
                    sim_measurements.push(meas);
                    break; // We know that only one station is in visibility at each time.
                }
            }
        }

        info!(
            "Generated {} measurements in {:.3} seconds",
            sim_measurements.len(),
            (Instant::now() - start).as_secs_f64()
        );

        // Build the ODP
        let nav = self.nav.stm_propagator();
        let kf = self.kf;
        let mut odp = ODProcess::ekf(
            nav,
            kf,
            self.stations.clone(),
            false,
            100_000,
            NumMsrEkfTrigger::init(self.ekf_msr_trigger),
        );

        odp.process_measurements_covar(&sim_measurements);

        // Save to output file if requested
        // Create the output file
        if let Some(fmtr) = &self.formatter {
            let mut wtr =
                csv::Writer::from_path(fmtr.filename.clone()).expect("could not create file");
            wtr.serialize(&fmtr.headers)
                .expect("could not write headers");
            info!("Saving output to {}", fmtr.filename);
            for est in &odp.estimates {
                wtr.serialize(fmtr.fmt(est))
                    .expect("could not format state");
            }
        };
    }
}
