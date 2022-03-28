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

extern crate csv;

pub use crate::cosmic::*;
use crate::dynamics::NyxError;
use crate::io::formatter::NavSolutionFormatter;
use crate::io::quantity::{parse_duration, ParsingError};
use crate::io::scenario::ScenarioSerde;
use crate::linalg::{Matrix2, Matrix6, Vector2, Vector6, U2, U3};
use crate::md::ui::MDProcess;
use crate::od::measurement::GroundStation;
use crate::od::ui::snc::SNC3;
use crate::od::ui::*;
use crate::od::MeasurementDevice;
use crate::propagators::Propagator;
use crate::time::{Duration, Unit};
use crate::Orbit;
use std::str::FromStr;
use std::sync::mpsc::channel;
use std::sync::Arc;
use std::time::Instant;

pub struct OdpScenario<'a> {
    truth: MDProcess<'a>,
    nav: MDProcess<'a>,
    ekf_msr_trigger: usize,
    ekf_disable_time: Duration,
    kf: KF<Orbit, U3, U2>,
    stations: Vec<GroundStation>,
    formatter: Option<NavSolutionFormatter>,
}

impl<'a> OdpScenario<'a> {
    #[allow(clippy::identity_op)]
    pub fn try_from_scenario(
        scenario: &ScenarioSerde,
        seq_name: String,
        cosm: Arc<Cosm>,
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
        let all_measurements = match scenario.measurements.as_ref() {
            Some(msr) => msr,
            None => {
                return Err(ParsingError::LoadingError(
                    "no measurements provided".to_string(),
                ))
            }
        };

        let all_estimates = match scenario.estimate.as_ref() {
            Some(est) => est,
            None => {
                return Err(ParsingError::LoadingError(
                    "no estimates provided".to_string(),
                ))
            }
        };

        let all_stations = match scenario.stations.as_ref() {
            Some(st) => st,
            None => {
                return Err(ParsingError::LoadingError(
                    "no stations provided".to_string(),
                ))
            }
        };

        if let Some(odp_seq) = odp.get(&seq_name.to_lowercase()) {
            // Get the measurement generation
            match all_measurements.get(&odp_seq.measurements.to_lowercase()) {
                None => unimplemented!("{}", &odp_seq.measurements),
                Some(msr) => {
                    // Get the IAU Earth frame
                    let iau_earth = cosm.frame("IAU Earth");
                    // Build the stations
                    let mut stations = Vec::with_capacity(5);
                    for device in &msr.msr_device {
                        match all_stations.get(&device.to_lowercase()) {
                            None => {
                                return Err(ParsingError::OD(format!(
                                    "station `{}` in sequence `{}` not found",
                                    device, seq_name
                                )))
                            }
                            Some(s) => {
                                let gs = if let Some(base) = &s.inherit {
                                    match base.to_lowercase().as_str() {
                                        "dss13" => GroundStation::dss13_goldstone(
                                            s.elevation,
                                            s.range_noise,
                                            s.range_rate_noise,
                                            cosm.clone(),
                                        ),
                                        "dss34" => GroundStation::dss34_canberra(
                                            s.elevation,
                                            s.range_noise,
                                            s.range_rate_noise,
                                            cosm.clone(),
                                        ),
                                        "dss65" => GroundStation::dss65_madrid(
                                            s.elevation,
                                            s.range_noise,
                                            s.range_rate_noise,
                                            cosm.clone(),
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
                                        cosm.clone(),
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
                        false,
                        cosm.clone(),
                    )?
                    .0;

                    // Build the initial estimate
                    if all_estimates
                        .get(&odp_seq.initial_estimate.to_lowercase())
                        .is_none()
                    {
                        return Err(ParsingError::OD(format!(
                            "estimate `{}` not found",
                            &odp_seq.initial_estimate
                        )));
                    }

                    let est_serde = &all_estimates[&odp_seq.initial_estimate];
                    if scenario
                        .state
                        .get(&est_serde.state.to_lowercase())
                        .is_none()
                    {
                        return Err(ParsingError::OD(format!(
                            "state `{}` not found",
                            &est_serde.state
                        )));
                    }

                    let est_init_state_serde = &scenario.state[&est_serde.state.to_lowercase()];
                    let state_frame =
                        cosm.frame(est_init_state_serde.frame.as_ref().unwrap().as_str());
                    let est_init_state = est_init_state_serde.as_state(state_frame)?;

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
                    let mut init_sc_state = md.init_state;
                    init_sc_state.orbit = est_init_state;
                    let initial_estimate = KfEstimate::from_covar(init_sc_state.orbit, cov);
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
                            // Disable SNC if there is more than 120 seconds between two measurements
                            let disable_time = match &odp_seq.snc_disable {
                                None => {
                                    warn!("No SNC disable time specified, assuming 120 seconds");
                                    2 * Unit::Minute
                                }
                                Some(snc_disable_dt) => match Duration::from_str(snc_disable_dt) {
                                    Ok(d) => d,
                                    Err(e) => {
                                        return Err(ParsingError::IllDefined(format!(
                                            "When parsing SNC duration: {}",
                                            e
                                        )))
                                    }
                                },
                            };

                            // Build the process noise
                            let process_noise = match &odp_seq.snc_decay {
                                None => SNC3::from_diagonal(disable_time, snc),
                                Some(decay_str) => {
                                    if decay_str.len() != 3 {
                                        return Err(ParsingError::OD(
                                            "SNC decay must have three components".to_string(),
                                        ));
                                    }
                                    let mut scn_decay_s = [0.0; 3];
                                    for (i, ds) in decay_str.iter().enumerate() {
                                        scn_decay_s[i] = parse_duration(ds)?.v();
                                    }
                                    SNC3::with_decay(disable_time, snc, &scn_decay_s)
                                }
                            };

                            info!("Using SNC: {}", process_noise);

                            // And build the filter
                            KF::new(initial_estimate, process_noise, measurement_noise)
                        }
                    };

                    // Build the estimation propagator
                    let estimator = MDProcess::try_from_scenario(
                        scenario,
                        odp_seq.navigation_prop.to_string(),
                        true,
                        cosm.clone(),
                    )?
                    .0;

                    // Get the formatter
                    let formatter = match &odp_seq.output {
                        Some(output) => match &scenario.output.get(&output.to_lowercase()) {
                            None => {
                                return Err(ParsingError::OD(format!(
                                    "output `{}` not found",
                                    output
                                )))
                            }
                            Some(output) => Some(output.to_nav_sol_formatter(cosm)?),
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
                        ekf_disable_time: match &odp_seq.ekf_disable_time {
                            Some(val) => *val,
                            None => 1 * Unit::Hour, // defaults to one hour
                        },
                        stations,
                        formatter,
                    });
                }
            }
        }

        Err(ParsingError::SequenceNotFound(seq_name))
    }

    /// Will generate the measurements and run the filter.
    pub fn execute(self) -> Result<(), NyxError> {
        // Generate the measurements.
        let prop_time = self.truth.prop_time.unwrap();

        // Create the output file for the truth
        let mut maybe_wtr = match &self.truth.formatter {
            Some(fmtr) => {
                let mut wtr =
                    csv::Writer::from_path(fmtr.filename.clone()).expect("could not create file");
                wtr.serialize(&fmtr.headers)
                    .expect("could not write headers");
                info!("Saving truth to {}", fmtr.filename);
                Some(wtr)
            }
            None => None,
        };

        let mut prop_setup = Propagator::default(self.truth.sc_dyn);
        prop_setup.set_tolerance(self.truth.prop_tol);
        let mut truth_prop = prop_setup.with(self.truth.init_state);

        // let mut truth_prop = self.truth.propagator();
        truth_prop.set_step(10.0 * Unit::Second, true);

        // Set up the channels
        let (tx, rx) = channel();

        let mut initial_state = Some(truth_prop.state);

        // Generate the measurements
        info!("Generating measurements over {} ", prop_time);

        let start = Instant::now();
        info!("Initial state: {}", truth_prop.state);

        truth_prop.for_duration_with_channel(prop_time, tx)?;

        info!(
            "Final state (computed in {:.3} seconds):",
            (Instant::now() - start).as_secs_f64()
        );
        info!("\t\t{}", truth_prop.state);
        info!("\t\t{:x}", truth_prop.state);

        let mut sim_measurements = Vec::with_capacity(10000);
        let start = Instant::now();
        while let Ok(rx_state) = rx.try_recv() {
            if let Some(wtr) = &mut maybe_wtr {
                if let Some(first_state) = initial_state {
                    wtr.serialize(
                        self.truth
                            .formatter
                            .as_ref()
                            .unwrap()
                            .fmt(&first_state.orbit),
                    )
                    .expect("could not format state");
                    initial_state = None;
                }
                wtr.serialize(self.truth.formatter.as_ref().unwrap().fmt(&rx_state.orbit))
                    .expect("could not format state");
            }
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
        let mut prop_setup = Propagator::default(self.nav.sc_dyn);
        prop_setup.set_tolerance(self.nav.prop_tol);
        // Make sure to set the STM
        let init_state = self.nav.init_state.with_stm();
        let mut nav = prop_setup.with(init_state);

        nav.set_step(10.0 * Unit::Second, true);

        let kf = self.kf;
        let trig = StdEkfTrigger::new(self.ekf_msr_trigger, self.ekf_disable_time);
        let mut odp = ODProcess::ekf(nav, kf, self.stations.clone(), trig);

        odp.process_measurements(&sim_measurements)?;

        odp.iterate(&sim_measurements, IterationConf::default())?;

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

        if let Some(final_estimate) = &odp.estimates.last() {
            println!("Final estimate:\n{}", final_estimate);
        }

        Ok(())
    }
}
