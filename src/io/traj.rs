/*
    Nyx, blazing fast astrodynamics
    Copyright (C) 2021 Christopher Rabotin <christopher.rabotin@gmail.com>

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
extern crate crossbeam;
extern crate csv;

use self::crossbeam::thread;
use crate::io::scenario::StateSerde;
use crate::md::Ephemeris;
use crate::{celestia::Cosm, NyxError};
use std::sync::{mpsc::channel, Arc};

/// Imports a trajectory from a CSV. This is experimental and is subject to change between minor versions.
/// The CSV must have a column for each of the following fields with that specific case
/// [frame, epoch, x, y, z, vx, vy, vz]
/// The Epoch should be formatted in the ISO standard and include the time system (UTC, TAI, TT, etc.)
/// You may also specify `unit_position` and `unit_velocity` in an extra column for every item
/// **Caveats:** only the frame of the first item is considered. Also, this is super slow!
pub fn traj_from_csv(filename: &str, cosm: Arc<Cosm>) -> Result<Ephemeris, NyxError> {
    // Build the CSV reader and iterate over each record.
    match csv::Reader::from_path(filename) {
        Ok(mut rdr) => {
            thread::scope(|s| {
                let (tx, rx) = channel();
                let mut first = true;
                let traj_thread;
                let this_frame;

                let first_result: Result<StateSerde, _> = rdr.deserialize().next().unwrap();
                match first_result {
                    Ok(record) => {
                        // Compute the frame
                        this_frame = cosm.try_frame(&record.frame)?;
                        match record.as_state(this_frame) {
                            Ok(state) => {
                                traj_thread =
                                    s.spawn(move |_| Ephemeris::new_bucket_as(state, 8, rx));
                            }
                            Err(_) => {
                                return Err(NyxError::LoadingError(
                                    "Could not parse state".to_string(),
                                ))
                            }
                        };
                    }
                    Err(e) => return Err(NyxError::LoadingError(e.to_string())),
                }

                for result in rdr.deserialize() {
                    let maybe_result: Result<StateSerde, _> = result;
                    match maybe_result {
                        Ok(record) => {
                            match record.as_state(this_frame) {
                                Ok(state) => {
                                    if first {
                                        first = false;
                                        continue; // Skip the first state since we passed it to the initialization
                                    }
                                    tx.send(state).unwrap()
                                }
                                Err(_) => {
                                    return Err(NyxError::LoadingError(
                                        "Could not parse state".to_string(),
                                    ))
                                }
                            };
                        }
                        Err(e) => return Err(NyxError::LoadingError(e.to_string())),
                    }
                }

                let traj = traj_thread.join().unwrap_or_else(|_| {
                    Err(NyxError::NoInterpolationData(
                        "Could not generate trajectory".to_string(),
                    ))
                })?;
                Ok(traj)
            })
            .unwrap()
        }
        Err(e) => Err(NyxError::LoadingError(e.to_string())),
    }
}
