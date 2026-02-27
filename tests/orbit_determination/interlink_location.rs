extern crate nyx_space as nyx;

use anise::constants::frames::IAU_MOON_FRAME;
use anise::astro::{Location, TerrainMask};
use indexmap::IndexSet;
use nyx::cosmic::Orbit;
use nyx::md::prelude::*;
use nyx::od::interlink::InterlinkTxSpacecraft;
use nyx::od::prelude::*;
use nyx::time::Epoch;
use nyx::Spacecraft;
use nyx::od::GroundStation;

use anise::prelude::Almanac;
use rstest::*;
use std::sync::Arc;

use crate::test_almanac_arcd;

#[fixture]
fn almanac() -> Arc<Almanac> {
    test_almanac_arcd()
}

#[rstest]
fn test_interlink_elevation_mask(almanac: Arc<Almanac>) {
    let moon_iau = almanac.frame_info(IAU_MOON_FRAME).unwrap();

    let epoch = Epoch::from_gregorian_tai(2021, 5, 29, 19, 51, 16, 852_000);

    // 1. Define the Location (Rover on Moon) with a 5 degree mask
    let location = Location {
        latitude_deg: 0.0,
        longitude_deg: 0.0,
        height_km: 0.0,
        frame: moon_iau.into(),
        terrain_mask: TerrainMask::from_flat_terrain(5.0),
        terrain_mask_ignored: false,
    };

    // 2. Setup GroundStation (at same location)
    let mut gs = GroundStation {
        name: "RoverGS".to_string(),
        location: location.clone(),
        measurement_types: { 
            let mut set = IndexSet::new();
            set.insert(MeasurementType::Range);
            set
        },
        integration_time: None,
        light_time_correction: false,
        timestamp_noise_s: None,
        stochastic_noises: None,
    };

    // 3. Define the Rover Spacecraft at that location
    let rover_orbit = gs.to_orbit(epoch, &almanac).unwrap();
    let rover_sc = Spacecraft::from(rover_orbit);
    let mut rover_traj = Trajectory::new();
    rover_traj.states.push(rover_sc);

    // 4. Setup InterlinkTxSpacecraft with same Location
    let mut interlink = InterlinkTxSpacecraft {
        traj: rover_traj,
        location: Some(location),
        measurement_types: gs.measurement_types.clone(),
        integration_time: None,
        timestamp_noise_s: None,
        stochastic_noises: None,
        ab_corr: None,
    };

    // 5. Test at low elevation (e.g., 2 degrees)
    // 1000 km range, but very small elevation
    let orbiter_low = Orbit::cartesian(
        1737.4 + 10.0, 1000.0, 0.0, // Almost tangent to the surface
        0.0, 0.0, 0.0,
        epoch,
        moon_iau,
    );
    let orbiter_sc_low = Spacecraft::from(orbiter_low);

    let msr_interlink_low = interlink.measure_instantaneous(orbiter_sc_low, None, almanac.clone()).unwrap();
    let msr_gs_low = gs.measure_instantaneous(orbiter_sc_low, None, almanac.clone()).unwrap();

    println!("Interlink measurement at 2 deg: {:?}", msr_interlink_low.is_some());
    println!("GroundStation measurement at 2 deg: {:?}", msr_gs_low.is_some());

    assert!(msr_interlink_low.is_none(), "Interlink should have masked the measurement");
    assert!(msr_gs_low.is_none(), "GroundStation should have masked the measurement");

    // 6. Test at high elevation (e.g., Zenith)
    let orbiter_high = Orbit::cartesian(
        1737.4 + 1000.0, 0.0, 0.0,
        0.0, 0.0, 0.0,
        epoch,
        moon_iau,
    );
    let orbiter_sc_high = Spacecraft::from(orbiter_high);

    let msr_interlink_high = interlink.measure_instantaneous(orbiter_sc_high, None, almanac.clone()).unwrap();
    let msr_gs_high = gs.measure_instantaneous(orbiter_sc_high, None, almanac.clone()).unwrap();

    println!("Interlink measurement at Zenith: {:?}", msr_interlink_high.is_some());
    println!("GroundStation measurement at Zenith: {:?}", msr_gs_high.is_some());

    assert!(msr_interlink_high.is_some());
    assert!(msr_gs_high.is_some());
}
