//! Regression tests for the reference-frame contract of trajectory Parquet files.
//! A Parquet trajectory currently has one frame per table, not one per state.

use anise::constants::frames::{EARTH_J2000, IAU_EARTH_FRAME};
use anise::prelude::{Frame, Orbit};
use nyx_space::io::ExportCfg;
use nyx_space::md::prelude::Traj;
use nyx_space::time::{Epoch, TimeUnits};
use nyx_space::Spacecraft;
use std::path::PathBuf;

fn spacecraft(epoch: Epoch, frame: Frame) -> Spacecraft {
    Spacecraft::from(Orbit::new(
        7000.0, 400.0, 100.0, 0.0, 7.6, 0.1, epoch, frame,
    ))
}

fn trajectory(second_frame: Frame) -> Traj<Spacecraft> {
    let epoch = Epoch::from_gregorian_utc_at_noon(2026, 11, 11);
    Traj {
        name: Some("parquet-frame-fidelity".to_owned()),
        states: vec![
            spacecraft(epoch, EARTH_J2000),
            spacecraft(epoch + 30.seconds(), second_frame),
        ],
    }
}

fn temporary_path(suffix: &str) -> PathBuf {
    std::env::temp_dir().join(format!(
        "nyx-parquet-frame-{}-{suffix}.parquet",
        std::process::id()
    ))
}

fn round_trip(source: &Traj<Spacecraft>, suffix: &str) -> Result<Traj<Spacecraft>, String> {
    let path = temporary_path(suffix);
    let result = source
        .to_parquet(&path, ExportCfg::default())
        .map_err(|error| error.to_string())
        .and_then(|exported| {
            Traj::<Spacecraft>::from_parquet(&exported).map_err(|error| error.to_string())
        });
    let _ = std::fs::remove_file(&path);
    result
}

#[test]
fn uniform_frame_round_trip_preserves_frames_and_coordinates() {
    let source = trajectory(EARTH_J2000);
    let recovered = round_trip(&source, "uniform").expect("uniform-frame round-trip");
    assert_eq!(recovered.states.len(), 2);
    for (before, after) in source.states.iter().zip(recovered.states.iter()) {
        assert_eq!(before.orbit.frame, after.orbit.frame);
        assert_eq!(before.orbit.radius_km, after.orbit.radius_km);
        assert_eq!(before.orbit.velocity_km_s, after.orbit.velocity_km_s);
    }
}

#[test]
fn mixed_frames_must_be_preserved_or_explicitly_rejected() {
    let mixed = trajectory(IAU_EARTH_FRAME);
    let uniform = trajectory(EARTH_J2000);
    assert_ne!(mixed.states[1].orbit.frame, uniform.states[1].orbit.frame);

    let uniform_imported = round_trip(&uniform, "uniform-comparison").unwrap();
    let mixed_imported = match round_trip(&mixed, "mixed-comparison") {
        Ok(recovered) => recovered,
        Err(error) => {
            assert!(
                error.contains("Cannot export mixed-frame trajectory to Parquet"),
                "mixed-frame rejection must be explicit, got: {error}"
            );
            return;
        }
    };

    assert_eq!(mixed_imported.states.len(), 2);
    assert_ne!(
        uniform_imported.states[1].orbit.frame,
        mixed_imported.states[1].orbit.frame,
        "distinct reference frames silently collapsed on import"
    );
    assert_eq!(mixed_imported.states[1].orbit.frame, mixed.states[1].orbit.frame);
}

#[test]
fn filtered_export_uses_the_frame_of_the_selected_rows() {
    let source = trajectory(IAU_EARTH_FRAME);
    let selected_epoch = source.states[1].orbit.epoch;
    let cfg = ExportCfg::builder()
        .start_epoch(selected_epoch)
        .end_epoch(selected_epoch)
        .step(30.seconds())
        .build();
    let path = temporary_path("filtered");
    let exported = source.to_parquet(&path, cfg).expect("filtered export");
    let recovered = Traj::<Spacecraft>::from_parquet(&exported).expect("filtered import");
    let _ = std::fs::remove_file(&exported);

    assert_eq!(recovered.states.len(), 1);
    assert_eq!(
        recovered.states[0].orbit.frame,
        source.states[1].orbit.frame,
        "Frame metadata must describe the exported row, not an excluded row"
    );
}
