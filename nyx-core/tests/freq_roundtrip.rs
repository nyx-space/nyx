use hifitime::prelude::*;
use nyx_space::od::msr::{MeasurementType, TrackingDataArc};
use nyx_space::od::GroundStation;
use nyx_space::io::ExportCfg;
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;
use anise::prelude::Almanac;
use std::sync::Arc;

#[test]
fn test_frequency_roundtrip() {
    let mut gs = GroundStation::default();
    gs.name = "DSS14".to_string();
    let rc = nyx_space::od::RadioConfig::new(2.1e9, 240.0 / 221.0, Duration::ZERO);
    gs.radio_config = Some(rc);

    // Create a fake Doppler measurement
    let epoch = Epoch::from_gregorian_utc_at_midnight(2023, 2, 22);
    let mut arc = TrackingDataArc::default();
    let mut msr = nyx_space::od::msr::measurement::Measurement::new(gs.name.clone(), epoch);
    msr.push(MeasurementType::Doppler, 0.1); // 100 m/s
    arc.measurements.insert(epoch, msr);

    // Convert to frequency data
    let arc_freq = arc.clone().to_frequency_data(rc);

    assert!(arc_freq.unique_types().contains(&MeasurementType::ReceiveFrequency));
    assert!(arc_freq.unique_types().contains(&MeasurementType::TransmitFrequency));

    // Export to TDM
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_roundtrip.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();

    let mut cfg = ExportCfg::default();
    cfg.metadata = Some(std::collections::HashMap::from([
        ("TURNAROUND_NUMERATOR".to_string(), "240".to_string()),
        ("TURNAROUND_DENOMINATOR".to_string(), "221".to_string()),
    ]));
    arc_freq.to_tdm_file("MySpacecraft".to_string(), None, &path, cfg).unwrap();

    // Read back from TDM
    let tdm_str = std::fs::read_to_string(&path).unwrap();
    println!("Exported TDM:\n{}", tdm_str);
    let arc_back = TrackingDataArc::from_tdm(&path, None).unwrap();

    // Check that Doppler is reconstructed correctly
    let msr_back = arc_back.measurements.get(&epoch).unwrap();
    let doppler_back = msr_back.data.get(&MeasurementType::Doppler).unwrap();

    assert!((doppler_back - 0.1).abs() < 1e-8, "Doppler mismatch: got {}, expected 0.1", doppler_back);

    // Check that original frequencies are also preserved
    assert!(msr_back.data.contains_key(&MeasurementType::ReceiveFrequency));
    assert!(msr_back.data.contains_key(&MeasurementType::TransmitFrequency));
}
