use anise::constants::SPEED_OF_LIGHT_KM_S;
use nyx_space::od::msr::{MeasurementType, TrackingDataArc};
use std::fs::File;
use std::io::Write;
use std::path::PathBuf;

#[test]
fn test_tdm_correction_range_ns() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = ns
  CORRECTION_RANGE = 5000.0
  CORRECTIONS_APPLIED = NO
META_STOP
DATA_START
  RANGE = 2023-02-22T19:18:27.160 10000.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_corr_ns.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 1);
    let msr = arc.measurements.first().unwrap();

    let range_km = msr
        .data
        .get(&MeasurementType::Range)
        .expect("Range measurement missing");

    // Corrected total round trip range in ns = 10000 + 5000 = 15000 ns
    // 1-way range in km = (15000 * 1e-9 * SPEED_OF_LIGHT_KM_S) / 2
    let expected_range_km = (15000.0 * 1e-9 * SPEED_OF_LIGHT_KM_S) / 2.0;
    assert!(
        (range_km - expected_range_km).abs() < 1e-10,
        "Range mismatch: got {range_km}, expected {expected_range_km}"
    );
}

#[test]
fn test_tdm_correction_range_m() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = m
  CORRECTION_RANGE = 5000.0
  CORRECTIONS_APPLIED = NO
META_STOP
DATA_START
  RANGE = 2023-02-22T19:18:27.160 10000.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_corr_m.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 1);
    let msr = arc.measurements.first().unwrap();

    let range_km = msr
        .data
        .get(&MeasurementType::Range)
        .expect("Range measurement missing");

    // Corrected total round trip range in m = 10000 + 5000 = 15000 m = 15 km
    // 1-way range in km = 15 / 2 = 7.5 km
    let expected_range_km = 7.5;
    assert!(
        (range_km - expected_range_km).abs() < 1e-10,
        "Range mismatch: got {range_km}, expected {expected_range_km}"
    );
}

#[test]
fn test_tdm_correction_range_km() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = km
  CORRECTION_RANGE = 5000.0
  CORRECTIONS_APPLIED = NO
META_STOP
DATA_START
  RANGE = 2023-02-22T19:18:27.160 10000.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_corr_km.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 1);
    let msr = arc.measurements.first().unwrap();

    let range_km = msr
        .data
        .get(&MeasurementType::Range)
        .expect("Range measurement missing");

    // Corrected total round trip range in km = 10000 + 5000 = 15000 km
    // 1-way range in km = 15000 / 2 = 7500 km
    let expected_range_km = 7500.0;
    assert!(
        (range_km - expected_range_km).abs() < 1e-10,
        "Range mismatch: got {range_km}, expected {expected_range_km}"
    );
}

#[test]
fn test_tdm_correction_applied_yes() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = ns
  CORRECTION_RANGE = 5000.0
  CORRECTIONS_APPLIED = YES
META_STOP
DATA_START
  RANGE = 2023-02-22T19:18:27.160 10000.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_corr_yes.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 1);
    let msr = arc.measurements.first().unwrap();

    let range_km = msr
        .data
        .get(&MeasurementType::Range)
        .expect("Range measurement missing");

    // Since CORRECTIONS_APPLIED = YES, correction is NOT added again.
    // 1-way range in km = (10000 * 1e-9 * SPEED_OF_LIGHT_KM_S) / 2
    let expected_range_km = (10000.0 * 1e-9 * SPEED_OF_LIGHT_KM_S) / 2.0;
    assert!(
        (range_km - expected_range_km).abs() < 1e-10,
        "Range mismatch: got {range_km}, expected {expected_range_km}"
    );
}

#[test]
fn test_tdm_correction_applied_default_omitted() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = ns
  CORRECTION_RANGE = 5000.0
META_STOP
DATA_START
  RANGE = 2023-02-22T19:18:27.160 10000.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_corr_omitted.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    let mut file = File::create(&path).unwrap();
    file.write_all(tdm_content.as_bytes()).unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    assert_eq!(arc.len(), 1);
    let msr = arc.measurements.first().unwrap();

    let range_km = msr
        .data
        .get(&MeasurementType::Range)
        .expect("Range measurement missing");

    // CCSDS TDM spec: CORRECTIONS_APPLIED defaults to NO when omitted.
    // So correction IS applied.
    let expected_range_km = (15000.0 * 1e-9 * SPEED_OF_LIGHT_KM_S) / 2.0;
    assert!(
        (range_km - expected_range_km).abs() < 1e-10,
        "Range mismatch: got {range_km}, expected {expected_range_km}"
    );
}
