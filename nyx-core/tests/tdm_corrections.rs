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

    let expected_range_km = (15000.0 * 1e-9 * SPEED_OF_LIGHT_KM_S) / 2.0;
    assert!(
        (range_km - expected_range_km).abs() < 1e-10,
        "Range mismatch: got {range_km}, expected {expected_range_km}"
    );
}

#[test]
fn test_tdm_correction_range_us_ms_s() {
    for (unit, val, corr, expected_1way) in [
        (
            "us",
            10000.0,
            5000.0,
            (15000.0 * 1e-6 * SPEED_OF_LIGHT_KM_S) / 2.0,
        ),
        ("ms", 10.0, 5.0, (15.0 * 1e-3 * SPEED_OF_LIGHT_KM_S) / 2.0),
        ("s", 2.0, 1.0, (3.0 * SPEED_OF_LIGHT_KM_S) / 2.0),
    ] {
        let tdm_content = format!(
            r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = {unit}
  CORRECTION_RANGE = {corr}
  CORRECTIONS_APPLIED = NO
META_STOP
DATA_START
  RANGE = 2023-02-22T19:18:27.160 {val}
DATA_STOP
"#
        );

        let path =
            PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("target/test_corr_{unit}.tdm"));
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        let mut file = File::create(&path).unwrap();
        file.write_all(tdm_content.as_bytes()).unwrap();

        let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
        let msr = arc.measurements.first().unwrap();
        let range_km = msr.data.get(&MeasurementType::Range).unwrap();

        assert!(
            (range_km - expected_1way).abs() < 1e-8,
            "Failed for unit {unit}: got {range_km}, expected {expected_1way}"
        );
    }
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
    let msr = arc.measurements.first().unwrap();
    let range_km = msr.data.get(&MeasurementType::Range).unwrap();
    assert_eq!(*range_km, 7.5);
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
    let msr = arc.measurements.first().unwrap();
    let range_km = msr.data.get(&MeasurementType::Range).unwrap();
    assert_eq!(*range_km, 7500.0);
}

#[test]
fn test_tdm_correction_range_invalid_units() {
    for unit in ["RU", "invalid_unit"] {
        let tdm_content = format!(
            r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = {unit}
  CORRECTION_RANGE = 5000.0
  CORRECTIONS_APPLIED = NO
META_STOP
DATA_START
  RANGE = 2023-02-22T19:18:27.160 10000.0
DATA_STOP
"#
        );

        let path = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join(format!("target/test_corr_err_{unit}.tdm"));
        std::fs::create_dir_all(path.parent().unwrap()).unwrap();
        let mut file = File::create(&path).unwrap();
        file.write_all(tdm_content.as_bytes()).unwrap();

        assert!(TrackingDataArc::from_tdm(&path, None).is_err());
    }
}

#[test]
fn test_tdm_correction_applied_yes_and_invalid() {
    // Test YES
    let tdm_yes = r#"CCSDS_TDM_VERS = 2.0
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
    File::create(&path)
        .unwrap()
        .write_all(tdm_yes.as_bytes())
        .unwrap();
    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    let range_km = *arc
        .measurements
        .first()
        .unwrap()
        .data
        .get(&MeasurementType::Range)
        .unwrap();
    let expected_no_corr = (10000.0 * 1e-9 * SPEED_OF_LIGHT_KM_S) / 2.0;
    assert_eq!(range_km, expected_no_corr);

    // Test invalid flag string (warns and defaults to false -> applies correction)
    let tdm_invalid_flag = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = km
  CORRECTION_RANGE = 5000.0
  CORRECTIONS_APPLIED = INVALID_FLAG
META_STOP
DATA_START
  RANGE = 2023-02-22T19:18:27.160 10000.0
DATA_STOP
"#;
    let path2 = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_corr_inv_flag.tdm");
    File::create(&path2)
        .unwrap()
        .write_all(tdm_invalid_flag.as_bytes())
        .unwrap();
    let arc2 = TrackingDataArc::from_tdm(&path2, None).unwrap();
    let range_km2 = *arc2
        .measurements
        .first()
        .unwrap()
        .data
        .get(&MeasurementType::Range)
        .unwrap();
    assert_eq!(range_km2, 7500.0);
}

#[test]
fn test_tdm_correction_other_types() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  CORRECTION_ANGLE_1 = 0.5
  CORRECTION_ANGLE_2 = -0.2
  CORRECTIONS_APPLIED = NO
META_STOP
DATA_START
  ANGLE_1 = 2023-02-22T19:18:27.160 45.0
  ANGLE_2 = 2023-02-22T19:18:27.160 30.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_corr_azel.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    File::create(&path)
        .unwrap()
        .write_all(tdm_content.as_bytes())
        .unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    let msr = arc.measurements.first().unwrap();
    let az = msr.data.get(&MeasurementType::Azimuth).unwrap();
    let el = msr.data.get(&MeasurementType::Elevation).unwrap();

    // Azimuth and Elevation are 1-way measurements, so unscaled by msr_divider
    assert_eq!(*az, 45.5);
    assert_eq!(*el, 29.8);
}

#[test]
fn test_tdm_correction_invalid_value_string() {
    let tdm_content = r#"CCSDS_TDM_VERS = 2.0
META_START
  TIME_SYSTEM = UTC
  PARTICIPANT_1 = DSS14
  PARTICIPANT_2 = MySpacecraft
  MODE = SEQUENTIAL
  PATH = 1,2,1
  RANGE_UNITS = km
  CORRECTION_RANGE = NOT_A_NUMBER
  CORRECTIONS_APPLIED = NO
META_STOP
DATA_START
  RANGE = 2023-02-22T19:18:27.160 10000.0
DATA_STOP
"#;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("target/test_corr_nan.tdm");
    std::fs::create_dir_all(path.parent().unwrap()).unwrap();
    File::create(&path)
        .unwrap()
        .write_all(tdm_content.as_bytes())
        .unwrap();

    let arc = TrackingDataArc::from_tdm(&path, None).unwrap();
    let msr = arc.measurements.first().unwrap();
    let range_km = msr.data.get(&MeasurementType::Range).unwrap();

    // Correction ignored due to parse error, so 10000 / 2 = 5000 km
    assert_eq!(*range_km, 5000.0);
}
